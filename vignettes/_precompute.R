# _precompute.R - fit the models the vignettes report, once
#
# The vignettes show the code that fits each model but do not run it when they are
# built: a fit takes minutes, and a built package cannot depend on a sampler being
# installed. This script fits every model the vignettes report and saves, for each
# vignette, only what it displays (tables of estimates and the tables the plots are
# drawn from), with the package version and the date. Everything else in a vignette
# runs when it is built.
#
# Each fit is checked before anything is saved: rank-normalised split-Rhat below
# 1.01, bulk and tail effective sample sizes above 400, and divergent transitions
# in at most 0.1% of the draws (4 of 4,000). A fit that fails stops the script, so a
# vignette cannot report a model that did not converge. A few divergent transitions
# among thousands of draws, with Rhat and the effective sample sizes in order, do
# not by themselves indicate a biased posterior, and removing every one of them by
# shortening the sampler's steps can make a model several times slower to fit.
# Their number is saved with the results, and each vignette reports it.
#
# Rerun this script from the package root after a change to the fitting or deriving
# functions, then rebuild the vignettes and read them against the new numbers. With
# no argument it fits the models of every vignette. With the numbers of one or more
# vignettes it fits only theirs and leaves the saved results of the others as they
# are (e.g., `Rscript vignettes/_precompute.R 2` for the second vignette). The third
# vignette takes much longer than the others (about six hours), because each of its
# four models is fitted and then refitted 12 times, once without each trial.
#
# Inputs:
#   data/static_assay.rda
#   data/stream_record.rda
#   data/ramp_assay.rda
#   vignettes/results-01.rds (the second vignette compares its results with it)
#
# Outputs:
#   vignettes/results-01.rds
#   vignettes/results-02.rds
#   vignettes/results-03.rds
#   vignettes/results-04.rds

# Setup ------------------------------------------------------------------------

library(posterior)
library(pkgload)

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from vignettes/", call. = FALSE)
}
load_all(quiet = TRUE)

stamp <- list(version = as.character(utils::packageVersion("rTFT")),
              date    = format(Sys.Date()))
seed  <- 20261006

# The vignettes whose models are fitted in this run: all of them, unless their
# numbers were given on the command line.
vignettes <- as.integer(commandArgs(trailingOnly = TRUE))
if (!length(vignettes)) vignettes <- 1:4
stopifnot("the vignettes are numbered 1 to 4" = all(vignettes %in% 1:4))

# Convergence, checked on every fit --------------------------------------------

# The convergence summary of the population-level parameters and the standard
# deviations of the group-level effects, after checking that the fit converged.
max_divergent <- 0.001      # the largest share of draws from divergent transitions
check_fit <- function(fit) {
  draws <- posterior::as_draws_array(fit)
  keep  <- grep("^(b|sd)_", posterior::variables(draws), value = TRUE)
  conv  <- posterior::summarise_draws(posterior::subset_draws(draws, variable = keep),
                                      "median", "sd", "rhat", "ess_bulk", "ess_tail")
  nuts  <- brms::nuts_params(fit)
  n_divergent <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  n_draws     <- posterior::ndraws(draws)
  stopifnot("a fit did not converge" = all(conv$rhat < 1.01),
            "a fit has too few effective draws" = all(conv$ess_bulk > 400, conv$ess_tail > 400),
            "a fit has too many divergent transitions" =
              n_divergent <= max_divergent * n_draws)
  if (n_divergent > 0)
    message(n_divergent, " of ", n_draws, " draws came from divergent transitions.")
  list(table = as.data.frame(conv), n_divergent = n_divergent, n_draws = n_draws)
}

# Warm periods of the stream record ----------------------------------------------

# Predictions are often wanted for a warm period and not for a whole record, so each
# vignette also predicts the mortality over the five days of the heat wave of the
# stream record and over its warmest day. A part of the record is taken with its
# time counted from its first reading, so that exposure starts with it.
warm_periods <- list(heat_wave = 16:20, warmest_day = 18)
stopifnot("day 18 is not the warmest day of the record" =
            stream_record$day[which.max(stream_record$temp)] == 18)
part_of_record <- function(days, record = stream_record) {
  part <- record[record$day %in% days, ]
  part$time <- part$time - part$time[1]
  part
}
# The mortality of a 15 g fish over each warm period, predicted by a fitted model.
# `record` and `levels` are for a model with a co-stressor, whose record may follow
# the co-stressor and whose levels are otherwise given.
warm_period_mortality <- function(fit, record = stream_record, levels = NULL) {
  lapply(warm_periods, \(days) suppressMessages(
    predict_survival(fit, part_of_record(days, record), newdata = data.frame(mass = 15),
                     levels = levels))$mortality)
}

# Vignette 1: getting started, on the static assay ------------------------------

if (1 %in% vignettes) {
  static_path <- data.frame(
    path_id = rep(c("T23", "T25", "T27"), each = 2),
    time    = rep(c(0, 7200), times = 3),
    temp    = rep(c(23, 25, 27), each = 2))
  static_checks <- expand.grid(
    path_id = c("T23", "T25", "T27"),
    time    = seq(0, 7200, by = 60))

  long_static <- build_hazard_data(static_assay, static_path, "interval",
                                   covariates = c("mass", "trial"),
                                   checks = static_checks)

  fit_static <- fit_hazard(
    long_static, m_form = "exp", baseline = "weibull",
    dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
    # a smaller step size than the default, as an example of an option passed on to
    # the sampler
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr", refresh = 0)

  # The same model with the memoryless baseline, for the comparison of the two
  # predictions along the stream record.
  fit_static_mem <- fit_hazard(
    long_static, m_form = "exp", baseline = "mem",
    dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr", refresh = 0)
  invisible(check_fit(fit_static_mem))

  # The checks and the prediction draw random numbers, so the seed is set before each.
  # The residuals are drawn 20 times, so the vignette can show how much what they
  # show owes to one draw; under the same seed the first draw is the one a single
  # draw gives.
  set.seed(seed)
  check_static <- check_survival(fit_static)
  set.seed(seed)
  residuals_static <- hazard_residuals(fit_static, n_randomisations = 20)

  results_01 <- list(
    stamp       = stamp,
    convergence = check_fit(fit_static),
    # the posterior probability that the hazard rises with exposure time (k > 1)
    p_shape_above_1 = mean(as.matrix(fit_static, variable = "b_logshape_Intercept") > 0),
    check       = check_static,
    residuals   = residuals_static,
    # thermal tolerance at two body masses, on the default grids
    tolerance   = suppressMessages(
      derive_tdt_ttl(fit_static, newdata = data.frame(mass = c(10, 20)), t_ref = 60)),
    # the CTmax a ramping trial would record, for three ramps from 14 degC
    dctmax      = suppressMessages(
      derive_dctmax(fit_static, rate = c(0.1, 0.2, 0.3), T_start = 14,
                    newdata = data.frame(mass = 15))),
    # survival of a 15 g fish along the stream record, with each day's share of the risk
    prediction  = suppressMessages(
      predict_survival(fit_static, stream_record, newdata = data.frame(mass = 15),
                       by = "day")),
    # the same prediction under the memoryless baseline: its mortality only
    mortality_mem = suppressMessages(
      predict_survival(fit_static_mem, stream_record,
                       newdata = data.frame(mass = 15)))$mortality,
    # the mortality over the heat wave and over its warmest day, under each baseline
    warm_periods = list(weibull    = warm_period_mortality(fit_static),
                        memoryless = warm_period_mortality(fit_static_mem)))

  saveRDS(results_01, "vignettes/results-01.rds", compress = "xz")
}

# Vignette 2: ramping assays ---------------------------------------------------

if (2 %in% vignettes) {
  # The trials run in normoxia; those run at low dissolved oxygen are left to the
  # vignette on co-stressors. Two sets of three ramp rates are fitted: the narrow
  # one (0.1, 0.2 and 0.3 degC per min, a 3-fold span), which the vignette works
  # through, and the wide one (0.025, 0.1 and 0.3, a 12-fold span), with which it
  # shows what the span of the rates allows the model to estimate.
  ramp_normoxia <- ramp_assay[ramp_assay$do == 10, ]
  ramp_sets     <- list(narrow = c(0.1, 0.2, 0.3), wide = c(0.025, 0.1, 0.3))

  # No checks and no time step are given, so the intervals at risk follow each ramp
  # in steps of 0.1 degC, the default.
  fit_ramps <- function(rates) {
    paths <- ramp_path(sprintf("r%s", rates), T_start = 14, rate = rates,
                       t_max = 15.5 / rates)
    long  <- build_hazard_data(ramp_normoxia[ramp_normoxia$rate %in% rates, ], paths,
                               "continuous", covariates = c("mass", "trial"))
    fit_hazard(
      long, m_form = "exp", baseline = "weibull",
      dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
      chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
      control = list(adapt_delta = 0.95),
      backend = "cmdstanr", refresh = 0)
  }

  # What the vignette reports of both fits: convergence, the parameters and sCTmax
  # of a 15 g fish, its dCTmax under the three ramps of the narrow set, and its
  # survival along the stream record.
  summarise_ramps <- function(fit) {
    tolerance <- suppressMessages(
      derive_tdt_ttl(fit, newdata = data.frame(mass = 15), t_ref = 60,
                     # the temperatures at which fish failed in the narrow set (the
                     # default would start at 14 degC, where every ramp began)
                     Ta = seq(25.5, 29.5, by = 0.1)))
    list(
      convergence = check_fit(fit),
      tolerance   = tolerance,
      dctmax      = suppressMessages(
        derive_dctmax(fit, rate = ramp_sets$narrow, T_start = 14,
                      newdata = data.frame(mass = 15))),
      prediction  = suppressMessages(
        predict_survival(fit, stream_record, newdata = data.frame(mass = 15))))
  }

  fit_narrow <- fit_ramps(ramp_sets$narrow)
  fit_wide   <- fit_ramps(ramp_sets$wide)

  set.seed(seed)
  check_ramp <- check_survival(fit_narrow)
  set.seed(seed)
  residuals_ramp <- hazard_residuals(fit_narrow)
  # The same checks of the wide set, which the vignette reports without showing.
  set.seed(seed)
  check_wide <- check_survival(fit_wide)
  set.seed(seed)
  residuals_wide <- hazard_residuals(fit_wide)

  results_02 <- c(
    list(stamp = stamp),
    summarise_ramps(fit_narrow),
    list(
      # the posterior probability that the hazard rises with exposure time (k > 1)
      p_shape_above_1 = mean(as.matrix(fit_narrow, variable = "b_logshape_Intercept") > 0),
      check     = check_ramp,
      residuals = residuals_ramp,
      # the dCTmax under ramps from a tenth of the slowest rate of the narrow set
      # to about three times its fastest
      dctmax_rates = suppressMessages(
        derive_dctmax(fit_narrow, rate = seq(0.01, 1, by = 0.01),
                      T_start = 14, newdata = data.frame(mass = 15))),
      # the same summaries for the wide set
      wide = summarise_ramps(fit_wide),
      # the mortality over the heat wave and over its warmest day, from each set
      warm_periods = list(narrow = warm_period_mortality(fit_narrow),
                          wide   = warm_period_mortality(fit_wide))))
  # The thermal death time curves and the landscape of the wide set are not shown.
  results_02$wide$tolerance <- results_02$wide$tolerance[c("params", "sctmax")]
  results_02$wide$check     <- check_wide
  results_02$wide$residuals <- residuals_wide

  # Two statements the vignette makes about these results without quoting a number:
  # the observed survival lies within the simulated band throughout, and the dCTmax
  # derived here and from the static assay of the first vignette have overlapping
  # intervals at the three rates.
  dctmax_static <- readRDS("vignettes/results-01.rds")$dctmax
  stopifnot(
    "an observed proportion surviving lies outside the simulated band" =
      all(check_ramp$observed >= check_ramp$lower,
          check_ramp$observed <= check_ramp$upper),
    "an observed proportion surviving lies outside the simulated band, wide set" =
      all(check_wide$observed >= check_wide$lower,
          check_wide$observed <= check_wide$upper),
    # the residuals of each fit are compatible with a uniform distribution
    "the residuals of the narrow set depart from the uniform distribution" =
      stats::ks.test(residuals_ramp$residual, "punif")$p.value > 0.05,
    "the residuals of the wide set depart from the uniform distribution" =
      stats::ks.test(residuals_wide$residual, "punif")$p.value > 0.05,
    "the two assays were not derived for the same ramps" =
      isTRUE(all.equal(dctmax_static$rate, results_02$dctmax$rate)),
    "the dCTmax intervals of the ramping and static assays do not overlap" =
      all(results_02$dctmax$lower <= dctmax_static$upper,
          dctmax_static$lower <= results_02$dctmax$upper))

  saveRDS(results_02, "vignettes/results-02.rds", compress = "xz")
}

# Vignette 3: comparing and averaging models -----------------------------------

if (3 %in% vignettes) {
  # The ramping trials in normoxia at 0.025, 0.1 and 0.3 degC per min (the set of
  # rates with the wider span in the second vignette), fitted with four models: the
  # exponential failure rate with the memoryless, Weibull and log-logistic
  # baselines, and the zero-power threshold failure rate with the memoryless one.
  rates_3 <- c(0.025, 0.1, 0.3)
  ramp_3  <- ramp_assay[ramp_assay$do == 10 & ramp_assay$rate %in% rates_3, ]
  paths_3 <- ramp_path(sprintf("r%s", rates_3), T_start = 14, rate = rates_3,
                       t_max = 15.5 / rates_3)
  long_3  <- build_hazard_data(ramp_3, paths_3, "continuous",
                               covariates = c("mass", "trial"))

  candidates <- list(
    "Exponential, memoryless"          = c("exp", "mem"),
    "Exponential, Weibull"             = c("exp", "weibull"),
    "Exponential, log-logistic"        = c("exp", "llogis"),
    "Zero-power threshold, memoryless" = c("zeropower_threshold", "mem"))
  fits_3 <- lapply(candidates, \(model) fit_hazard(
    long_3, m_form = model[1], baseline = model[2],
    dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr", refresh = 0))
  convergence_3 <- lapply(fits_3, check_fit)

  # Each model against the data, by the two checks of the first vignette: the number
  # of times at which the proportion observed to survive lies outside the simulated
  # band, in each ramp, and whether the residuals depart from the uniform
  # distribution (Kolmogorov-Smirnov test).
  fit_checks <- do.call(rbind, lapply(names(fits_3), \(m) {
    set.seed(seed)
    check <- check_survival(fits_3[[m]])
    set.seed(seed)
    residuals <- hazard_residuals(fits_3[[m]])
    outside <- !(check$observed >= check$lower & check$observed <= check$upper)
    data.frame(model = m, n_times = nrow(check), n_outside = sum(outside),
               t(tapply(outside, check$path_id, sum)), check.names = FALSE,
               ks_p = suppressWarnings(stats::ks.test(residuals$residual, "punif")$p.value))
  }))

  # Leave-one-trial-out cross-validation by importance sampling. Its warnings are
  # silenced here and the number of trials they concern is recorded, because the
  # vignette shows what the function says when the approximation cannot be trusted.
  psis <- lapply(fits_3, \(fit) suppressWarnings(suppressMessages(loo_hazard(fit))))
  psis_table <- data.frame(
    model     = names(psis),
    elpd      = vapply(psis, \(x) x$estimates["elpd_loo", "Estimate"], numeric(1)),
    se        = vapply(psis, \(x) x$estimates["elpd_loo", "SE"], numeric(1)),
    max_k     = vapply(psis, \(x) max(x$diagnostics$pareto_k), numeric(1)),
    n_high_k  = vapply(psis, \(x) sum(x$diagnostics$pareto_k > 0.7), numeric(1)),
    row.names = NULL)

  # The same by refitting: each model is refitted 12 times, once without each trial.
  cv_3 <- lapply(fits_3, \(fit) suppressMessages(
    kfold_hazard(fit, cores = 4, seed = seed, refresh = 0)))
  stopifnot("a trial was not left out on its own" =
              all(vapply(cv_3, \(x) x$K == 12 && all(table(x$folds) == 1), logical(1))))
  # The comparison as a table with the models in a `model` column, which is where
  # loo has them from version 2.10.0; before it they were the row names of a matrix.
  compare_3 <- as.data.frame(loo::loo_compare(cv_3))
  if (is.null(compare_3$model)) {
    compare_3 <- data.frame(model = rownames(compare_3), compare_3, row.names = NULL)
  }

  weights_3 <- suppressMessages(hazard_weights(cv_3))
  set.seed(seed)
  weights_pbma <- suppressMessages(hazard_weights(cv_3, method = "pseudobma"))

  # What each model gives on its own, and what their average gives, for a 15 g
  # fish: sCTmax for 60 min, the thermal death time curves over the temperatures at
  # which most fish failed (below about 25 degC the times to failure of the models
  # differ by many orders of magnitude), dCTmax at a rate of the assay and at two
  # outside it, and survival along the stream record.
  fish  <- data.frame(mass = 15)
  Ta_3  <- seq(25.5, 29.5, by = 0.1)
  each_tolerance <- lapply(fits_3, \(fit) suppressMessages(
    derive_tdt_ttl(fit, newdata = fish, t_ref = 60, Ta = Ta_3))[c("sctmax", "tdt_curves")])
  each_dctmax <- lapply(fits_3, \(fit) suppressMessages(
    derive_dctmax(fit, rate = c(0.01, 0.3, 1), T_start = 14, newdata = fish)))
  each_prediction <- lapply(fits_3, \(fit) suppressMessages(
    predict_survival(fit, stream_record, newdata = fish)))
  average_tolerance <- suppressMessages(
    average_tdt_ttl(fits_3, weights_3, newdata = fish, t_ref = 60, Ta = Ta_3))
  average_prediction <- suppressMessages(
    average_survival(fits_3, weights_3, stream_record, newdata = fish))

  # Further summaries of the same fits. The parameters of each model, on their
  # natural scale; the survival along the record averaged with the pseudo-BMA
  # weights, to set beside that averaged with the stacking weights; and the
  # mortality predicted for two shorter parts of the record, the five days of the
  # heat wave and its warmest day, by each model and averaged, to show how the
  # disagreement among the models grows with the length of the record.
  each_params <- lapply(fits_3, \(fit) suppressMessages(
    derive_tdt_ttl(fit, newdata = fish, t_ref = 60, Ta = 27))$params)
  average_prediction_pbma <- suppressMessages(
    average_survival(fits_3, weights_pbma, stream_record, newdata = fish))[
      c("mortality", "weights")]
  shorter_records <- lapply(warm_periods, \(days) {
    part <- part_of_record(days)
    list(hours = max(part$time) / 60,
         each = lapply(fits_3, \(fit) suppressMessages(
           predict_survival(fit, part, newdata = fish))$mortality),
         average = suppressMessages(
           average_survival(fits_3, weights_3, part, newdata = fish))$mortality)
  })

  # The dCTmax averaged over the models, under the ramps it is derived for from each
  # model, and the sCTmax averaged with the pseudo-BMA weights, to set beside that
  # averaged with the stacking weights.
  average_dctmax_3 <- suppressMessages(
    average_dctmax(fits_3, weights_3, rate = c(0.01, 0.3, 1), T_start = 14,
                   newdata = fish))$dctmax
  average_sctmax_pbma <- suppressMessages(
    average_tdt_ttl(fits_3, weights_pbma, newdata = fish, t_ref = 60, Ta = Ta_3))$sctmax

  results_03 <- list(
    stamp       = stamp,
    convergence = convergence_3,
    fit_checks  = fit_checks,
    psis        = psis_table,
    kfold       = data.frame(
      model = names(cv_3),
      elpd  = vapply(cv_3, \(x) x$estimates["elpd_kfold", "Estimate"], numeric(1)),
      se    = vapply(cv_3, \(x) x$estimates["elpd_kfold", "SE"], numeric(1)),
      row.names = NULL),
    folds       = cv_3[[1]]$folds,
    compare     = compare_3,
    weights     = weights_3,
    weights_pseudobma  = weights_pbma,
    each_tolerance     = each_tolerance,
    each_dctmax        = each_dctmax,
    each_prediction    = each_prediction,
    average_tolerance  = average_tolerance,
    average_prediction = average_prediction,
    each_params        = each_params,
    average_prediction_pseudobma = average_prediction_pbma,
    shorter_records    = shorter_records,
    average_dctmax     = average_dctmax_3,
    average_sctmax_pseudobma = average_sctmax_pbma)

  saveRDS(results_03, "vignettes/results-03.rds", compress = "xz")
}

# Vignette 4: temperature and a co-stressor ------------------------------------

if (4 %in% vignettes) {
  # The ramping trials at the three rates of the wide span (0.025, 0.1 and 0.3 degC
  # per min) at both levels of oxygen, with dissolved oxygen as a co-stressor. Its
  # reference is normoxia (10 mg/L), so the coefficient is per mg/L of oxygen deficit.
  # Two models are fitted: oxygen with its main effect alone, and with its interaction
  # with temperature as well.
  rates_4 <- c(0.025, 0.1, 0.3)
  ramp_4  <- ramp_assay[ramp_assay$rate %in% rates_4, ]
  stopifnot(nrow(ramp_4) == 240, all(table(ramp_4$rate, ramp_4$do) == 40))
  paths_4 <- ramp_path(sprintf("r%s", rates_4), T_start = 14, rate = rates_4,
                       t_max = 15.5 / rates_4)
  oxygen <- list(
    "Main effect" = costressor("DO", column = "do", ref = 10, scale = 1, direction = -1),
    "Main effect and interaction" = costressor("DO", column = "do", ref = 10, scale = 1,
                                               direction = -1, interaction = TRUE))
  # The table of intervals at risk is the same for both: it holds the oxygen
  # deficit of each row whichever way oxygen enters the model.
  long_4 <- build_hazard_data(ramp_4, paths_4, "continuous",
                              covariates = c("mass", "trial"),
                              costressors = list(oxygen[[1]]))

  fits_4 <- lapply(oxygen, \(DO) fit_hazard(
    long_4, m_form = "exp", baseline = "weibull",
    dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
    costressors = list(DO),
    chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr", refresh = 0))
  convergence_4 <- lapply(fits_4, check_fit)

  # Each model against the data, for each ramp and level of oxygen.
  fit_checks_4 <- do.call(rbind, lapply(names(fits_4), \(m) {
    set.seed(seed)
    check <- check_survival(fits_4[[m]], by = "z_DO")
    set.seed(seed)
    residuals <- hazard_residuals(fits_4[[m]])
    outside <- !(check$observed >= check$lower & check$observed <= check$upper)
    data.frame(model = m, n_times = nrow(check), n_outside = sum(outside),
               ks_p = suppressWarnings(stats::ks.test(residuals$residual, "punif")$p.value))
  }))

  fish <- data.frame(mass = 15)
  Ta_4 <- seq(25.5, 29.5, by = 0.1)
  do_levels <- c(10, 5)

  # The stream record with an oxygen series that falls as the stream warms: 10 mg/L
  # at 15 degC and below, and 0.5 mg/L less for each degree above it.
  record_4 <- stream_record
  record_4$do <- pmin(10, 10 - 0.5 * (record_4$temp - 15))

  # What the vignette reports of each model: its parameters; sCTmax, the thermal
  # death time curves, the thermal tolerance landscape and dCTmax at each level of
  # oxygen; the effect of 1 mg/L of
  # oxygen deficit on the log failure rate across temperatures; and mortality along
  # the stream record, over the 30 days and over the warm periods, with oxygen at
  # 10 mg/L throughout and with oxygen following the series above.
  summarise_oxygen <- function(fit) {
    tolerance <- lapply(do_levels, \(level) suppressMessages(
      derive_tdt_ttl(fit, newdata = fish, t_ref = 60, levels = c(DO = level), Ta = Ta_4)))
    dctmax <- lapply(do_levels, \(level) suppressMessages(
      derive_dctmax(fit, rate = c(0.1, 0.3), T_start = 14, newdata = fish,
                    levels = c(DO = level))))
    temperatures <- seq(16, 29, by = 0.5)
    rate_at <- function(level) {
      suppressMessages(M_draws(fit, Ta = temperatures, newdata = fish,
                               levels = c(DO = level), n_draws = 1000))
    }
    normoxia <- rate_at(10)
    deficit  <- rate_at(9)
    stopifnot(identical(normoxia[c(".draw", "Ta")], deficit[c(".draw", "Ta")]))
    effect <- log(deficit$M / normoxia$M)
    oxygen_effect <- do.call(rbind, lapply(temperatures, \(Ta) {
      x <- effect[normoxia$Ta == Ta]
      data.frame(Ta = Ta, median = stats::median(x),
                 lower = unname(stats::quantile(x, 0.025)),
                 upper = unname(stats::quantile(x, 0.975)), p_negative = mean(x < 0))
    }))
    whole <- function(record, levels) suppressMessages(
      predict_survival(fit, record, newdata = fish, levels = levels))
    list(
      params = tolerance[[1]]$params,
      sctmax = stats::setNames(lapply(tolerance, `[[`, "sctmax"), do_levels),
      tdt_curves = stats::setNames(lapply(tolerance, `[[`, "tdt_curves"), do_levels),
      ttl_surface = stats::setNames(lapply(tolerance, `[[`, "ttl_surface"), do_levels),
      dctmax = stats::setNames(dctmax, do_levels),
      oxygen_effect = oxygen_effect,
      prediction = list(normoxia = whole(stream_record, c(DO = 10)),
                        following = whole(record_4, NULL)),
      warm_periods = list(normoxia  = warm_period_mortality(fit, levels = c(DO = 10)),
                          following = warm_period_mortality(fit, record = record_4)))
  }

  results_04 <- list(
    stamp       = stamp,
    convergence = convergence_4,
    fit_checks  = fit_checks_4,
    models      = lapply(fits_4, summarise_oxygen),
    # the posterior probability that the interaction is positive
    p_interaction_positive =
      mean(as.matrix(fits_4[[2]], variable = "b_dTDO_Intercept") > 0),
    # the reference temperature and time of the fitted models
    Tr = unique(long_4$Tr), tau_ref = unique(long_4$tau_ref))

  # The interaction is estimated in the direction of the model the data come from,
  # and its credible interval holds the simulated value.
  truth_4    <- attr(ramp_assay, "truth")
  interaction_4 <- with(results_04$models[[2]], params[params$quantity == "dTDO", ])
  stopifnot(results_04$p_interaction_positive > 0.95,
            interaction_4$lower <= truth_4$d_do, truth_4$d_do <= interaction_4$upper)

  saveRDS(results_04, "vignettes/results-04.rds", compress = "xz")
}
