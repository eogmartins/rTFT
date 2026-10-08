# _precompute.R - fit the models the vignettes report, once
#
# The vignettes show the code that fits each model but do not run it when they are
# built: a fit takes minutes, and a built package cannot depend on a sampler being
# installed. This script fits every model the vignettes report and saves, for each
# vignette, only what it displays (tables of estimates and the tables the plots are
# drawn from), with the package version and the date. Everything else in a vignette
# runs when it is built.
#
# Each fit is checked before anything is saved: rank-normalized split-Rhat below
# 1.01, bulk and tail effective sample sizes above 400, and no divergent
# transitions. A fit that fails stops the script, so a vignette cannot report a
# model that did not converge.
#
# Rerun this script from the package root after a change to the fitting or deriving
# functions, then rebuild the vignettes and read them against the new numbers. With
# no argument it fits the models of every vignette. With the numbers of one or more
# vignettes it fits only theirs and leaves the saved results of the others as they
# are (e.g., `Rscript vignettes/_precompute.R 2` for the second vignette).
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
if (!length(vignettes)) vignettes <- 1:2
stopifnot("the vignettes are numbered 1 and 2" = all(vignettes %in% 1:2))

# Convergence, checked on every fit --------------------------------------------

# The convergence summary of the population-level parameters and the standard
# deviations of the group-level effects, after checking that the fit converged.
check_fit <- function(fit) {
  draws <- posterior::as_draws_array(fit)
  keep  <- grep("^(b|sd)_", posterior::variables(draws), value = TRUE)
  conv  <- posterior::summarise_draws(posterior::subset_draws(draws, variable = keep),
                                      "median", "sd", "rhat", "ess_bulk", "ess_tail")
  nuts  <- brms::nuts_params(fit)
  n_divergent <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  stopifnot("a fit did not converge" = all(conv$rhat < 1.01),
            "a fit has too few effective draws" = all(conv$ess_bulk > 400, conv$ess_tail > 400),
            "a fit has divergent transitions" = n_divergent == 0)
  list(table = as.data.frame(conv), n_divergent = n_divergent,
       n_draws = posterior::ndraws(draws))
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
  set.seed(seed)
  check_static <- check_survival(fit_static)
  set.seed(seed)
  residuals_static <- hazard_residuals(fit_static)

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
                       newdata = data.frame(mass = 15)))$mortality)

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
      wide = summarise_ramps(fit_wide)))
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
