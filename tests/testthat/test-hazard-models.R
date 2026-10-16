# ===========================================================================
# 1. Helpers
# ===========================================================================

# Three fitted-model fixtures. `exp` and `zeroexp` are two models of the same 40
# individuals, with body mass on log a and no group-level term. `fit_trials` has
# eight trials of six individuals, with a random intercept per trial on log a. Each
# has two chains and 600 draws.
fits <- list(exp     = readRDS(test_path("fixtures", "hazard-fit-exp.rds")),
             zeroexp = readRDS(test_path("fixtures", "hazard-fit-zeroexp.rds")))
fit_trials <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))

# The table each was fitted to, as far as fit_hazard_models() reads it: the fitted
# data with the individual of each row.
long_of <- function(fit) cbind(fit$data, id = fit$rtft$id)
long        <- long_of(fits$exp)
long_trials <- long_of(fit_trials)

two_models <- build_hazard_models(c("exp", "zeroexp"),
                                  dpar_formulas = list(loga = ~ 1 + mass))
trial_models <- build_hazard_models("exp", c("mem", "weibull"),
                                    dpar_formulas = list(loga = ~ 1 + (1 | trial)))

# Fitting takes minutes, so these tests replace fit_hazard() with a stand-in that
# records what it was given and returns a stored fit: that of the failure-rate
# function asked for, or `fit_trials` where one is given. The power function has no
# stored fit and stands for a model whose fit fails.
stand_in_fit <- function(record, always = NULL) {
  function(long, m_form, baseline, dpar_formulas, observation, costressors, ...) {
    record$calls[[length(record$calls) + 1]] <- list(
      m_form = m_form, baseline = baseline, dpar_formulas = dpar_formulas,
      observation = observation, costressors = costressors, args = list(...))
    if (m_form == "power") stop("the sampler stopped")
    if (is.null(always)) fits[[m_form]] else always
  }
}
new_record <- function() {
  record <- new.env()
  record$calls <- list()
  record
}

# The stand-in for a refit of test-loo-hazard.R: it returns, for each unit left out,
# minus the unit's position among all the units, and fails from the call numbered
# `fail_from` on.
stand_in_refit <- function(record, fail_from = Inf) {
  function(fit, data, train, unit_values, ...) {
    record$calls[[length(record$calls) + 1]] <- list(args = list(...))
    if (length(record$calls) >= fail_from) stop("the refit stopped")
    held <- unique(as.character(unit_values[!train]))
    list(elpd = stats::setNames(-match(held, unique(as.character(unit_values))), held),
         fit = fit)
  }
}

# Run an expression without the package's messages and warnings, which the tests
# that are not about them do not need.
quietly <- function(expr) suppressWarnings(suppressMessages(expr))

# The number of conditions of a class that an expression gives.
count_conditions <- function(expr, class) {
  n <- 0
  suppressWarnings(suppressMessages(
    withCallingHandlers(expr, condition = function(cnd) {
      if (inherits(cnd, class)) n <<- n + 1
    })))
  n
}


# ===========================================================================
# 2. A table of models (build_hazard_models)
# ===========================================================================

test_that("build_hazard_models combines every failure-rate function with every baseline", {
  models <- build_hazard_models(c("exp", "zeropower_threshold"), c("mem", "weibull"))
  expect_s3_class(models, "data.frame")
  expect_named(models, c("model", "m_form", "baseline", "dpar_formulas", "costressors"))
  # each failure-rate function in turn, with every baseline
  expect_equal(models$m_form, rep(c("exp", "zeropower_threshold"), each = 2))
  expect_equal(models$baseline, rep(c("mem", "weibull"), 2))
  expect_equal(models$model, c("exp_mem", "exp_weibull", "zeropower_threshold_mem",
                               "zeropower_threshold_weibull"))
  # with no formula and no co-stressor, each model has an empty list of each
  expect_true(all(lengths(models$dpar_formulas) == 0))
  expect_true(all(lengths(models$costressors) == 0))
  # the baseline defaults to the memoryless one
  expect_equal(build_hazard_models("exp")$baseline, "mem")
})


test_that("a formula is given to the models that have its parameter, and to no other", {
  f_loga  <- ~ 1 + (1 | trial)
  f_shape <- ~ 1 + mass
  f_Tc    <- ~ 1 + acclim
  models <- build_hazard_models(c("exp", "zeropower_threshold"), c("mem", "weibull"),
                                dpar_formulas = list(loga = f_loga, logshape = f_shape,
                                                     Tc = f_Tc))
  given <- stats::setNames(lapply(models$dpar_formulas, names), models$model)
  expect_equal(given$exp_mem, "loga")
  expect_equal(given$exp_weibull, c("loga", "logshape"))
  expect_equal(given$zeropower_threshold_mem, c("loga", "Tc"))
  expect_equal(given$zeropower_threshold_weibull, c("loga", "logshape", "Tc"))
  expect_equal(models$dpar_formulas[[4]]$logshape, f_shape)
  # a parameter that none of the models has is an error, not a formula dropped
  expect_error(build_hazard_models("exp", dpar_formulas = list(logshape = f_shape)),
               class = "rTFT_formula_for_no_model")
  expect_error(build_hazard_models("exp", dpar_formulas = list(logA = f_loga)),
               class = "rTFT_formula_for_no_model")
})


test_that("with cross = FALSE the functions and baselines are paired, and models can be named", {
  models <- build_hazard_models(c("exp", "exp", "zeropower_threshold"),
                                c("mem", "llogis", "mem"), cross = FALSE,
                                model_names = c("first", "second", "third"))
  expect_equal(models$m_form, c("exp", "exp", "zeropower_threshold"))
  expect_equal(models$baseline, c("mem", "llogis", "mem"))
  expect_equal(models$model, c("first", "second", "third"))
  # a single value is repeated
  one_baseline <- build_hazard_models(c("exp", "zeroexp"), "weibull", cross = FALSE)
  expect_equal(one_baseline$baseline, c("weibull", "weibull"))
  expect_error(build_hazard_models(c("exp", "zeroexp", "const"), c("mem", "weibull"),
                                   cross = FALSE), class = "rTFT_bad_models")
})


test_that("a co-stressor is given to every model, and its coefficients can have formulas", {
  do <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                   interaction = TRUE)
  models <- build_hazard_models(c("exp", "zeroexp"), costressors = list(do),
                                dpar_formulas = list(cDO = ~ 1 + mass))
  expect_equal(models$costressors, list(list(do), list(do)))
  expect_equal(lapply(models$dpar_formulas, names), list("cDO", "cDO"))
  # without the co-stressor there is no such parameter
  expect_error(build_hazard_models("exp", dpar_formulas = list(cDO = ~ 1 + mass)),
               class = "rTFT_formula_for_no_model")
})


test_that("build_hazard_models refuses what is not a set of models", {
  expect_error(build_hazard_models("nope"), class = "rTFT_unknown_m_form")
  expect_error(build_hazard_models("exp", "nope"), class = "rTFT_unknown_baseline")
  expect_error(build_hazard_models(character(0)), class = "rTFT_bad_models")
  expect_error(build_hazard_models(1), class = "rTFT_bad_models")
  expect_error(build_hazard_models("exp", NA_character_), class = "rTFT_bad_models")
  expect_error(build_hazard_models("exp", cross = NA), class = "rTFT_bad_models")
  # the same model twice
  expect_error(build_hazard_models(c("exp", "exp")), class = "rTFT_bad_models")
  expect_error(build_hazard_models(c("exp", "exp"), c("mem", "mem"), cross = FALSE),
               class = "rTFT_bad_models")
  # formulas that are not a named list of formulas
  expect_error(build_hazard_models("exp", dpar_formulas = list(~ 1 + mass)),
               class = "rTFT_bad_models")
  expect_error(build_hazard_models("exp", dpar_formulas = list(loga = "1 + mass")),
               class = "rTFT_bad_models")
  expect_error(build_hazard_models("exp", dpar_formulas = ~ 1 + mass),
               class = "rTFT_bad_models")
  # model names: one per model, each its own
  expect_error(build_hazard_models(c("exp", "zeroexp"), model_names = "one"),
               class = "rTFT_bad_models")
  expect_error(build_hazard_models(c("exp", "zeroexp"), model_names = c("a", "a")),
               class = "rTFT_bad_models")
  expect_error(build_hazard_models(c("exp", "zeroexp"), model_names = c("a", "")),
               class = "rTFT_bad_models")
})


test_that(".check_models accepts a table that was bound or cut, and refuses a broken one", {
  bound <- rbind(two_models, build_hazard_models("exp", "weibull"))
  expect_no_error(.check_models(bound))
  expect_no_error(.check_models(two_models[2, ]))
  expect_error(.check_models(rbind(two_models, two_models)), class = "rTFT_bad_models")
  expect_error(.check_models(two_models[0, ]), class = "rTFT_bad_models")
  expect_error(.check_models(two_models[c("model", "m_form")]), class = "rTFT_bad_models")
  expect_error(.check_models(list(model = "a")), class = "rTFT_bad_models")
  no_lists <- as.data.frame(two_models[c("model", "m_form", "baseline")])
  no_lists$dpar_formulas <- "~ 1"
  no_lists$costressors   <- "none"
  expect_error(.check_models(no_lists), class = "rTFT_bad_models")
  edited <- two_models
  edited$m_form[2] <- "nope"
  expect_error(.check_models(edited), class = "rTFT_unknown_m_form")
  edited <- two_models
  edited$baseline[1] <- "nope"
  expect_error(.check_models(edited), class = "rTFT_unknown_baseline")
  # a table read back with its names as factors is refused, not read by position
  as_factors <- as.data.frame(two_models)
  as_factors$m_form <- factor(as_factors$m_form)
  expect_error(.check_models(as_factors), class = "rTFT_bad_models")
  as_factors <- as.data.frame(two_models)
  as_factors$baseline <- factor(as_factors$baseline)
  expect_error(.check_models(as_factors), class = "rTFT_bad_models")
})


# ===========================================================================
# 3. The convergence rule
# ===========================================================================

test_that(".fit_convergence reads the diagnostics of the parameters, not of each group's effect", {
  conv <- .fit_convergence(fit_trials)
  expect_named(conv, c("chains", "draws", "rhat", "ess_bulk", "ess_tail", "divergent",
                       "rhat_group", "ess_group"))
  expect_equal(conv$chains, 2)
  expect_equal(conv$draws, 600)
  # by hand: the two intercepts and the standard deviation among trials
  draws <- posterior::subset_draws(
    posterior::as_draws_array(fit_trials),
    variable = c("b_loga_Intercept", "b_logb_Intercept", "sd_trial__loga_Intercept"))
  by_hand <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
  expect_equal(conv$rhat, max(by_hand$rhat))
  expect_equal(conv$ess_bulk, min(by_hand$ess_bulk))
  expect_equal(conv$ess_tail, min(by_hand$ess_tail))
  nuts <- brms::nuts_params(fit_trials)
  expect_equal(conv$divergent, sum(nuts$Value[nuts$Parameter == "divergent__"]))
  # the effect of each trial, read beside the rule and not in it
  effects <- grep("^r_trial", posterior::variables(posterior::as_draws_array(fit_trials)),
                  value = TRUE)
  expect_length(effects, 8)
  by_trial <- posterior::summarise_draws(
    posterior::subset_draws(posterior::as_draws_array(fit_trials), variable = effects),
    "rhat", "ess_bulk", "ess_tail")
  expect_equal(conv$rhat_group, max(by_trial$rhat))
  expect_equal(conv$ess_group, min(by_trial$ess_bulk, by_trial$ess_tail))
  # a model without group-level effects has none to report
  no_groups <- .fit_convergence(fits$exp)
  expect_true(is.na(no_groups$rhat_group))
  expect_true(is.na(no_groups$ess_group))
})


test_that("the rule needs every diagnostic within its limit, and none missing", {
  good <- list(chains = 4, draws = 4000, rhat = 1.004, ess_bulk = 900, ess_tail = 700,
               divergent = 4L)
  rule <- function(conv) .is_converged(conv, rhat = 1.01, ess = 400, divergent = 0.001)
  expect_true(rule(good))
  expect_false(rule(utils::modifyList(good, list(rhat = 1.01))))        # at the limit
  expect_false(rule(utils::modifyList(good, list(ess_bulk = 400))))
  expect_false(rule(utils::modifyList(good, list(ess_tail = 399))))
  expect_false(rule(utils::modifyList(good, list(divergent = 5L))))     # 5 of 4000
  expect_true(.is_converged(utils::modifyList(good, list(divergent = 5L)), 1.01, 400, 0.002))
  # a diagnostic that could not be computed is not a diagnostic that passed
  expect_false(rule(utils::modifyList(good, list(rhat = NA_real_))))
  expect_false(rule(utils::modifyList(good, list(ess_tail = NA_real_))))
  expect_false(rule(utils::modifyList(good, list(divergent = NA_integer_))))
})


test_that("the limits of the rule are checked", {
  expect_no_error(.check_rule(1.01, 400, 0.001))
  expect_error(.check_rule(1, 400, 0.001), class = "rTFT_bad_rule")
  expect_error(.check_rule("1.01", 400, 0.001), class = "rTFT_bad_rule")
  expect_error(.check_rule(1.01, 0, 0.001), class = "rTFT_bad_rule")
  expect_error(.check_rule(1.01, c(400, 400), 0.001), class = "rTFT_bad_rule")
  expect_error(.check_rule(1.01, 400, 2), class = "rTFT_bad_rule")
  expect_error(.check_rule(1.01, 400, NA_real_), class = "rTFT_bad_rule")
})


test_that("checking convergence without the posterior package is an error that says so", {
  local_mocked_bindings(.has_posterior = function() FALSE)
  expect_error(.fit_convergence(fits$exp), class = "rTFT_needs_posterior")
})


# ===========================================================================
# 4. Fitting the models of a table (fit_hazard_models)
# ===========================================================================

test_that("fit_hazard_models fits each model with what the table and the call give it", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  out <- quietly(fit_hazard_models(two_models, long, chains = 2, seed = 7))
  expect_named(out, c("fits", "loo", "kfold", "summary"))
  # one fit per model, each given its own row of the table and the arguments of the call
  expect_length(record$calls, 2)
  expect_equal(vapply(record$calls, \(cl) cl$m_form, character(1)), c("exp", "zeroexp"))
  expect_equal(vapply(record$calls, \(cl) cl$baseline, character(1)), c("mem", "mem"))
  for (cl in record$calls) {
    expect_equal(cl$dpar_formulas, two_models$dpar_formulas[[1]])
    expect_equal(cl$observation, "interval")
    expect_equal(cl$args, list(chains = 2, seed = 7))
  }
  # the fits and their cross-validation results, named by model
  expect_named(out$fits, c("exp_mem", "zeroexp_mem"))
  expect_identical(out$fits$exp_mem, fits$exp)
  expect_identical(out$fits$zeroexp_mem, fits$zeroexp)
  expect_named(out$loo, c("exp_mem", "zeroexp_mem"))
  expect_equal(out$loo$exp_mem$pointwise, loo_hazard(fits$exp)$pointwise)
  expect_null(out$kfold)
  # what is returned is what the weights and the averages read
  weights <- suppressMessages(hazard_weights(out$loo))
  expect_setequal(weights$model, c("exp_mem", "zeroexp_mem"))
  expect_no_error(suppressMessages(
    average_tdt_ttl(out$fits, weights, newdata = data.frame(mass = 10), t_ref = 30,
                    n_draws = 50)))
})


test_that("the summary has one row per model, with its diagnostics and its cross-validation", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  out <- quietly(fit_hazard_models(two_models, long, chains = 2))
  s <- out$summary
  expect_named(s, c("model", "m_form", "baseline", "status", "message", "chains", "draws",
                    "rhat", "ess_bulk", "ess_tail", "divergent", "converged",
                    "rhat_group", "ess_group", "elpd_loo", "se_loo", "high_k", "seconds"))
  expect_equal(s$model, two_models$model)
  expect_equal(s$status, c("fitted", "fitted"))
  expect_true(all(is.na(s$message)))
  expect_equal(s$chains, c(2, 2))
  expect_equal(s$draws, c(600, 600))
  conv <- lapply(fits, .fit_convergence)
  expect_equal(s$rhat, unname(vapply(conv, \(x) x$rhat, numeric(1))))
  expect_equal(s$ess_bulk, unname(vapply(conv, \(x) x$ess_bulk, numeric(1))))
  expect_equal(s$ess_tail, unname(vapply(conv, \(x) x$ess_tail, numeric(1))))
  expect_equal(s$divergent, unname(vapply(conv, \(x) x$divergent, integer(1))))
  expect_equal(s$elpd_loo, unname(vapply(out$loo, \(x) x$estimates["elpd_loo", "Estimate"],
                                         numeric(1))))
  expect_equal(s$se_loo, unname(vapply(out$loo, \(x) x$estimates["elpd_loo", "SE"],
                                       numeric(1))))
  expect_equal(s$high_k, c(0L, 0L))
  expect_true(all(is.na(s$rhat_group) & is.na(s$ess_group)))   # no group-level terms
  expect_true(all(s$seconds >= 0))
})


test_that("a fit that does not meet the convergence rule is kept, marked and named", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  # the stored fits are short, 300 draws a chain, and do not meet the default rule
  expect_warning(
    out <- suppressMessages(fit_hazard_models(two_models, long, chains = 2)),
    class = "rTFT_not_converged")
  expect_equal(out$summary$converged, c(FALSE, FALSE))
  expect_named(out$fits, c("exp_mem", "zeroexp_mem"))
  # under a rule they meet nothing is said, and under one that one of them meets
  # only the other is marked
  expect_no_warning(
    loose <- suppressMessages(fit_hazard_models(two_models, long, chains = 2, rhat = 1.5,
                                                ess = 50, divergent = 0.01)),
    class = "rTFT_not_converged")
  expect_equal(loose$summary$converged, c(TRUE, TRUE))
  between <- mean(loose$summary$rhat)
  one <- quietly(fit_hazard_models(two_models, long, chains = 2, rhat = between, ess = 50))
  expect_equal(one$summary$converged, loose$summary$rhat < between)
  expect_equal(sum(one$summary$converged), 1)
})


test_that("a model whose fit fails does not stop the others, and is reported", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  three <- rbind(
    build_hazard_models(c("exp", "power"), dpar_formulas = list(loga = ~ 1 + mass)),
    two_models[2, ])
  expect_warning(
    out <- suppressMessages(suppressWarnings(
      fit_hazard_models(three, long, chains = 2), classes = "rTFT_not_converged")),
    class = "rTFT_models_failed")
  expect_length(record$calls, 3)                         # the third was still fitted
  expect_equal(out$summary$status, c("fitted", "failed", "fitted"))
  expect_match(out$summary$message[2], "the sampler stopped")
  expect_true(all(is.na(out$summary$message[-2])))
  # the model that failed is in the summary, with nothing in place of what it lacks,
  # and in neither list
  failed <- out$summary[2, ]
  expect_true(all(is.na(failed[c("chains", "draws", "rhat", "ess_bulk", "ess_tail",
                                 "divergent", "converged", "rhat_group", "ess_group",
                                 "elpd_loo", "se_loo", "high_k")])))
  expect_named(out$fits, c("exp_mem", "zeroexp_mem"))
  expect_named(out$loo, c("exp_mem", "zeroexp_mem"))
})


test_that("a fit that lost a chain is refused", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  # the stored fits have two chains: asked for four, two were lost
  expect_error(.fit_one_model(two_models[1, ], long, "interval", "id", chains = 4,
                              dots = list()), class = "rTFT_lost_chain")
  expect_no_error(.fit_one_model(two_models[1, ], long, "interval", "id", chains = 2,
                                 dots = list()))
  # four chains is what brms runs when none is asked for
  expect_warning(out <- suppressMessages(fit_hazard_models(two_models, long)),
                 class = "rTFT_models_failed")
  expect_equal(out$summary$status, c("failed", "failed"))
  expect_match(out$summary$message, "Only 2 of the 4 chains")
  expect_length(out$fits, 0)
  expect_length(out$loo, 0)
})


test_that("the unit left out is settled once, for all the models, before any is fitted", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record, always = fit_trials))
  # the individuals of a trial share its effect, so trials are left out, and that is
  # said once and not once per model
  expect_equal(count_conditions(
    out <- fit_hazard_models(trial_models, long_trials, chains = 2), "rTFT_loo_unit"), 1)
  expect_equal(vapply(out$loo, \(x) x$unit, character(1)),
               c(exp_mem = "trial", exp_weibull = "trial"))
  expect_equal(nrow(out$loo$exp_mem$pointwise), 8)
  # a unit that is given is used, with one warning that the trial's effect stays in
  expect_equal(count_conditions(
    by_id <- fit_hazard_models(trial_models, long_trials, unit = "id", chains = 2),
    "rTFT_loo_leak"), 1)
  expect_equal(nrow(by_id$loo$exp_mem$pointwise), 48)
  # a warning about a later model is given even when the first has nothing to say
  later <- rbind(build_hazard_models("exp", "weibull"), trial_models[1, ])
  expect_equal(count_conditions(.models_unit(later, long_trials, "interval", "id"),
                                "rTFT_loo_leak"), 1)
  expect_equal(.models_unit(two_models, long, "interval", NULL), "id")
})


test_that("models that would leave out different units are refused before any is fitted", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record, always = fit_trials))
  mixed_units <- rbind(trial_models[1, ], build_hazard_models("exp", "weibull"))
  expect_error(suppressMessages(fit_hazard_models(mixed_units, long_trials, chains = 2)),
               class = "rTFT_units_differ")
  expect_length(record$calls, 0)
  # with the unit given they can be compared
  expect_no_error(
    quietly(fit_hazard_models(mixed_units, long_trials, unit = "id", chains = 2)))
  expect_length(record$calls, 2)
})


test_that("units for which importance sampling cannot be trusted are counted and reported", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record(), always = fit_trials))
  expect_warning(
    out <- suppressMessages(suppressWarnings(
      fit_hazard_models(trial_models, long_trials, chains = 2),
      classes = "rTFT_not_converged")),
    class = "rTFT_loo_unreliable")
  k <- out$loo$exp_mem$diagnostics$pareto_k
  expect_gt(sum(k > .k_threshold(600)), 0)
  expect_equal(out$summary$high_k, rep(sum(k > .k_threshold(600)), 2))
  # the trial effects are reported for a model that has them
  conv <- .fit_convergence(fit_trials)
  expect_equal(out$summary$rhat_group, rep(conv$rhat_group, 2))
  expect_equal(out$summary$ess_group, rep(conv$ess_group, 2))
  # the warning points to refitting, unless refitting was asked for
  said <- function(kfold) tryCatch(
    .report_models(out$summary, 1.5, 50, 0.01, kfold = kfold),
    rTFT_loo_unreliable = function(w) conditionMessage(w))
  expect_match(said(FALSE), "kfold = TRUE", fixed = TRUE)
  expect_no_match(said(TRUE), "kfold = TRUE", fixed = TRUE)
})


test_that("fit_hazard_models says which model it is fitting", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  expect_equal(count_conditions(fit_hazard_models(two_models, long, chains = 2),
                                "rTFT_models_progress"), 2)
})


test_that("fit_hazard_models refuses bad arguments before fitting anything", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  expect_error(fit_hazard_models(two_models[0, ], long), class = "rTFT_bad_models")
  expect_error(fit_hazard_models(two_models, long, kfold = "yes"), class = "rTFT_bad_kfold")
  expect_error(fit_hazard_models(two_models, long$Ta), class = "rTFT_missing_column")
  expect_error(fit_hazard_models(two_models, long[names(long) != "id"]),
               class = "rTFT_missing_column")
  expect_error(fit_hazard_models(two_models, long, observation = "nope"),
               class = "rTFT_bad_observation")
  expect_error(fit_hazard_models(two_models, long, rhat = 0.9), class = "rTFT_bad_rule")
  expect_error(fit_hazard_models(two_models, long, unit = "mass"), class = "rTFT_bad_unit")
  # what the table sets for each model cannot be set again for all
  expect_error(fit_hazard_models(two_models, long, m_form = "exp"),
               class = "rTFT_set_in_models")
  expect_error(
    fit_hazard_models(two_models, long, dpar_formulas = list(), costressors = list()),
    class = "rTFT_set_in_models")
  expect_length(record$calls, 0)
})


# ===========================================================================
# 5. Cross-validation by refitting, with the same folds for every model
# ===========================================================================

test_that("with kfold = TRUE each model is refitted, with the same folds", {
  refits <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record(), always = fit_trials),
                        .refit_fold = stand_in_refit(refits))
  set.seed(3)
  out <- quietly(fit_hazard_models(trial_models, long_trials, kfold = TRUE, K = 3,
                                   chains = 2, cores = 2, sd_prior = "exponential(1)"))
  expect_named(out$kfold, c("exp_mem", "exp_weibull"))
  expect_length(refits$calls, 6)                         # three folds for each model
  for (cv in out$kfold) {
    expect_s3_class(cv, c("kfold", "loo"))
    expect_equal(cv$unit, "trial")
    expect_equal(cv$K, 3)
  }
  # the folds are assigned once, so the two models are scored on the same splits
  expect_identical(out$kfold$exp_mem$folds, out$kfold$exp_weibull$folds)
  expect_equal(sort(as.vector(table(out$kfold$exp_mem$folds))), c(2, 3, 3))
  # a refit is given what brms takes, and keeps the priors of the fit it starts from
  for (cl in refits$calls) expect_equal(cl$args, list(chains = 2, cores = 2))
  # the summary gains the estimate by refitting, beside that by importance sampling
  expect_equal(names(out$summary)[17:20], c("high_k", "elpd_kfold", "se_kfold", "seconds"))
  expect_equal(out$summary$elpd_kfold, rep(-36, 2))      # minus 1 to minus 8, summed
  expect_equal(out$summary$se_kfold, rep(sqrt(8 * var(-(1:8))), 2))
  # and the results are those the weights read
  expect_no_error(suppressMessages(hazard_weights(out$kfold)))
})


test_that("one fold per unit is the default, and a unit given is left out by refitting too", {
  refits <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record(), always = fit_trials),
                        .refit_fold = stand_in_refit(refits))
  out <- quietly(fit_hazard_models(trial_models[1, ], long_trials, kfold = TRUE, chains = 2))
  expect_length(refits$calls, 8)                         # one refit per trial
  expect_equal(out$kfold$exp_mem$K, 8)
})


test_that("a number of folds that cannot be is refused before any model is fitted", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record, always = fit_trials))
  for (bad in list(1, 9, 2.5))
    expect_error(
      suppressMessages(fit_hazard_models(trial_models, long_trials, kfold = TRUE, K = bad,
                                         chains = 2)),
      class = "rTFT_bad_K")
  # 40 individuals: there cannot be 41 folds of them
  expect_error(suppressMessages(fit_hazard_models(two_models, long, kfold = TRUE, K = 41)),
               class = "rTFT_bad_K")
  expect_length(record$calls, 0)
})


test_that("a model whose refits fail is left out of the results by refitting, and reported", {
  refits <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record(), always = fit_trials),
                        .refit_fold = stand_in_refit(refits, fail_from = 3))
  set.seed(3)
  expect_warning(
    out <- suppressMessages(suppressWarnings(
      fit_hazard_models(trial_models, long_trials, kfold = TRUE, K = 2, chains = 2),
      classes = c("rTFT_not_converged", "rTFT_loo_unreliable"))),
    class = "rTFT_models_failed")
  expect_named(out$kfold, "exp_mem")
  expect_named(out$fits, c("exp_mem", "exp_weibull"))    # its fit is kept
  expect_equal(out$summary$status, c("fitted", "fitted"))
  expect_match(out$summary$message[2], "the refit stopped")
  expect_true(is.na(out$summary$message[1]))
  expect_equal(is.na(out$summary$elpd_kfold), c(FALSE, TRUE))
  expect_equal(is.na(out$summary$se_kfold), c(FALSE, TRUE))
})

