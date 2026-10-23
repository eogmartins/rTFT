# hazard-models.R

# A table of candidate models, and the fitting of all of them. Comparing models by
# cross-validation needs every candidate fitted to the same table with the same
# settings, checked for convergence by the same rule, and cross-validated with the
# same unit left out. Doing that by hand, one model after another, is where the
# candidates come to differ in ways that were not meant. Here the candidates are
# rows of a table, and one call fits them all and returns what the functions that
# weight and average models read.

# ===========================================================================
# 1. A table of models
# ===========================================================================

# .model_dpars() / .check_models() ---------------------------------------------

# The parameters of a model that can be given a formula: those of its failure-rate
# function, those of its baseline and the coefficients of its co-stressors.
.model_dpars <- function(m_form, baseline, costressors = list()) {
  c(.m_form_spec(m_form)$dpars, .baseline_spec(baseline)$dpars,
    .costressor_dpars(costressors))
}

# Check that `models` is a table of models as build_hazard_models() returns it. A
# user may bind two such tables, drop rows or edit a formula, so the table is
# checked again where it is used.
.check_models <- function(models) {
  columns <- c("model", "m_form", "baseline", "dpar_formulas", "costressors")
  if (!is.data.frame(models) || !all(columns %in% names(models)) || !nrow(models))
    cli::cli_abort(c(
      "{.arg models} must be a table of models, with at least one row.",
      "i" = "Build it with {.fn build_hazard_models}."
    ), class = "rTFT_bad_models")
  name <- models$model
  if (!is.character(name) || anyNA(name) || any(!nzchar(name)) || anyDuplicated(name))
    cli::cli_abort(c(
      "Each model of {.arg models} must have a name of its own in {.field model}.",
      "x" = if (anyDuplicated(name)) "Repeated: {.val {unique(name[duplicated(name)])}}."
    ), class = "rTFT_bad_models")
  if (!is.character(models$m_form) || !is.character(models$baseline))
    cli::cli_abort(c(
      "{.field m_form} and {.field baseline} must hold the names of the functions as text.",
      "i" = "A table read from a file may hold them as factors; convert them with {.fn as.character}."
    ), class = "rTFT_bad_models")
  if (!is.list(models$dpar_formulas) || !is.list(models$costressors))
    cli::cli_abort(
      "{.field dpar_formulas} and {.field costressors} must hold a list for each model.",
      class = "rTFT_bad_models")
  # plain loops, so that an unknown name gives its own error and not one wrapped
  # in the error of a mapping function
  for (m in unique(models$m_form)) .m_form_spec(m)
  for (b in unique(models$baseline)) .baseline_spec(b)
  invisible(models)
}


# build_hazard_models() --------------------------------------------------------

#' A table of failure-rate models to fit and compare
#'
#' Builds a table with one row per model, a failure-rate function with a baseline
#' and the formulas of its parameters, for [fit_hazard_models()] to fit. Each model
#' is assembled once here, so a misspelt name or a formula for a parameter no model
#' has is reported at once, before anything is fitted.
#'
#' With `cross = TRUE` every failure-rate function is combined with every baseline.
#' With `cross = FALSE` they are paired in the order given, for a set of models that
#' is not every combination. The table is an ordinary data frame: two tables can be
#' bound by rows, a row dropped, or the formulas of one model replaced, as long as
#' each model keeps a name of its own.
#'
#' The models differ in the parameters they have (a threshold function has `Tc`, a
#' Weibull baseline has `logshape`), so a formula in `dpar_formulas` is given to
#' each model that has the parameter and to no other. A formula for `loga`, which
#' every failure-rate function has, is therefore shared by all the models, and one
#' for `logshape` reaches only those with a Weibull or a log-logistic baseline.
#'
#' @param m_form A character vector, the failure-rate functions, each a key of
#'   [M_catalogue].
#' @param baseline A character vector, the baselines, each a key of [h0_catalogue].
#'   Defaults to `"mem"`, the memoryless baseline.
#' @param dpar_formulas A named list of one-sided formulas, one per parameter that
#'   depends on covariates or has a group-level effect, as in [build_hazard_spec()].
#'   Each is given to the models that have the parameter.
#' @param costressors An optional list of [costressor] objects, given to every
#'   model. Defaults to none.
#' @param cross A single logical value. `TRUE` (the default) combines every
#'   failure-rate function with every baseline; `FALSE` pairs them in the order
#'   given, a single value being repeated.
#' @param model_names An optional character vector, a name for each model, in the order
#'   of the rows. By default a model is named by its failure-rate function and its
#'   baseline, as in `"exp_weibull"`.
#' @return A data frame with one row per model and the columns `model` (its name),
#'   `m_form`, `baseline`, `dpar_formulas` and `costressors`, the last two holding a
#'   list for each model.
#' @seealso [fit_hazard_models()], which fits the models of the table.
#' @examples
#' # Two failure-rate functions with each of two baselines, all with a random
#' # intercept per trial on log a. The Weibull shape also depends on body mass, in
#' # the two models that have a Weibull baseline.
#' build_hazard_models(
#'   m_form = c("exp", "zeropower_threshold"), baseline = c("mem", "weibull"),
#'   dpar_formulas = list(loga = ~ 1 + (1 | trial), logshape = ~ 1 + mass))
#'
#' # A set that is not every combination, with names of its own.
#' build_hazard_models(
#'   m_form = c("exp", "exp", "zeropower_threshold"),
#'   baseline = c("mem", "llogis", "mem"), cross = FALSE,
#'   model_names = c("Exponential, memoryless", "Exponential, log-logistic",
#'             "Zero-power threshold, memoryless"))
#' @export
build_hazard_models <- function(m_form, baseline = "mem", dpar_formulas = list(),
                                costressors = list(), cross = TRUE,
                                model_names = NULL) {
  if (!is.character(m_form) || !length(m_form) || anyNA(m_form) ||
      !is.character(baseline) || !length(baseline) || anyNA(baseline))
    cli::cli_abort(
      "{.arg m_form} and {.arg baseline} must each be a character vector of at least one name.",
      class = "rTFT_bad_models")
  if (!is.logical(cross) || length(cross) != 1L || is.na(cross))
    cli::cli_abort("{.arg cross} must be {.code TRUE} or {.code FALSE}.",
                   class = "rTFT_bad_models")
  # look each name up for the error an unknown one raises; the entries are not kept
  for (m in unique(m_form)) .m_form_spec(m)
  for (b in unique(baseline)) .baseline_spec(b)
  .check_costressors(costressors)

  if (cross) {
    # each failure-rate function in turn, with every baseline
    grid <- expand.grid(baseline = baseline, m_form = m_form, stringsAsFactors = FALSE)
  } else {
    sizes <- c(length(m_form), length(baseline))
    if (sizes[1] != sizes[2] && !1L %in% sizes)
      cli::cli_abort(c(
        "With {.code cross = FALSE}, {.arg m_form} and {.arg baseline} are paired, so they must have the same length, or one of them a single value.",
        "x" = "{.arg m_form} has {length(m_form)} and {.arg baseline} has {length(baseline)}."
      ), class = "rTFT_bad_models")
    grid <- data.frame(baseline = baseline, m_form = m_form)
  }
  if (anyDuplicated(grid))
    cli::cli_abort(
      "A failure-rate function is given with the same baseline more than once.",
      class = "rTFT_bad_models")

  if (!is.list(dpar_formulas) ||
      (length(dpar_formulas) && (is.null(names(dpar_formulas)) ||
                                 any(!nzchar(names(dpar_formulas))))) ||
      !all(vapply(dpar_formulas, inherits, logical(1), "formula")))
    cli::cli_abort(c(
      "{.arg dpar_formulas} must be a named list of one-sided formulas.",
      "i" = "For example {.code list(loga = ~ 1 + (1 | trial))}."
    ), class = "rTFT_bad_models")
  dpars <- purrr::map2(grid$m_form, grid$baseline,
                       \(m, b) .model_dpars(m, b, costressors))
  for_none <- setdiff(names(dpar_formulas), unlist(dpars))
  if (length(for_none))
    cli::cli_abort(c(
      "{.arg dpar_formulas} has a formula for {.val {for_none}}, a parameter that none of the models has.",
      "i" = "The parameters of these models are {.val {unique(unlist(dpars))}}."
    ), class = "rTFT_formula_for_no_model")

  if (is.null(model_names)) model_names <- paste(grid$m_form, grid$baseline, sep = "_")
  if (!is.character(model_names) || length(model_names) != nrow(grid))
    cli::cli_abort(c(
      "{.arg model_names} must be a character vector with a name for each model.",
      "x" = "There are {nrow(grid)} model{?s} and {length(model_names)} name{?s}."
    ), class = "rTFT_bad_models")

  given  <- names(dpar_formulas)
  models <- tibble::tibble(
    model         = model_names,
    m_form        = grid$m_form,
    baseline      = grid$baseline,
    dpar_formulas = purrr::map(dpars, \(d) dpar_formulas[given %in% d]),
    costressors   = rep(list(costressors), nrow(grid)))
  .check_models(models)
  # Assemble each model as fit_hazard() will, and keep nothing: the call is there
  # for the error it raises when a model cannot be assembled, so that the error
  # comes now and not from inside a fit.
  for (i in seq_len(nrow(models)))
    build_hazard_spec(models$m_form[i], models$baseline[i], models$dpar_formulas[[i]],
                      costressors = models$costressors[[i]])
  models
}


# ===========================================================================
# 2. Fitting the models of a table
# ===========================================================================

# .has_posterior() / .check_rule() / .fit_convergence() / .meets_rule() --------

# Whether the posterior package is installed. .fit_convergence() is the one function
# that needs it, so it is a suggested package rather than a required one.
.has_posterior <- function() requireNamespace("posterior", quietly = TRUE)

# Check that `rhat`, `ess` and `divergent`, the three limits of the sampling rule,
# are single numbers a rule can use, and return TRUE invisibly.
.check_rule <- function(rhat, ess, divergent) {
  one <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!one(rhat) || rhat <= 1 || !one(ess) || ess <= 0 ||
      !one(divergent) || divergent < 0 || divergent > 1)
    cli::cli_abort(c(
      "The sampling rule needs {.arg rhat} above 1, a positive {.arg ess}, and {.arg divergent} as a proportion from 0 to 1.",
      "x" = "You supplied {.code rhat = {rhat}}, {.code ess = {ess}} and {.code divergent = {divergent}}."
    ), class = "rTFT_bad_rule")
  invisible(TRUE)
}

# Read from `fit` what the sampling rule needs, and return it as a list: the
# number of chains and of draws, the largest rank-normalised split-Rhat and the
# smallest bulk and tail effective sample sizes over the population-level
# coefficients and the standard deviations and correlations of the group-level
# effects, and the number of draws that came from a divergent transition. Also
# return, outside the rule, the largest Rhat and the smallest bulk or tail effective
# sample size over the group-level effects themselves (each trial's, say), or NA for
# a model that has none.
.fit_convergence <- function(fit) {
  if (!.has_posterior())
    cli::cli_abort(c(
      "The {.pkg posterior} package is needed to check the convergence of the fits.",
      "i" = "Install it with {.code install.packages(\"posterior\")}."
    ), class = "rTFT_needs_posterior")
  draws <- posterior::as_draws_array(fit)
  summary_of <- function(pattern) {
    chosen <- grep(pattern, posterior::variables(draws), value = TRUE)
    if (!length(chosen)) return(NULL)
    posterior::summarise_draws(posterior::subset_draws(draws, variable = chosen),
                               "rhat", "ess_bulk", "ess_tail")
  }
  conv  <- summary_of("^(b|sd|cor)_")
  group <- summary_of("^r_")
  nuts  <- brms::nuts_params(fit)
  list(chains     = brms::nchains(fit),
       draws      = brms::ndraws(fit),
       rhat       = max(conv$rhat),
       ess_bulk   = min(conv$ess_bulk),
       ess_tail   = min(conv$ess_tail),
       divergent  = as.integer(sum(nuts$Value[nuts$Parameter == "divergent__"])),
       rhat_group = if (is.null(group)) NA_real_ else max(group$rhat),
       ess_group  = if (is.null(group)) NA_real_
                    else min(group$ess_bulk, group$ess_tail))
}

# Whether the diagnostics in `conv`, as .fit_convergence() returns them, meet the
# rule set by `rhat`, `ess` and `divergent`: TRUE or FALSE. A diagnostic that could
# not be computed is missing, and a fit with a missing diagnostic has not been shown
# to meet the rule, so it gives FALSE.
.meets_rule <- function(conv, rhat, ess, divergent) {
  values <- c(conv$rhat, conv$ess_bulk, conv$ess_tail, conv$divergent, conv$draws)
  !anyNA(values) && conv$rhat < rhat && conv$ess_bulk > ess && conv$ess_tail > ess &&
    conv$divergent <= divergent * conv$draws
}


# .models_unit() / .say_again() / .models_plan() / .fit_one_model() / .model_row()

# The name of the unit to leave out in cross-validation, settled for all the models
# of `models` before any is fitted: from the groupings in their formulas and the
# table `long` they are fitted to, as loo_hazard() would choose it, or `unit` when it
# is given. The models are to be compared, so they must share it. What loo_hazard()
# would say about the unit for each model, its message and its warning, is held
# back and said once, whichever model it came from.
.models_unit <- function(models, long, observation, unit) {
  groups_of <- function(i) {
    spec <- build_hazard_spec(models$m_form[i], models$baseline[i],
                              models$dpar_formulas[[i]], observation,
                              costressors = models$costressors[[i]])
    .grouping_terms(list(formula = spec$formula))
  }
  held <- list()
  hold <- function(cnd) {
    held[[length(held) + 1L]] <<- cnd
    invokeRestart(if (inherits(cnd, "message")) "muffleMessage" else "muffleWarning")
  }
  chosen <- withCallingHandlers(
    vapply(seq_len(nrow(models)),
           \(i) .choose_unit(long, long$id, groups_of(i), unit)$name, character(1)),
    rTFT_loo_unit = hold, rTFT_loo_leak = hold)
  if (length(unique(chosen)) > 1L)
    cli::cli_abort(c(
      "The models do not share the unit to leave out in cross-validation, so they could not be compared.",
      "x" = "By their group-level terms it would be {.val {unique(chosen)}}.",
      "i" = "Give {.arg unit}, for every model."
    ), class = "rTFT_units_differ")
  # each distinct message or warning once, passed on as cli made it
  said <- vapply(held, conditionMessage, character(1))
  purrr::walk(held[!duplicated(said)], .say_again)
  chosen[1]
}

# Pass on `cnd`, a message or a warning that cli made and that was held back, as it
# was, and return nothing. rlang ends a message with a line ending when it prints it
# itself, but base message() prints a message's text as it is, so the line ending is
# added here; otherwise whatever is printed next would run on from the same line.
.say_again <- function(cnd) {
  if (!inherits(cnd, "message")) return(invisible(warning(cnd)))
  if (!endsWith(cnd$message, "\n")) cnd$message <- paste0(cnd$message, "\n")
  message(cnd)
}

# Check the arguments of a run of the models of `models` on the table `long`, settle
# what every model of the run shares, before anything is fitted, and return it as a
# list: the observation mode, the unit left out in cross-validation, the units of
# `long` (the distinct values of that unit) and the number of chains asked for in
# `dots`. With `kfold`, the number of folds `K` is checked against those units, so
# that a bad `K` stops the run now and not after the fits; without it, a `K` given
# is refused, since it would be ignored.
.models_plan <- function(models, long, observation, unit, kfold, K, rhat, ess,
                         divergent, dots) {
  .check_models(models)
  if (!is.data.frame(long) || !"id" %in% names(long))
    cli::cli_abort(c(
      "{.arg long} must be the table of intervals at risk, with its {.field id} column.",
      "i" = "Build it with {.fn build_hazard_data}."
    ), class = "rTFT_missing_column")
  observation <- .match_observation(observation)
  .check_rule(rhat, ess, divergent)
  if (!is.logical(kfold) || length(kfold) != 1L || is.na(kfold))
    cli::cli_abort("{.arg kfold} must be {.code TRUE} or {.code FALSE}.",
                   class = "rTFT_bad_kfold")
  if (!kfold && !is.null(K))
    cli::cli_abort(c(
      "{.arg K} was given with {.code kfold = FALSE}, so it would have no effect.",
      "i" = "Give {.code kfold = TRUE} for cross-validation by refitting in {.arg K} folds, or leave {.arg K} out."
    ), class = "rTFT_K_without_kfold")
  of_table <- intersect(names(dots),
                        c("m_form", "baseline", "dpar_formulas", "costressors"))
  if (length(of_table))
    cli::cli_abort(c(
      "{.arg {of_table}} cannot be given here: {?it is/they are} set for each model in {.arg models}.",
      "i" = "Change the table, or build it again with {.fn build_hazard_models}."
    ), class = "rTFT_set_in_models")
  unit  <- .models_unit(models, long, observation, unit)
  left_out <- if (unit == "id") long$id else .grouping_values(long, unit)
  units    <- unique(as.character(left_out))
  if (kfold) .check_K(K, length(units))
  list(observation = observation, unit = unit, units = units,
       chains = if (is.null(dots$chains)) 4 else dots$chains)
}

# Fit the model of `row`, one row of the table, to `long` with the arguments in
# `dots`; refuse it if fewer than `chains` chains finished; cross-validate it by
# importance sampling, leaving out `unit`; and return the fit, its cross-validation
# result, its convergence diagnostics and the number of units for which importance
# sampling cannot be trusted. brms summarises the chains that are left when one
# fails, so a fit that lost a chain looks complete unless the chains are counted.
.fit_one_model <- function(row, long, observation, unit, chains, dots) {
  fit <- do.call(fit_hazard, c(
    list(long = long, m_form = row$m_form, baseline = row$baseline,
         dpar_formulas = row$dpar_formulas[[1]], observation = observation,
         costressors = row$costressors[[1]]),
    dots))
  if (brms::nchains(fit) < chains)
    cli::cli_abort(
      "Only {brms::nchains(fit)} of the {chains} chains finished.",
      class = "rTFT_lost_chain")
  cv <- withCallingHandlers(
    loo_hazard(fit, unit = unit),
    rTFT_loo_leak = function(w) invokeRestart("muffleWarning"),
    rTFT_loo_unreliable = function(w) invokeRestart("muffleWarning"))
  list(fit = fit, loo = cv, convergence = .fit_convergence(fit),
       high_k = sum(cv$diagnostics$pareto_k > .k_threshold(brms::ndraws(fit))))
}

# One row of the summary, as a table of one row, for the model of `row`: a model that
# was fitted, when `res` holds what .fit_one_model() returned, or one that failed,
# when `res` is the error. The rule set by `rhat`, `ess` and `divergent` decides
# `rule_met`, and `seconds` is the time the model took.
.model_row <- function(row, res, seconds, rhat, ess, divergent) {
  failed <- inherits(res, "error")
  conv <- if (failed) list(chains = NA_integer_, draws = NA_integer_, rhat = NA_real_,
                           ess_bulk = NA_real_, ess_tail = NA_real_,
                           divergent = NA_integer_, rhat_group = NA_real_,
                           ess_group = NA_real_)
          else res$convergence
  elpd <- if (failed) c(NA_real_, NA_real_) else res$loo$estimates["elpd_loo", ]
  # The rule is applied here and not inside the table below, where `rhat` and
  # `divergent` would be read as the columns of those names and not as the limits.
  met <- if (failed) NA else .meets_rule(conv, rhat, ess, divergent)
  tibble::tibble(
    model      = row$model,
    m_form     = row$m_form,
    baseline   = row$baseline,
    status     = if (failed) "failed" else "fitted",
    message    = if (failed) cli::ansi_strip(conditionMessage(res)) else NA_character_,
    chains     = conv$chains,
    draws      = conv$draws,
    rhat       = conv$rhat,
    ess_bulk   = conv$ess_bulk,
    ess_tail   = conv$ess_tail,
    divergent  = conv$divergent,
    rule_met   = met,
    rhat_group = conv$rhat_group,
    ess_group  = conv$ess_group,
    elpd_loo   = elpd[[1]],
    se_loo     = elpd[[2]],
    high_k     = if (failed) NA_integer_ else as.integer(res$high_k),
    seconds    = seconds)
}


# fit_hazard_models() ----------------------------------------------------------

#' Fit, check and cross-validate the models of a table
#'
#' Fits each model of a table built with [build_hazard_models()] to the same data
#' with [fit_hazard()], checks the sampling of each fit by one rule, and
#' cross-validates each with [loo_hazard()], leaving out the same unit. With
#' `kfold = TRUE` each is also cross-validated by refitting, with
#' [kfold_hazard()] and the same folds for every model. What it returns is what
#' [hazard_weights()] and the functions that average models read.
#'
#' The unit left out is settled before any model is fitted, as [loo_hazard()]
#' chooses it, from the group-level terms of the models. Models that are compared
#' must be cross-validated with the same unit, so the models of the table must share
#' it; where their group-level terms would give different units, give `unit`.
#'
#' A model whose fit fails does not stop the others. It is left out of the fits
#' returned and is reported, with the reason, in the summary and in a warning. A fit
#' that lost one of its chains is treated as failed: the diagnostics of the chains
#' that are left would make it look complete.
#'
#' The sampling rule is applied to the population-level coefficients and to the standard
#' deviations and correlations of the group-level effects: the rank-normalised
#' split-Rhat of each below `rhat`, the bulk and the tail effective sample size of each
#' above `ess`, and at most the proportion `divergent` of the draws from divergent
#' transitions. The three check different things: Rhat whether the chains agree, the
#' effective sample sizes how precise the draws are, and the divergent transitions
#' whether the sampler could reach every region; the rule is not called convergence,
#' which is only the first of them. A fit that does not meet it is kept, since it is
#' often the one to look at, and is marked in the summary and named in a warning.
#' Whether to refit it, drop it or keep it is left to the user; zero divergent
#' transitions and a rule that is met are not, by themselves, evidence that a model was
#' sampled well. The limit of 400 on the effective sample sizes is 100 per chain for the
#' four chains of the default (Vehtari et al. 2021); with more chains, raise it in
#' proportion. The tolerance on divergent transitions exists so that a stray one does
#' not fail a fit that is otherwise well sampled; a divergent transition is still
#' evidence of a region the sampler could not enter, so where any remain, look at where
#' those draws fall against the group-level standard deviations before reading the fit,
#' and set `divergent = 0` to refuse any.
#'
#' The group-level effects themselves, such as the effect of each trial, are not
#' part of the rule: a model with many levels would take its worst case over many
#' parameters and exceed the limits by chance. They are reported beside it, as the
#' largest Rhat and the smallest effective sample size over them, because the
#' cross-validation by importance sampling, [check_survival()] and
#' [hazard_residuals()] use each trial's fitted effect. The derived quantities do
#' not: they are computed at the level of the population.
#'
#' @param models A data frame, the table of models from [build_hazard_models()].
#' @param long A data frame, the table of intervals at risk from
#'   [build_hazard_data()], to which every model is fitted.
#' @param observation A character string, how the failure times were observed, as in
#'   [fit_hazard()].
#' @param unit An optional character string, the unit to leave out in
#'   cross-validation, as in [loo_hazard()]. `NULL` (the default) chooses it as
#'   [loo_hazard()] does.
#' @param kfold A single logical value. `TRUE` also cross-validates each model by
#'   refitting, which takes one further fit per fold and model. Defaults to `FALSE`.
#' @param K An optional single whole number, the number of folds when `kfold` is
#'   `TRUE`, as in [kfold_hazard()]. `NULL` (the default) gives one fold per unit.
#'   With fewer folds than units the units are assigned to folds at random, once for
#'   all the models; call [set.seed()] first to get the same folds, and the same
#'   estimates by refitting, again (see [kfold_hazard()]).
#' @param rhat A single number above 1, the limit below which the Rhat of every
#'   parameter must lie. Defaults to 1.01.
#' @param ess A single positive number, the limit above which the bulk and the tail
#'   effective sample size of every parameter must lie. Defaults to 400.
#' @param divergent A single number from 0 to 1, the largest proportion of the draws
#'   that may come from divergent transitions. Defaults to 0.001.
#' @param ... Further arguments for every fit, passed to [fit_hazard()] and through
#'   it to [brms::brm()], such as `chains`, `iter`, `cores`, `seed` and `backend`,
#'   or the priors that [fit_hazard()] takes. Those that [brms::brm()] takes are
#'   given to the refits of `kfold = TRUE` as well.
#' @return A list with four elements. `fits` is a named list of the fitted models
#'   (`brmsfit`), without those that failed. `loo` is a named list of their
#'   cross-validation results from [loo_hazard()]. `kfold` is a named list of their
#'   results from [kfold_hazard()], or `NULL` when `kfold` is `FALSE`. `summary` is a
#'   data frame with one row per model of the table: its name, failure-rate function
#'   and baseline; `status` (`"fitted"` or `"failed"`) and, for a model that failed,
#'   the reason in `message`; the number of `chains` and `draws`; the largest `rhat`,
#'   the smallest `ess_bulk` and `ess_tail` and the number of `divergent` draws;
#'   `rule_met`, whether the sampling rule is met; `rhat_group` and `ess_group`, the
#'   largest Rhat and the smallest bulk or tail effective sample size over the
#'   group-level effects, outside the rule (`NA` for a model without them); `elpd_loo`
#'   with its standard error `se_loo`, and `high_k`, the number of units for which
#'   importance sampling cannot be trusted; `elpd_kfold` and `se_kfold` when `kfold` is
#'   `TRUE`; and `seconds`, the time the fit and its cross-validation took.
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: an improved Rhat for assessing
#' convergence of MCMC. *Bayesian Analysis* 16, 667-718.
#' \doi{10.1214/20-BA1221}
#' @seealso [hazard_weights()] to weight the models, and [average_tdt_ttl()],
#'   [average_dctmax()] and [average_survival()] to average them.
#' @examples
#' \dontrun{
#' # a small ramping trial: 20 individuals warmed from 20 to 35 degrees over 60
#' # minutes and checked every minute, in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1,
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = "trial",
#'                             checks = checks)
#'
#' # the exponential failure rate with two baselines, each with a random intercept
#' # per trial on log a
#' models <- build_hazard_models("exp", c("mem", "weibull"),
#'                               dpar_formulas = list(loga = ~ 1 + (1 | trial)))
#' out <- fit_hazard_models(models, long, chains = 4, iter = 2000, cores = 4, seed = 1)
#' out$summary
#'
#' # weight the models by how well each predicts a trial it was not fitted to, and
#' # average what is derived from them
#' weights <- hazard_weights(out$loo)
#' average_tdt_ttl(out$fits, weights, t_ref = 60)$sctmax
#' }
#' @export
fit_hazard_models <- function(models, long,
                              observation = c("interval", "continuous", "mixed"),
                              unit = NULL, kfold = FALSE, K = NULL,
                              rhat = 1.01, ess = 400, divergent = 0.001, ...) {
  dots <- list(...)
  plan <- .models_plan(models, long, observation, unit, kfold, K, rhat, ess, divergent,
                       dots)
  observation <- plan$observation
  unit        <- plan$unit
  chains      <- plan$chains

  n <- nrow(models)
  results <- vector("list", n)
  rows    <- vector("list", n)
  for (i in seq_len(n)) {
    cli::cli_inform("Fitting model {i} of {n}, {.val {models$model[i]}}.",
                    class = "rTFT_models_progress")
    started <- proc.time()[["elapsed"]]
    results[[i]] <- tryCatch(
      .fit_one_model(models[i, ], long, observation, unit, chains, dots),
      error = function(e) e)
    seconds   <- proc.time()[["elapsed"]] - started
    rows[[i]] <- .model_row(models[i, ], results[[i]], seconds, rhat, ess, divergent)
  }
  model_summary <- purrr::list_rbind(rows)
  fitted  <- model_summary$status == "fitted"
  named   <- function(what) stats::setNames(purrr::map(results[fitted], what),
                                            models$model[fitted])
  out <- list(fits = named("fit"), loo = named("loo"), kfold = NULL,
              summary = model_summary)

  if (kfold) {
    cv <- .kfold_models(out$fits, unit, K, dots)
    out$kfold   <- cv$kfold
    out$summary <- .kfold_columns(out$summary, cv$kfold, cv$errors)
  }

  .report_models(out$summary, rhat, ess, divergent, kfold)
  out
}


# .refit_args() / .kfold_models() / .kfold_columns() / .report_models() ---------

# The arguments of `dots`, those given for every fit, that are given to a refit: those
# brms takes. Those that only fit_hazard() takes, the priors and the steepness of a
# threshold, are left out, since a refit keeps them from the fit it starts from.
.refit_args <- function(dots) {
  dots[setdiff(names(dots), c("prior", "intercept_priors", "slope_prior", "sd_prior",
                              "k_sw"))]
}

# Cross-validate each model of `fits` by refitting, leaving out `unit` in `K` folds,
# and return a list of two elements: `kfold`, the results by model, and `errors`, the
# reason each model whose refits failed was left out. The folds are assigned once,
# for the first model cross-validated, and given to every model after it. The refits
# are given the arguments of `dots` that .refit_args() keeps.
.kfold_models <- function(fits, unit, K, dots) {
  refit_args <- .refit_args(dots)
  quietly <- function(expr) withCallingHandlers(
    expr, rTFT_loo_leak = function(w) invokeRestart("muffleWarning"))
  folds  <- NULL
  kfold  <- list()
  errors <- character(0)
  for (m in names(fits)) {
    cli::cli_inform("Cross-validating {.val {m}} by refitting.",
                    class = "rTFT_models_progress")
    res <- tryCatch({
      setup <- quietly(.kfold_setup(fits[[m]], unit, K, folds))
      do.call(.kfold_refits, c(list(fits[[m]], setup), refit_args))
    }, error = function(e) e)
    if (inherits(res, "error")) {
      errors[[m]] <- cli::ansi_strip(conditionMessage(res))
    } else {
      kfold[[m]] <- res
      folds <- res$folds
    }
  }
  list(kfold = kfold, errors = errors)
}

# Add to `model_summary` the estimate by refitting of each model in `kfold`, a named
# list of results in the layout of kfold_hazard(), as the columns `elpd_kfold` and
# `se_kfold` after `high_k` (missing for a model without a result), and put in
# `message` the reason each model named in `errors`, a named character vector, was
# not cross-validated by refitting, with a warning that names those models. Return
# the summary.
.kfold_columns <- function(model_summary, kfold, errors) {
  estimate <- function(m, column) {
    if (is.null(kfold[[m]])) return(NA_real_)
    kfold[[m]]$estimates["elpd_kfold", column]
  }
  model_summary$elpd_kfold <- vapply(model_summary$model, estimate, numeric(1),
                                     "Estimate", USE.NAMES = FALSE)
  model_summary$se_kfold   <- vapply(model_summary$model, estimate, numeric(1), "SE",
                                     USE.NAMES = FALSE)
  model_summary <- dplyr::relocate(model_summary, "elpd_kfold", "se_kfold",
                                   .after = "high_k")
  not_refitted <- model_summary$model %in% names(errors)
  model_summary$message[not_refitted] <-
    unname(errors[model_summary$model[not_refitted]])
  if (length(errors))
    cli::cli_warn(c(
      "!" = "{length(errors)} model{?s} could not be cross-validated by refitting: {.val {names(errors)}}.",
      "i" = "The reason for each is in {.field message} of the summary."
    ), class = "rTFT_models_failed")
  model_summary
}

# Say, from `model_summary`, which models failed, which did not meet the convergence
# rule set by `rhat`, `ess` and `divergent`, and for which importance sampling cannot
# be trusted, each in one warning that names them, and return `model_summary`
# invisibly. The last warning points to cross-validation by refitting, or to its
# results when `kfold` says they were asked for.
.report_models <- function(model_summary, rhat, ess, divergent, kfold = FALSE) {
  fitted <- model_summary$status == "fitted"
  failed <- model_summary$model[model_summary$status == "failed"]
  if (length(failed))
    cli::cli_warn(c(
      "!" = "{length(failed)} of the {nrow(model_summary)} model{?s} could not be fitted: {.val {failed}}.",
      "i" = "The reason for each is in {.field message} of the summary.",
      # not said of a collection from a cluster with models still missing
      "i" = if (any(fitted) && sum(fitted) + length(failed) == nrow(model_summary))
        "The other models were fitted."
    ), class = "rTFT_models_failed")
  unmet <- model_summary$model[fitted & !model_summary$rule_met]
  if (length(unmet))
    cli::cli_warn(c(
      "!" = "{length(unmet)} model{?s} did not meet the sampling rule: {.val {unmet}}.",
      "i" = "The rule is Rhat below {rhat}, bulk and tail effective sample sizes above {ess}, and at most {divergent * 100}% of the draws from divergent transitions.",
      "i" = "{cli::qty(length(unmet))}{?It is/They are} kept in the fits returned; see {.field rhat}, {.field ess_bulk}, {.field ess_tail} and {.field divergent} in the summary."
    ), class = "rTFT_rule_not_met")
  high <- model_summary$model[fitted & model_summary$high_k > 0]
  if (length(high))
    cli::cli_warn(c(
      "!" = "PSIS cannot be trusted for some of the units left out, in {length(high)} model{?s}: {.val {high}}.",
      "i" = "The number of such units is in {.field high_k} of the summary.",
      "i" = if (kfold) "The estimates by refitting, in {.field kfold}, do not rely on it."
            else "Use {.code kfold = TRUE} for an estimate by refitting, without the approximation."
    ), class = "rTFT_loo_unreliable")
  invisible(model_summary)
}
