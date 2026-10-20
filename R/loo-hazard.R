# loo-hazard.R

# Cross-validation for a failure-rate model fitted one row per individual and
# interval at risk. A row is a slice of one individual's exposure, not an
# independent observation, so leaving out a row asks how well the model predicts
# the rest of an individual's exposure from its beginning, which is nearly always
# easy. The unit left out has to be at least the individual, and where individuals
# share a group-level effect in the model (the individuals of one trial, say), the
# group: leaving out one individual of a trial still leaves the trial's effect,
# estimated from the others, in the fit, and the prediction looks better than it is.

# ===========================================================================
# 1. Helpers
# ===========================================================================

# .grouping_terms() / .grouping_values() --------------------------------------

# The groupings of the model's group-level terms, named as brms reads them: `trial`
# for `(1 | trial)`; `trial` and `trial:tank` for `(1 | trial/tank)`, tanks nested
# in trials; `trial:tank` for `(1 | trial:tank)`; and for `(1 | gr(trial, by =
# sex))` only `trial`, since `by` sets a separate spread per sex and groups nothing.
# A label `p` in `(1 | p | trial)` lies to the left of the outer bar and is not a
# grouping.
.grouping_terms <- function(fit) {
  label <- function(e) {
    if (is.name(e)) return(as.character(e))
    if (is.call(e) && identical(e[[1]], quote(gr))) return(label(e[[2]]))
    if (is.call(e) && identical(e[[1]], quote(`:`)))
      return(paste(label(e[[2]]), label(e[[3]]), sep = ":"))
    paste(all.vars(e), collapse = ":")
  }
  # `a/b` is `a` and `a:b`; `a/b/c` adds `a:b:c`.
  expand <- function(e) {
    if (is.call(e) && identical(e[[1]], quote(`/`))) {
      outer <- expand(e[[2]])
      return(c(outer, paste(outer[length(outer)], label(e[[3]]), sep = ":")))
    }
    label(e)
  }
  terms <- purrr::map(fit$formula$pforms, \(f) {
    pieces <- .rhs_terms(f)
    bars   <- pieces[vapply(pieces, .is_bar, logical(1))]
    unlist(purrr::map(bars, \(b) expand(b[[2]][[3]])))
  })
  unique(unlist(terms, use.names = FALSE))
}

# The value of a grouping on each row of `data`. For a grouping of several columns
# (`trial:tank`) it is their combination, so that tank 1 of trial A and tank 1 of
# trial B are different groups, as brms estimates them, even where the tanks of
# each trial reuse the same names.
.grouping_values <- function(data, grouping) {
  cols <- strsplit(grouping, ":", fixed = TRUE)[[1]]
  if (length(cols) == 1L) return(data[[cols]])
  do.call(paste, c(unname(as.list(data[cols])), sep = ":"))
}


# .choose_unit() ---------------------------------------------------------------

# Choose the unit to leave out, from the fitted data `data`, the individual of each
# row `ids`, and the model's groupings `groups` (as .grouping_terms() names them);
# returns its name, the name messages use for it, and its value on each row.
#
# A grouping variable is shared when one of its levels holds more than one
# individual. With no unit given, the unit is the coarsest shared grouping variable
# that keeps each individual within one level (the trial, when individuals share a
# trial's effect), with a message saying so, and otherwise the individual. Nested
# groupings (e.g. multiple tanks within trials) are handled by that: leaving out a 
# trial takes its tanks with it. Crossed groupings are not (e.g., families spread 
# across  trials): no single unit removes the effects of both, and which to leave out
# depends on what is to be predicted, new trials or new families, so the unit must 
# then be given.
#
# A unit can be given by name, or as "id" for the individual. Where one level of a
# shared grouping variable spans more than one unit, that variable's effect stays in
# the fit, estimated from the rest of the level, and a warning says so: appropriate
# if the data to be predicted come from the same levels, optimistic if from new ones.
.choose_unit <- function(data, ids, groups, unit = NULL) {
  value <- stats::setNames(purrr::map(groups, \(g) .grouping_values(data, g)), groups)
  in_one_level <- function(values) {
    all(tapply(values, ids, \(v) length(unique(v))) == 1)
  }
  # TRUE when each level of `a` lies within one level of `b`.
  nested_in <- function(a, b) {
    all(tapply(value[[b]], value[[a]], \(v) length(unique(v))) == 1)
  }
  shared <- purrr::keep(groups, \(g) {
    any(tapply(ids, value[[g]], \(i) length(unique(i))) > 1)
  })

  if (is.null(unit)) {
    nesting <- purrr::keep(shared, \(g) in_one_level(value[[g]]))
    if (length(nesting)) {
      levels_of <- vapply(nesting, \(g) length(unique(value[[g]])), integer(1))
      unit <- nesting[which.min(levels_of)]
      crossed <- purrr::discard(setdiff(nesting, unit), \(g) nested_in(g, unit))
      if (length(crossed))
        cli::cli_abort(c(
          "The model has crossed groupings, {.field {c(unit, crossed)}}: leaving out a level of one leaves the effects of the others in the fit.",
          "i" = "Choose {.arg unit} by what you want to predict: new levels of {paste0(c(unit, crossed), ' (unit = \"', c(unit, crossed), '\")', collapse = ' or of ')}."
        ), class = "rTFT_crossed_groups")
      cli::cli_inform(c(
        "i" = "Leaving out one {.field {unit}} at a time: the individuals of each {.field {unit}} share its effect in the model, so leaving out one individual would leave that effect in the fit.",
        "i" = "Use {.code unit = \"id\"} to leave out individuals instead."
      ), class = "rTFT_loo_unit")
    } else {
      unit <- "id"
    }
  }

  # Only the individual or a grouping variable of the model can be a unit: a
  # covariate that happens to be constant within each individual (body mass, say)
  # would group individuals by its value, which is not a unit anyone leaves out.
  if (!is.character(unit) || length(unit) != 1L || !(unit == "id" || unit %in% groups))
    cli::cli_abort(c(
      "{.arg unit} must be {.val id} or one of the model's grouping variables.",
      "x" = "You supplied {.val {unit}}.",
      "i" = if (length(groups)) "The model's grouping variables are {.val {groups}}."
            else "The model has no grouping variables, so the unit can only be the individual."
    ), class = "rTFT_bad_unit")
  values <- if (unit == "id") ids else value[[unit]]
  if (unit != "id" && !in_one_level(values))
    cli::cli_abort(c(
      "Each individual must belong to one level of {.field {unit}} to be left out with it.",
      "x" = "Some individuals have rows in more than one level."
    ), class = "rTFT_bad_unit")

  # The unit as messages name it: "individual" rather than the column name `id`.
  label <- if (unit == "id") "individual" else unit
  # A shared grouping variable's effect stays in the fit when one of its levels spans
  # more than one unit: it is still estimated from the units kept. A coarser unit
  # takes whole levels with it, and leaves nothing behind.
  leaky <- purrr::keep(shared, \(g) {
    any(tapply(values, value[[g]], \(v) length(unique(v))) > 1)
  })
  if (length(leaky)) {
    # Suggest the coarsest leaking grouping: with nested groupings, leaving out its
    # levels takes the finer ones with them, so one change stops every leak.
    widest <- leaky[which.min(vapply(leaky, \(g) length(unique(value[[g]])), integer(1)))]
    cli::cli_warn(c(
      "!" = "Leaving out one {label} at a time leaves the effect of {.field {leaky}} in the fit, estimated from the rest of each level.",
      "i" = "That is appropriate if the data you want to predict come from the same levels of {.field {leaky}}, and optimistic if they come from new ones; then leave out whole levels, with {.code unit = \"{widest}\"}."
    ), class = "rTFT_loo_leak")
  }

  list(name = unit, label = label, values = values)
}


# .sum_by_unit() / .k_threshold() / .loo_units() -------------------------------

# Sum a draws-by-rows matrix of pointwise log-likelihood within each unit, giving a
# draws-by-units matrix with the units in order of first appearance. Given the
# parameters the likelihood factorises over the rows of a unit, so the sum is that
# unit's log-likelihood exactly, for every draw.
.sum_by_unit <- function(ll, values) {
  unit <- factor(as.character(values), levels = unique(as.character(values)))
  out <- vapply(split(seq_len(ncol(ll)), unit),
                \(j) rowSums(ll[, j, drop = FALSE]), numeric(nrow(ll)))
  matrix(out, nrow = nrow(ll), dimnames = list(NULL, levels(unit)))
}

# The largest Pareto k at which a PSIS estimate can be trusted, given `S` draws:
# min(1 - 1/log10(S), 0.7) (Vehtari et al. 2024), the threshold the loo package
# reports against, computed here because loo does not export it.
.k_threshold <- function(S) min(1 - 1 / log10(S), 0.7)

# PSIS-LOO over the units of a draws-by-units log-likelihood matrix, with a warning
# naming how many units it cannot be trusted for. Leaving out a large unit, such as
# a whole trial, can change the posterior too much for reweighting the draws of one
# fit to stand in for refitting without it.
.loo_units <- function(ll_unit, chain_id, unit_name, unit_label = unit_name) {
  # The relative efficiency is that of the likelihood, exp(ll). A whole unit's
  # log-likelihood can lie below about -745, where exp() is zero in double
  # precision, so each unit's values are centred on their largest first: the
  # efficiency is a ratio, which a constant factor per unit leaves unchanged.
  # Both on one core, whatever `mc.cores` says: with more, loo copies the process
  # with forks, which can hang for good once the threads of a loaded sampler are
  # running, and a matrix of one column per unit is quick on one core anyway.
  centred <- sweep(ll_unit, 2, apply(ll_unit, 2, max))
  r_eff   <- loo::relative_eff(exp(centred), chain_id = chain_id, cores = 1)
  out     <- suppressWarnings(loo::loo(ll_unit, r_eff = r_eff, cores = 1))
  limit <- .k_threshold(nrow(ll_unit))
  high  <- sum(out$diagnostics$pareto_k > limit)
  if (high)
    cli::cli_warn(c(
      "!" = "PSIS cannot be trusted for {high} of the {ncol(ll_unit)} units left out (one {unit_label} each): Pareto k is above {round(limit, 2)}.",
      "i" = "Leaving them out changes the fit too much to be approximated from it. Use {.fn kfold_hazard} for an exact estimate by refitting."
    ), class = "rTFT_loo_unreliable")
  out$unit <- unit_name
  out
}


# .check_fit_rows() ------------------------------------------------------------

# Check that the fit's record of which individual each row belongs to still lines
# up with its data. `st$id` has one entry per row, each individual repeated once for
# each of its rows. A different length means rows were dropped after the record was
# made (brms drops rows with missing values in any variable the model uses), and
# rows would be paired with the wrong individuals.
.check_fit_rows <- function(st, data) {
  if (length(st$id) != nrow(data))
    cli::cli_abort(c(
      "The fit's record of which individual each row belongs to does not match its data.",
      "x" = "{length(st$id)} row{?s} recorded for {nrow(data)} row{?s} in the fitted data.",
      "i" = "brms drops rows with missing values; remove them before fitting with {.fn fit_hazard}."
    ), class = "rTFT_rows_mismatch")
  invisible(st)
}


# ===========================================================================
# 2. Leave-one-unit-out for a fitted model
# ===========================================================================

# loo_hazard() -----------------------------------------------------------------

#' Leave-one-out cross-validation that leaves out whole individuals or groups
#'
#' Estimates how well a model fitted with [fit_hazard()] predicts data it was not
#' fitted to, by Pareto-smoothed importance sampling (PSIS-LOO), leaving out one
#' unit at a time: an individual, or a group of individuals.
#'
#' The model is fitted one row per individual and interval at risk, so
#' `brms::loo()` would leave out one interval, and ask how well the model predicts
#' the rest of an individual's exposure from its beginning, which is nearly always
#' easy. Here the pointwise log-likelihood is summed within each unit, which given
#' the parameters is that unit's log-likelihood exactly, and the unit is left out
#' whole.
#'
#' Where individuals share a group-level effect in the model, such as a random
#' intercept per trial, leaving out one individual still leaves its trial's effect,
#' estimated from the other individuals of the trial, in the fit, and the
#' prediction looks better than it would for a new trial. So by default the unit is
#' the coarsest grouping variable whose levels individuals share, with a message
#' saying so, and otherwise the individual. Nested groupings, such as tanks within
#' trials, are handled by this: leaving out a trial takes its tanks with it. Asking
#' for a finer unit gives a warning that the coarser grouping's effect stays in the
#' fit, which is appropriate when the data to be predicted come from the same
#' levels of it, and optimistic when they come from new ones.
#'
#' Crossed groupings, such as families spread across trials, have no default: no
#' single unit removes the effects of both, and leaving out trials answers how well
#' the model predicts new trials of known families, while leaving out families
#' answers how well it predicts new families. The unit must then be given.
#'
#' Leaving out a large unit can change the posterior too much for PSIS to stand in
#' for refitting. Where the Pareto k diagnostic of any unit is above the threshold
#' for the number of draws, a warning says so, and [kfold_hazard()] gives an exact
#' estimate by refitting.
#'
#' Models compared with [loo::loo_compare()] must be cross-validated with the same
#' unit.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param unit An optional character string, the unit to leave out: `"id"` for the
#'   individual, or one of the groupings of the model's group-level terms, such as
#'   `"trial"`, or `"trial:tank"` for the tanks of `(1 | trial/tank)`. `NULL` (the
#'   default) chooses as described above.
#' @return A `loo` object (see [loo::loo()]) with one pointwise value per unit, and
#'   an extra element `unit` naming the unit.
#' @references
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation
#' using leave-one-out cross-validation and WAIC. *Statistics and Computing* 27,
#' 1413-1432.
#'
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y. and Gabry, J. (2024). Pareto
#' smoothed importance sampling. *Journal of Machine Learning Research* 25(72),
#' 1-58.
#' @examples
#' \dontrun{
#' # a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
#' # degrees over 60 minutes and checked every minute, in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
#'                             checks = checks)
#'
#' fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + (1 | trial)))
#' loo_hazard(fit)                 # leaves out one trial at a time
#' }
#' @export
loo_hazard <- function(fit, unit = NULL) {
  st   <- .fit_settings(fit)
  data <- fit$data
  .check_fit_rows(st, data)
  u <- .choose_unit(data, st$id, .grouping_terms(fit), unit)

  ll <- brms::log_lik(fit)
  chain_id <- rep(seq_len(brms::nchains(fit)), each = nrow(ll) / brms::nchains(fit))
  .loo_units(.sum_by_unit(ll, u$values), chain_id, u$name, u$label)
}


# ===========================================================================
# 3. Leave-one-unit-out by refitting
# ===========================================================================

# .log_mean_exp() / .assign_folds() / .refit_fold() / .kfold_result() ----------
# .check_K() / .kfold_setup() / .kfold_refits()

# log(mean(exp(x))) for each column of a matrix, the log of a unit's predictive
# density averaged over the draws. The column's largest value is taken out first, so
# the exponentials neither overflow nor vanish for the very negative values a whole
# trial's log-likelihood can reach.
.log_mean_exp <- function(x) {
  top <- apply(x, 2, max)
  top + log(colMeans(exp(sweep(x, 2, top))))
}

# Assign each unit to one of K folds: a fold per unit when K is the number of units,
# and otherwise at random, as evenly as possible. Returns the fold of each unit,
# named by unit. Call set.seed() first for the same folds again.
.assign_folds <- function(units, K) {
  n <- length(units)
  fold <- if (K == n) seq_len(n) else sample(rep_len(seq_len(K), n))
  stats::setNames(fold, units)
}

# Refit the model to the rows in `train`, and return the log predictive density of
# each unit left out, averaged over the refit's draws, together with the refit. A
# unit's rows are summed within each draw before averaging, so the value is the
# unit's joint predictive density, as loo_hazard() estimates it. A unit's
# group-level effect is drawn from the fitted population distribution when the
# unit's level was left out entirely, so a left-out trial is predicted as a new
# trial.
#
# The refit is returned so the next fold can be refitted from it: a model loaded
# from a saved file has lost its compiled code and is compiled again on refitting,
# while a refit made in this session carries its compiled code, so refitting from
# it compiles nothing.
.refit_fold <- function(fit, data, train, unit_values, ...) {
  refit <- stats::update(fit, newdata = data[train, , drop = FALSE], refresh = 0, ...)
  ll <- brms::log_lik(refit, newdata = data[!train, , drop = FALSE],
                      allow_new_levels = TRUE, sample_new_levels = "gaussian")
  list(elpd = .log_mean_exp(.sum_by_unit(ll, unit_values[!train])), fit = refit)
}

# Assemble the result in the layout of a loo kfold object, so that
# loo::loo_compare() can use it: the summed estimate with its standard error over 
# the units, and one pointwise value per unit.
.kfold_result <- function(elpd, unit, K, folds) {
  n  <- length(elpd)
  se <- sqrt(n * stats::var(elpd))
  out <- list(
    estimates = cbind(Estimate = c(elpd_kfold = sum(elpd), kfoldic = -2 * sum(elpd)),
                      SE = c(elpd_kfold = se, kfoldic = 2 * se)),
    pointwise = cbind(elpd_kfold = elpd, kfoldic = -2 * elpd),
    unit = unit, K = K, folds = folds)
  class(out) <- c("kfold", "loo")
  out
}

# The number of folds for `n` units: one per unit unless given, and otherwise a whole
# number from 2 to `n`.
.check_K <- function(K, n) {
  if (is.null(K)) return(n)
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K != round(K) ||
      K < 2 || K > n)
    cli::cli_abort(c(
      "{.arg K} must be a whole number of folds from 2 to the number of units, {n}.",
      "x" = "You supplied {.val {K}}."
    ), class = "rTFT_bad_K")
  K
}

# Settle what the refits of a model need, before any is made, and return it: the
# fitted data, the unit left out with its value on each row, the units, the number
# of folds and the fold of each unit. The folds are assigned here unless `folds`
# gives them, as it does when several models are to be cross-validated with the same
# folds; the number of folds is then theirs, and `K` is not read.
.kfold_setup <- function(fit, unit = NULL, K = NULL, folds = NULL) {
  st   <- .fit_settings(fit)
  data <- fit$data
  .check_fit_rows(st, data)
  u <- .choose_unit(data, st$id, .grouping_terms(fit), unit)

  units <- unique(as.character(u$values))
  n <- length(units)
  if (!is.null(folds)) {
    if (!setequal(names(folds), units))
      cli::cli_abort(c(
        "The folds given are not those of this model's units.",
        "i" = "Models that share their folds must be fitted to the same table and leave out the same unit."
      ), class = "rTFT_bad_folds")
    # The refits run through the folds by number, from 1 to their count. A fold
    # whose number is skipped would be refitted with nothing left out, and the units
    # of a fold beyond the count would never be scored.
    if (!is.numeric(folds) || anyNA(folds) || !setequal(folds, seq_len(max(folds))))
      cli::cli_abort(
        "The folds given must be numbered from 1, with no number skipped.",
        class = "rTFT_bad_folds")
    return(list(data = data, unit = u, units = units, K = length(unique(folds)),
                folds = folds[units]))
  }
  K <- .check_K(K, n)
  list(data = data, unit = u, units = units, K = K, folds = .assign_folds(units, K))
}

# Refit the model once per fold of `setup`, each time without the units of the fold,
# and assemble the result. The first fold is refitted from the fit given and each
# later one from the refit before it.
.kfold_refits <- function(fit, setup, ...) {
  u <- setup$unit
  n <- length(setup$units)
  row_fold <- unname(setup$folds[as.character(u$values)])   # the fold of each row
  # A unit has no score until a refit gives it one, so its slot starts as missing:
  # a unit left unscored then makes the total missing, and cannot pass for a result.
  elpd <- stats::setNames(rep(NA_real_, n), setup$units)
  base <- fit                                     # refitted from, fold by fold
  for (k in seq_len(setup$K)) {
    cli::cli_inform(
      "Refitting without fold {k} of {setup$K}: {sum(setup$folds == k)} of the {n} units left out, one {u$label} each.",
      class = "rTFT_kfold_progress")
    fold <- .refit_fold(base, setup$data, row_fold != k, u$values, ...)
    elpd[names(fold$elpd)] <- fold$elpd
    base <- fold$fit
  }
  .kfold_result(elpd, u$name, setup$K, setup$folds)
}

# kfold_hazard() ---------------------------------------------------------------

#' Cross-validation by refitting, leaving out whole individuals or groups
#'
#' Estimates how well a model fitted with [fit_hazard()] predicts data it was not
#' fitted to, exactly, by refitting it with each unit left out in turn and
#' predicting the unit from the refit. It answers the same question as
#' [loo_hazard()], with the same choice of unit, but without its approximation: use
#' it where [loo_hazard()] warns that PSIS cannot be trusted, which is typical when
#' whole trials with a strong trial effect are left out.
#'
#' The units are split into `K` folds, and the model is refitted once per fold
#' without that fold's units. A unit's rows are summed within each draw of the
#' refit, so its value is its joint log predictive density, averaged over the
#' draws. Where the unit's level of a group-level term was left out entirely (a
#' whole trial), its effect is drawn from the fitted population distribution, so
#' the trial is predicted as a new trial would be.
#'
#' The refits take time: one model fit per fold. With `K` smaller than the number
#' of units, the units are assigned to folds at random; call [set.seed()] first to
#' get the same folds again. Models that are to be compared must be refitted with
#' the same folds: with fewer folds than units, call [set.seed()] with the same seed
#' before each call, or use [fit_hazard_models()], which assigns the folds once for
#' all the models.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param unit An optional character string, the unit to leave out, as in
#'   [loo_hazard()].
#' @param K An optional single whole number, the number of folds, from 2 to the
#'   number of units. `NULL` (the default) gives one fold per unit, leaving out each
#'   unit on its own.
#' @param ... Further arguments for the refits, passed to [brms::brm()] through
#'   [stats::update()], such as `chains`, `iter` or `cores`.
#' @return A `kfold` object in the layout of [loo::kfold()], with one pointwise value
#'   per unit, which [loo::loo_compare()] accepts, and extra elements `unit` (its
#'   name), `K`, and `folds` (the fold of each unit).
#' @examples
#' \dontrun{
#' # a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
#' # degrees over 60 minutes and checked every minute, in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
#'                             checks = checks)
#'
#' fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + (1 | trial)))
#' cv  <- kfold_hazard(fit, cores = 4)       # one refit per trial
#' }
#' @export
kfold_hazard <- function(fit, unit = NULL, K = NULL, ...) {
  .kfold_refits(fit, .kfold_setup(fit, unit, K), ...)
}
