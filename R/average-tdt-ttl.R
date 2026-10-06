# average-tdt-ttl.R

# Model averaging for failure-rate models. With several candidate failure-rate and
# baseline functions fitted to the same data, cross-validation often does not single
# one out. Rather than report the best alone, the models are weighted by how well
# each predicts data it was not fitted to, and what is derived from them (sCTmax, the
# thermal death time curves, the tolerance landscape) is averaged by pooling their
# posterior draws in proportion to the weights. The pooled draws carry both the
# uncertainty within each model and the uncertainty about which model holds.

# ===========================================================================
# 1. Model weights
# ===========================================================================

# .pointwise_elpd() / .cv_weights() --------------------------------------------

# The log predictive density of each unit from a cross-validation result, named by
# unit, whether the result is from loo_hazard() or kfold_hazard().
.pointwise_elpd <- function(x) {
  column <- intersect(c("elpd_loo", "elpd_kfold"), colnames(x$pointwise))
  if (!inherits(x, "loo") || !length(column) || is.null(rownames(x$pointwise)))
    return(NULL)
  x$pointwise[, column[1]]
}

# The weights of the models in the columns of `lpd`, a units-by-models matrix of log
# predictive densities: by stacking, which finds the combination of models that
# predicts best, or as pseudo-BMA weights with the Bayesian bootstrap, which follow
# each model's own predictive density and its uncertainty over the units.
.cv_weights <- function(lpd, method) {
  as.numeric(switch(method,
    stacking  = loo::stacking_weights(lpd),
    pseudobma = loo::pseudobma_weights(lpd)))
}

# hazard_weights() -------------------------------------------------------------

#' Weights for averaging failure-rate models
#'
#' Weights several models fitted to the same data by how well each predicts data
#' it was not fitted to, from their cross-validation results. The models can be any
#' set: any failure-rate functions with any baselines. Which models are averaged is
#' chosen by which results are given.
#'
#' With `method = "stacking"` the weights are those of the combination of models
#' that predicts best. It is the method Yao et al. (2018) recommend, and the
#' default here. Stacking can give a model a weight of zero. This happens when two
#' models make nearly the same predictions: including both adds nothing to the
#' combination, so one takes the weight and the other gets none, even though it
#' would predict almost as well on its own. With `method = "pseudobma"` each model is
#' weighted by its own predictive density, allowing for its uncertainty over the
#' units by the Bayesian bootstrap; these weights are smoother, and steadier when the
#' units are few. The Bayesian bootstrap is random: call [set.seed()] first for the
#' same weights again.
#'
#' With `min_weight` below 1, only the highest-weighted models are kept, as many as
#' it takes for their weights to add up to at least `min_weight`, and the weights of
#' those kept are then computed again among themselves. Recomputing matters for
#' stacking, whose weights are found jointly for a set of models and are not simply
#' rescaled when a model leaves the set.
#'
#' The weights rest on one value per unit left out, so they are less certain the
#' fewer units there are, as when the unit is one of a handful of trials. A message
#' says how many units there were.
#'
#' Results from [loo_hazard()] and [kfold_hazard()] can be mixed, as when the
#' leave-one-out approximation is unreliable for one model and that model alone is
#' cross-validated by refitting. The two agree in what they estimate when each fold
#' leaves out a single unit. With fewer folds than units, each refit sees less of the
#' data than a leave-one-out fit does, so its predictions score slightly lower and
#' its model is at a small disadvantage. A message names the models this applies to.
#'
#' @param cv A named list of cross-validation results, one per model, from
#'   [loo_hazard()] or [kfold_hazard()]. All must have the same unit and the same
#'   units left out. The names identify the models.
#' @param method A character string, `"stacking"` (the default) or `"pseudobma"`.
#' @param min_weight A single number above 0 and at most 1, the cumulative weight
#'   the models kept must reach. The default, 1, keeps every model.
#' @return A tibble with one row per model, highest weight first: `model`; `elpd`,
#'   its summed log predictive density; `weight_all`, its weight among all the
#'   models, and `cumulative`, the running total; `kept`, whether it is kept; and
#'   `weight`, its weight among the models kept, which is the weight to average
#'   with (0 for a model not kept). The unit, the number of units and the method
#'   are attached as attributes.
#' @references
#' Yao, Y., Vehtari, A., Simpson, D. and Gelman, A. (2018). Using stacking to
#' average Bayesian predictive distributions. *Bayesian Analysis* 13, 917-1007.
#' @examples
#' \dontrun{
#' fits <- list("exp" = fit_exp, "zeroexp" = fit_zeroexp)
#' cv   <- lapply(fits, loo_hazard)
#' hazard_weights(cv)
#' hazard_weights(cv, min_weight = 0.95)
#' }
#' @export
hazard_weights <- function(cv, method = c("stacking", "pseudobma"), min_weight = 1) {
  method <- match.arg(method)
  if (!is.list(cv) || length(cv) < 2 || is.null(names(cv)) ||
      any(!nzchar(names(cv))) || anyDuplicated(names(cv)))
    cli::cli_abort(c(
      "{.arg cv} must be a named list of at least two cross-validation results.",
      "i" = "For example {.code list(exp = loo_hazard(fit_a), zeroexp = loo_hazard(fit_b))}."
    ), class = "rTFT_bad_cv")
  elpd <- purrr::map(cv, .pointwise_elpd)
  not_cv <- names(cv)[vapply(elpd, is.null, logical(1))]
  if (length(not_cv))
    cli::cli_abort(c(
      "Every element of {.arg cv} must be a result of {.fn loo_hazard} or {.fn kfold_hazard}.",
      "x" = "Not one: {.val {not_cv}}."
    ), class = "rTFT_bad_cv")
  if (!is.numeric(min_weight) || length(min_weight) != 1L || !is.finite(min_weight) ||
      min_weight <= 0 || min_weight > 1)
    cli::cli_abort(c(
      "{.arg min_weight} must be a single number above 0 and at most 1.",
      "x" = "You supplied {.val {min_weight}}."
    ), class = "rTFT_bad_min_weight")

  # Models can be weighted against each other only on the same predictions: the
  # same unit, and the same units left out.
  units <- purrr::map(cv, \(x) x$unit)
  same_units <- all(vapply(elpd, \(e) setequal(names(e), names(elpd[[1]])), logical(1)))
  if (length(unique(units)) != 1L || !same_units)
    cli::cli_abort(c(
      "The cross-validation results do not leave out the same units, so the models cannot be weighted against each other.",
      "x" = "Units: {paste(names(cv), unlist(units), sep = ' = ', collapse = ', ')}.",
      "i" = "Cross-validate every model with the same {.arg unit}, on fits of the same data."
    ), class = "rTFT_cv_mismatch")

  lpd <- vapply(elpd, \(e) e[names(elpd[[1]])], numeric(length(elpd[[1]])))
  n_units <- nrow(lpd)

  weight_all <- .cv_weights(lpd, method)
  ranked <- order(weight_all, decreasing = TRUE)
  # Keep the highest-weighted models until their weights reach min_weight (with a
  # small allowance for rounding in the weights).
  n_kept <- which(cumsum(weight_all[ranked]) >= min_weight - 1e-8)[1]
  kept   <- ranked[seq_len(n_kept)]
  weight <- numeric(length(cv))
  weight[kept] <- if (n_kept == length(cv)) weight_all[kept]
                  else if (n_kept == 1L) 1
                  else .cv_weights(lpd[, kept, drop = FALSE], method)

  out <- tibble::tibble(model = names(cv), elpd = unname(colSums(lpd)),
                        weight_all = weight_all, cumulative = NA_real_,
                        kept = seq_along(cv) %in% kept, weight = weight)[ranked, ]
  out$cumulative <- cumsum(out$weight_all)
  attr(out, "unit")    <- units[[1]]
  attr(out, "n_units") <- n_units
  attr(out, "method")  <- method

  # The number of fits each model's predictions came from: one per fold when
  # refitting, and in effect one per unit for leave-one-out. Models cross-validated
  # with fewer were predicted from less of the data than the others.
  folds <- vapply(cv, \(x) if (is.null(x$K)) n_units else as.numeric(x$K), numeric(1))
  fewer <- names(cv)[folds < max(folds)]
  if (length(fewer))
    cli::cli_inform(c(
      "!" = "The models were not cross-validated with the same number of folds.",
      "i" = "With fewer folds, each fit sees less of the data, which lowers a model's score slightly: {.val {fewer}}."
    ), class = "rTFT_mixed_cv")
  cli::cli_inform(c(
    "i" = "Weights by {method}, from {n_units} units left out (one {if (units[[1]] == 'id') 'individual' else units[[1]]} each)."
  ), class = "rTFT_weights")
  out
}


# ===========================================================================
# 2. Averaging what is derived
# ===========================================================================

# .allocate_draws() / .pool_draws() --------------------------------------------

# Share `n` draws among models in proportion to their weights, as whole numbers
# that add up to `n`: each model gets the whole part of its share, and the draws
# left over go to the models with the largest remainders.
.allocate_draws <- function(weight, n) {
  share <- weight / sum(weight) * n
  whole <- floor(share)
  extra <- order(share - whole, decreasing = TRUE)[seq_len(n - sum(whole))]
  whole[extra] <- whole[extra] + 1
  stats::setNames(as.integer(whole), names(weight))
}

# Pool the draws of several models, as .derive_draws() returns them, into one set in
# the same layout: for each row of newdata and each quantity, the models' draws one
# after another. The parameters are left out, since the parameters of different
# failure-rate functions are not the same quantities.
.pool_draws <- function(per_model) {
  first <- per_model[[1]]
  rows <- purrr::map(seq_along(first$rows), \(r) {
    quantity <- function(name, i) purrr::map(per_model, \(m) m$rows[[r]][[name]][[i]])
    list(
      sctmax  = purrr::map(seq_along(first$p),  \(i) unlist(quantity("sctmax", i),
                                                           use.names = FALSE)),
      tdt     = purrr::map(seq_along(first$p),  \(i) do.call(rbind, quantity("tdt", i))),
      surface = purrr::map(seq_along(first$Ta), \(j) do.call(rbind, quantity("surface", j))))
  })
  list(newdata = first$newdata, p = first$p, t_ref = first$t_ref, Ta = first$Ta,
       tau = first$tau, rows = rows)
}

# .average_setup() / .each_model() ---------------------------------------------

# Check the models and weights given for averaging and decide how many posterior
# draws each model contributes. Returns the models that carry weight, in the order
# of the weights table, the number of draws from each, and the weights table with
# that number added.
.average_setup <- function(fits, weights, n_draws) {
  if (!is.list(fits) || is.null(names(fits)) || anyDuplicated(names(fits)))
    cli::cli_abort("{.arg fits} must be a named list of fitted models.",
                   class = "rTFT_bad_fits")
  if (!is.data.frame(weights) || !all(c("model", "weight") %in% names(weights)))
    cli::cli_abort("{.arg weights} must be the table {.fn hazard_weights} returns.",
                   class = "rTFT_bad_weights")
  used <- weights[weights$weight > 0, , drop = FALSE]
  no_fit <- setdiff(used$model, names(fits))
  if (length(no_fit))
    cli::cli_abort(c(
      "{.arg fits} has no model named {.val {no_fit}}.",
      "i" = "The weights are for {.val {weights$model}}; {.arg fits} holds {.val {names(fits)}}."
    ), class = "rTFT_bad_fits")
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1)
    cli::cli_abort("{.arg n_draws} must be a single positive number of draws.",
                   class = "rTFT_bad_n_draws")
  fits <- fits[used$model]

  # Averaging is over models of the same data: the same individuals in the same rows.
  ids <- purrr::map(fits, \(f) .fit_settings(f)$id)
  if (!all(vapply(ids, identical, logical(1), ids[[1]])))
    cli::cli_abort(c(
      "The models were not fitted to the same data, so they cannot be averaged.",
      "i" = "Fit every model to the same table from {.fn build_hazard_data}."
    ), class = "rTFT_fits_mismatch")

  # Each model contributes its share of the pool. Where a model has fewer draws than
  # its share, the pool is made smaller so the shares stay in proportion.
  available <- vapply(fits, brms::ndraws, numeric(1))
  share     <- used$weight / sum(used$weight)
  n_pool    <- min(floor(n_draws), floor(min(available / share)))
  if (n_pool < floor(n_draws)) {
    short <- names(fits)[which.min(available / share)]
    cli::cli_inform(c(
      "!" = "Averaging over {n_pool} draws, not the {floor(n_draws)} asked for.",
      "i" = "Model {.val {short}} holds {available[[short]]} draws, fewer than its share of {floor(n_draws)} ({round(share[names(fits) == short] * floor(n_draws))}).",
      "i" = "Refit with more iterations for a larger pool."
    ), class = "rTFT_pool_reduced")
  }
  n_each    <- .allocate_draws(stats::setNames(used$weight, used$model), n_pool)

  weights$n_draws <- unname(n_each[weights$model])
  weights$n_draws[is.na(weights$n_draws)] <- 0L
  list(fits = fits[n_each > 0], n_each = n_each[n_each > 0], weights = weights)
}

# Compute something for each model in turn, giving it the model and its number of
# draws. A warning about `levels` the models do not use would be given by every
# model; such warnings are held back while the models are worked through, and each
# distinct one is then given once.
.each_model <- function(setup, compute) {
  held <- list()
  out <- withCallingHandlers(
    purrr::map(names(setup$fits), \(m) compute(setup$fits[[m]], setup$n_each[[m]])),
    rTFT_unused_levels = function(w) {
      held[[length(held) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  said <- vapply(held, conditionMessage, character(1))
  purrr::walk(held[!duplicated(said)], warning)
  out
}

# average_tdt_ttl() ------------------------------------------------------------

#' Thermal tolerance averaged over several fitted models
#'
#' Derives the static critical thermal maximum sCTmax, the thermal death time
#' curves and the thermal tolerance landscape as [derive_tdt_ttl()] does, averaged
#' over several models fitted to the same data. Each model's posterior draws of a
#' quantity are pooled in proportion to its weight, and the pooled draws are
#' summarised, so the median and credible interval reflect both the uncertainty
#' within each model and the uncertainty about which model holds.
#'
#' The parameters are not averaged. The parameters of different failure-rate
#' functions are different quantities under the same names, so a pooled a or b
#' would mean nothing; they are reported for each model by [derive_tdt_ttl()].
#'
#' The result has the layout of [derive_tdt_ttl()]'s, without its `params` table,
#' so [plot_tdt()], [plot_ttl()] and [plot_ttl_3d()] draw it as they draw a single
#' model's.
#'
#' @param fits A named list of `brmsfit`s, models fitted with [fit_hazard()] to the
#'   same data.
#' @param weights The weights from [hazard_weights()], whose model names are names
#'   of `fits`. Only the models it kept are averaged.
#' @inheritParams derive_tdt_ttl
#' @param n_draws A single number, the number of posterior draws in the pool. Each
#'   model contributes its share of them, evenly spaced through its posterior.
#'   Where a model holds fewer draws than its share, the pool is made smaller so
#'   that the shares stay in proportion to the weights, and a message says so.
#' @return A list of four tibbles: `sctmax`, `tdt_curves` and `ttl_surface` as from
#'   [derive_tdt_ttl()], and `weights`, the table of weights with the number of
#'   draws each model contributed.
#' @examples
#' \dontrun{
#' fits <- list("exp" = fit_exp, "zeroexp" = fit_zeroexp)
#' w    <- hazard_weights(lapply(fits, loo_hazard))
#' avg  <- average_tdt_ttl(fits, w, newdata = data.frame(mass = 10), t_ref = 60)
#' plot_tdt(avg)
#' }
#' @export
average_tdt_ttl <- function(fits, weights, newdata = NULL, p = c(0.1, 0.5, 0.9), t_ref,
                            Ta = NULL, tau = NULL, levels = NULL, prob = 0.95,
                            n_draws = 1000) {
  if (missing(t_ref)) .stop_no_t_ref()
  .check_prob(prob)
  setup <- .average_setup(fits, weights, n_draws)

  # The grids are those of the first model unless given, and the same for all.
  grid <- .fit_settings(setup$fits[[1]])
  if (is.null(Ta))  Ta  <- seq(grid$Ta_range[1], grid$Ta_range[2], length.out = 50)
  if (is.null(tau)) tau <- exp(seq(log(grid$tau_range[1]), log(grid$tau_range[2]),
                                   length.out = 50))

  per_model <- .each_model(setup, \(fit, n)
    .derive_draws(fit, newdata, p, t_ref, Ta, tau, levels, n))
  out <- .summarise_derived(.pool_draws(per_model), prob, params = FALSE)
  .report_unreached(out)
  c(out, list(weights = setup$weights))
}
