# check-hazard.R

# Functions to check a fitted failure-rate model against the data it was fitted
# to. Convergence says the sampler found the posterior; it does not say the model
# describes the data. Two checks do. check_survival() sets the proportion of
# individuals observed to survive through time, in each group, against the same
# proportion in datasets simulated from the fitted model (a posterior predictive
# check). hazard_residuals() gives one residual per individual, the fraction the
# model predicts to have failed by the time it did, uniform between 0 and 1 when the
# model is right, to be examined against whatever might explain a
# misfit: temperature, time, a covariate or the trial.
#
# Both follow each individual through its own exposure, as recorded on the fit by
# fit_hazard(): the individual and the temperature trajectory of every row, and the
# elapsed time at which its interval ends. In both, the hazard accrued in an
# interval is exp(eta), with eta the model's linear predictor for that interval, so
# the probability of failing in it, given survival to its start, is
# 1 - exp(-exp(eta)). The group-level effects estimated for the individual's own
# groups are included.

# ===========================================================================
# 1. Following individuals through the fitted data
# ===========================================================================

# .fit_record() / .individual_columns()    --------------------------------------

# Read from a fit what the checks need, or stop when the fit has no such record (the
# table it was fitted to lacked the trajectory or the times of its intervals): for
# each row of the fitted data, the individual, the temperature trajectory and the
# elapsed time at the end of the interval; and for each individual, its
# trajectory, its rows in order of time, the time at which it failed or was last
# seen, and whether it failed.
.fit_record <- function(fit) {
  st <- .fit_settings(fit)
  if (is.null(st$tau_end) || is.null(st$path_id))
    cli::cli_abort(c(
      "This model carries no record of each interval's time and trajectory, which the checks need.",
      "i" = "Fit it with {.fn fit_hazard} to a table from {.fn build_hazard_data}, which holds both."
    ), class = "rTFT_no_check_record")
  data <- fit$data
  .check_fit_rows(st, data)
  rows <- split(seq_along(st$id), factor(st$id, levels = unique(st$id)))
  rows <- purrr::map(rows, \(r) r[order(st$tau_end[r])])
  last <- vapply(rows, \(r) r[length(r)], integer(1))
  list(st = st, data = data, rows = rows,
       individuals = tibble::tibble(
         id      = unique(st$id),
         path_id = as.character(st$path_id[last]),
         time    = st$tau_end[last],
         event   = as.integer(data$failed[last] == 1),
         last    = unname(last)))
}

# The columns of the fitted data that hold one value per individual (a covariate, a
# trial, a reference), as against those that change along its exposure (the
# temperature, the elapsed time, the response).
.individual_columns <- function(data, rows) {
  constant <- function(col) all(vapply(rows, \(r) length(unique(data[[col]][r])) == 1L,
                                       logical(1)))
  purrr::keep(setdiff(names(data), "failed"), constant)
}


# ===========================================================================
# 2. A posterior predictive check of survival
# ===========================================================================

# .km_survival() / .replicate_survival() ---------------------------------------

# The Kaplan-Meier estimate of survival at the times `at`, from each individual's
# time of failure or of last observation. A failure known only to an interval is
# placed at the end of its interval.
.km_survival <- function(time, event, at) {
  failures <- sort(unique(time[event == 1]))
  factor_at <- vapply(failures, \(t) 1 - sum(time == t & event == 1) / sum(time >= t),
                      numeric(1))
  vapply(at, \(t) prod(factor_at[failures <= t]), numeric(1))
}

# Simulate, for each posterior draw, whether each individual of one temperature
# trajectory is still unaffected at the end of each interval of the trajectory.
# `p_fail` holds the probability of failing in each interval given survival to its
# start, a matrix with a row per draw and the individuals' intervals side by side
# (all intervals of the first individual, then of the second, and so on). Returns an
# array of draws by individuals by intervals.
.replicate_survival <- function(p_fail, n_individuals, n_intervals) {
  n_draws <- nrow(p_fail)
  failed  <- array(stats::runif(length(p_fail)) < p_fail,
                   dim = c(n_draws, n_intervals, n_individuals))
  alive <- array(TRUE, dim = c(n_draws, n_individuals, n_intervals))
  still <- matrix(TRUE, n_draws, n_individuals)
  for (j in seq_len(n_intervals)) {       # each interval depends on the one before
    still <- still & !matrix(failed[, j, ], n_draws, n_individuals)
    alive[, , j] <- still
  }
  alive
}

# check_survival() -------------------------------------------------------------

#' Check a fitted model against the survival observed in each group
#'
#' Compares, for each group of individuals, the proportion observed to survive
#' through time with the proportion surviving in datasets simulated from the fitted
#' model: a posterior predictive check. Where the model describes the data, the
#' observed proportions lie within the band of the simulated ones.
#'
#' A group is the set of individuals that followed one temperature trajectory (one
#' `path_id` of the data), or, with `by`, a subset of them that share a value of
#' another variable, such as a trial. Individuals on different trajectories are
#' never pooled, since at the same elapsed time they are at different temperatures.
#'
#' The observed proportion is the Kaplan-Meier estimate, which allows for
#' individuals no longer observed; a failure known only to an interval is placed at
#' the end of the interval. For the simulated proportion, each posterior draw is used
#' to simulate every individual from the start of its trajectory to the last time
#' any individual on that trajectory was observed, with its own covariates and the
#' group-level effects estimated for its groups. The band therefore includes the
#' variation expected among samples of this size, which is wide for a group of ten.
#' The simulation is random: call [set.seed()] first for the same result again.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param by An optional character string naming a variable of the fitted model
#'   with one value per individual (a grouping such as a trial, or a covariate with
#'   few values), to split the individuals of each trajectory by.
#' @param prob A single number between 0 and 1, the coverage of the band of
#'   simulated proportions. Defaults to 0.95.
#' @param n_draws A single number, the number of posterior draws to simulate from,
#'   evenly spaced through the posterior. Defaults to 200.
#' @return A tibble with a row for each group and time: `path_id`, the `by` variable
#'   when given, `n` (the number of individuals in the group), `time`, `Ta` (the
#'   temperature of the trajectory in the interval ending at that time), `observed`
#'   (the Kaplan-Meier estimate of survival), and `median`, `lower` and `upper` of
#'   the simulated proportion surviving. Each group starts with survival 1 at time 0.
#' @seealso [plot_check()] to draw it, [hazard_residuals()] for residuals by
#'   individual.
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
#' fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#'
#' chk <- check_survival(fit)
#' plot_check(chk, time_unit = "min")
#'
#' # each trial on its own
#' plot_check(check_survival(fit, by = "trial"))
#' }
#' @export
check_survival <- function(fit, by = NULL, prob = 0.95, n_draws = 200) {
  .check_prob(prob)
  rec  <- .fit_record(fit)
  data <- rec$data
  ind  <- rec$individuals
  per_individual <- .individual_columns(data, rec$rows)
  if (!is.null(by) && !(is.character(by) && length(by) == 1L && by %in% per_individual))
    cli::cli_abort(c(
      "{.arg by} must name one variable of the model with one value per individual.",
      "x" = "You supplied {.val {by}}.",
      "i" = "Those in this model: {.val {setdiff(per_individual, c('Tr', 'tau_ref', 'Ta', 'tau', 'log_tau', 'exposure', 'log_exposure'))}}."
    ), class = "rTFT_bad_by")
  ind$.by <- if (is.null(by)) "" else as.character(data[[by]][ind$last])
  ids <- .thin_ids(brms::ndraws(fit), n_draws)
  probs <- c(0.5, (1 - prob) / 2, 1 - (1 - prob) / 2)

  out <- purrr::map(unique(ind$path_id), \(path) {
    members  <- which(ind$path_id == path)
    # The trajectory is read from the individual observed longest on it.
    template <- rec$rows[[members[which.max(ind$time[members])]]]
    times    <- rec$st$tau_end[template]
    n_int    <- length(template)

    # Every individual of the trajectory over all of its intervals: the intervals of
    # the trajectory, with the individual's own covariates and groups.
    nd <- data[rep(template, times = length(members)), , drop = FALSE]
    own <- rep(ind$last[members], each = n_int)
    nd[per_individual] <- data[own, per_individual, drop = FALSE]
    eta   <- brms::posterior_linpred(fit, newdata = nd, draw_ids = ids)
    alive <- .replicate_survival(1 - exp(-exp(eta)), length(members), n_int)

    purrr::map(unique(ind$.by[members]), \(level) {
      in_group <- which(ind$.by[members] == level)
      simulated <- apply(alive[, in_group, , drop = FALSE], c(1, 3), mean)
      q <- apply(simulated, 2, stats::quantile, probs = probs, names = FALSE)
      tibble::tibble(
        path_id  = path, .by = level, n = length(in_group),
        time     = c(0, times),
        Ta       = c(NA_real_, data$Ta[template]),
        observed = c(1, .km_survival(ind$time[members][in_group],
                                     ind$event[members][in_group], times)),
        median   = c(1, q[1, ]), lower = c(1, q[2, ]), upper = c(1, q[3, ]),
        prob     = prob)
    })
  })
  out <- dplyr::bind_rows(purrr::list_flatten(out))
  if (is.null(by)) out$.by <- NULL else names(out)[names(out) == ".by"] <- by
  out
}


# ===========================================================================
# 3. Residuals by individual
# ===========================================================================

# hazard_residuals() -----------------------------------------------------------

#' Residuals of a fitted model, one per individual
#'
#' Computes a randomised residual for each individual, in the manner of the
#' simulation-based residuals of the DHARMa package: the position of what was
#' observed within what the fitted model predicts for that individual. When the
#' model describes the data, the residuals are uniform between 0 and 1.
#'
#' The residual is the fraction of individuals that the model predicts to have
#' failed by the time this one did, 1 - S(t), with the individual's own covariates
#' and the group-level effects estimated for its groups, averaged over the posterior
#' draws. A residual near 0 is an individual that failed earlier than the model
#' predicts, and one near 1 an individual that lasted longer. What is observed is not
#' always a single time, so the residual is drawn at random within what was
#' observed:
#'
#' * a failure at a known time has the residual 1 - S(t);
#' * a failure known only to an interval has a residual drawn uniformly between
#'   that fraction at the start of the interval and at its end;
#' * a survivor, known only to have outlasted its last observation, has a residual
#'   drawn uniformly between that fraction at its last observation and 1.
#'
#' The draws are random: call [set.seed()] first for the same residuals again. One
#' residual per individual is used, and not one per interval, because the intervals
#' of an individual are slices of one exposure and are not separate observations.
#'
#' The residuals are to be examined against what was set before the outcome was
#' known: a covariate, a treatment, a trial. They should not be examined against the
#' time of failure, or the temperature at failure in a ramp. An individual that
#' failed late has, by that fact alone, a high predicted fraction failed at its time
#' of failure, so the residuals rise with the time of failure under any model, right
#' or wrong.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param n_draws A single number, the number of posterior draws to average the
#'   predicted survival over, evenly spaced through the posterior. Defaults to 400.
#' @return A tibble with one row per individual: `id`; `path_id`; the variables of
#'   the model that have one value per individual (covariates and groupings);
#'   `time`, the time at which it failed or was last observed; `Ta`, the temperature
#'   at that time; `event`, 1 for a failure and 0 for a survivor; and `residual`.
#' @references
#' Hartig, F. (2024). DHARMa: Residual diagnostics for hierarchical (multi-level /
#' mixed) regression models. R package version 0.4.7.
#' \doi{10.32614/CRAN.package.DHARMa}
#' @seealso [plot_residuals()] to draw them, [check_survival()] for a check by
#'   group.
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
#' fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#'
#' set.seed(1)
#' res <- hazard_residuals(fit)
#' plot_residuals(res)                      # against the uniform distribution
#' plot_residuals(res, against = "trial")   # by trial
#' }
#' @export
hazard_residuals <- function(fit, n_draws = 400) {
  rec  <- .fit_record(fit)
  data <- rec$data
  ind  <- rec$individuals
  ids  <- .thin_ids(brms::ndraws(fit), n_draws)

  # The hazard accrued in each interval, a row per draw and a column per interval,
  # summed within individual: over all of its intervals, and over all but its last.
  accrued  <- exp(brms::posterior_linpred(fit, draw_ids = ids))
  position <- match(rec$st$id, ind$id)
  total    <- t(rowsum(t(accrued), position, reorder = TRUE))
  before   <- total - accrued[, ind$last, drop = FALSE]
  # The fraction predicted to have failed by the last observation, and by the start
  # of the last interval, averaged over the draws.
  F_end    <- 1 - colMeans(exp(-total))
  F_start  <- 1 - colMeans(exp(-before))

  # Whether each individual's time is known exactly or to an interval.
  exact <- switch(rec$st$observation,
    interval   = rep(FALSE, nrow(ind)),
    continuous = rep(TRUE, nrow(ind)),
    mixed      = data$exact[ind$last] == 1)
  u <- stats::runif(nrow(ind))
  residual <- ifelse(ind$event == 0, F_end + u * (1 - F_end),
              ifelse(exact, F_end, F_start + u * (F_end - F_start)))

  # The variables to carry over: those with one value per individual, less the
  # model's own columns and any term brms stored under its expression, such as
  # `I(mass - 15)`, whose variable is carried in its own column.
  per_individual <- setdiff(.individual_columns(data, rec$rows),
                            c("Tr", "tau_ref", "exact", "Ta", "tau", "log_tau",
                              "exposure", "log_exposure"))
  per_individual <- per_individual[make.names(per_individual) == per_individual]
  dplyr::bind_cols(
    ind[c("id", "path_id")],
    tibble::as_tibble(data[ind$last, per_individual, drop = FALSE]),
    tibble::tibble(time = ind$time, Ta = data$Ta[ind$last], event = ind$event,
                   residual = residual))
}
