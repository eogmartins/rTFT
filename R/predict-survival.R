# predict-survival.R

# Functions to predict survival along a measured record of temperature from a fitted
# failure-rate model, with the posterior uncertainty carried through. The assay
# estimates how the failure rate depends on temperature; the record is a measurement
# of temperatures individuals of the study species potenially encounter and for how 
# long. Between two readings the cumulative hazard grows by M(T) [H0(tau_1) - H0(tau_0)],
# the failure rate at the temperature of that interval times the growth of the 
# cumulative baseline hazard over it, and survival to a time is exp(-cumulative hazard
# up to it). With the memoryless baseline the growth of H0 is the length of the interval,
# and the cumulative hazard is the failure rate summed over the record.
#
# The path is predict_survival(), which calls .record_draws() to compute, for each
# posterior draw, the hazard accrued in each interval of the record, and
# .summarise_record() to reduce the draws to survival through time, the probability
# of failing by the end, and each period's share of the risk. average_survival()
# pools the draws of several models by their weights before summarising.

# ===========================================================================
# 1. The record and the hazard accrued along it
# ===========================================================================

# .check_record() / .record_intervals() ----------------------------------------

# Check a record of temperature through time: the columns it needs, at least two
# readings, no missing values, and times that increase.
.check_record <- function(record, columns) {
  if (!is.data.frame(record))
    cli::cli_abort("{.arg record} must be a data frame of times and temperatures.",
                   class = "rTFT_bad_record")
  absent <- setdiff(columns, names(record))
  if (length(absent))
    cli::cli_abort(c(
      "{.arg record} is missing column{?s} {.field {absent}}.",
      "i" = "It needs {.field time} and {.field temp}, as the {.arg path} of {.fn build_hazard_data} does."
    ), class = "rTFT_bad_record")
  if (nrow(record) < 2L)
    cli::cli_abort("{.arg record} needs at least two readings.", class = "rTFT_bad_record")
  numeric_ok <- vapply(record[c("time", "temp")], \(x) is.numeric(x) && all(is.finite(x)),
                       logical(1))
  if (!all(numeric_ok) || anyNA(record[columns]))
    cli::cli_abort(c(
      "{.arg record} has missing or non-numeric values in the columns the prediction uses.",
      "i" = "Fill gaps in the record before predicting; a gap cannot be told from an interval at risk."
    ), class = "rTFT_bad_record")
  if (any(diff(record$time) <= 0))
    cli::cli_abort("The times in {.arg record} must increase from one reading to the next.",
                   class = "rTFT_bad_record")
  invisible(record)
}

# The intervals between consecutive readings of a record: the elapsed exposure at
# their start and end, counted from the first reading, and the value at their
# midpoint of the temperature and of each further column asked for, taken halfway
# between the two readings as build_hazard_data() does for the fitted data.
.record_intervals <- function(record, columns = character(0)) {
  n   <- nrow(record)
  mid <- function(x) (x[-n] + x[-1]) / 2
  c(list(tau_start = record$time[-n] - record$time[1],
         tau_end   = record$time[-1] - record$time[1],
         Ta        = mid(record$temp)),
    purrr::map(stats::setNames(columns, columns), \(col) mid(record[[col]])))
}

# .refine_intervals() ----------------------------------------------------------

# Divide the intervals between the readings of a record into parts over which
# temperature changes by at most `grid_dT`, as build_hazard_data() divides the
# intervals at risk: the hazard is taken as constant within a part, at its midpoint,
# which is accurate only while temperature changes little within it. Between two
# readings the temperature, and any co-stressor, is taken to change linearly.
# Returns the parts in the layout of .record_intervals(), with `parent`, the
# interval between readings each part belongs to.
.refine_intervals <- function(record, columns = character(0), grid_dT = 0.1) {
  n      <- nrow(record)
  start  <- record$time[-n] - record$time[1]
  end    <- record$time[-1] - record$time[1]
  parts  <- pmax(1L, ceiling(abs(diff(record$temp)) / grid_dT - 1e-9))
  cut    <- .subdivide(start, end, parts)
  parent <- rep(seq_len(n - 1L), parts)
  # how far through its parent interval the midpoint of each part lies
  along  <- ((cut$ts + cut$te) / 2 - start[parent]) / (end - start)[parent]
  at_mid <- function(x) x[-n][parent] + along * diff(x)[parent]
  c(list(tau_start = cut$ts, tau_end = cut$te, Ta = at_mid(record$temp), parent = parent),
    purrr::map(stats::setNames(columns, columns), \(col) at_mid(record[[col]])))
}

# .record_hazard() -------------------------------------------------------------

# The hazard accrued in each interval of a record, for each posterior draw: a matrix
# with a row per draw and a column per interval, M(T) [H0(tau_end) - H0(tau_start)].
# `st` holds the settings of the fit, `d` one vector of draws per parameter,
# `intervals` what .record_intervals() returns, and `levels` the level of each
# co-stressor the record does not follow through time.
#
# The failure rates that are zero at the reference temperature (zeroexp, zeropower)
# are taken as zero below it too, which is what being zero at the reference means for
# a record that goes colder than the assay began. The power function has a positive
# rate at the reference and no defined one below, so a record that goes below it is
# refused.
.record_hazard <- function(st, d, intervals, levels = NULL) {
  nd <- length(d$loga)
  nT <- length(intervals$Ta)
  cs <- st$costressors
  Ta <- intervals$Ta
  if (st$m_form == "power" && any(Ta < st$Tr))
    cli::cli_abort(c(
      "The {.val power} failure rate is not defined below the reference temperature.",
      "x" = "The record goes down to {.val {round(min(Ta), 2)}}; the reference is {.val {st$Tr}}."
    ), class = "rTFT_below_Tr")
  if (st$m_form %in% c("zeroexp", "zeropower")) Ta <- pmax(Ta, st$Tr)

  k     <- if (M_catalogue[[st$m_form]]$k_sw != 0) st$k_sw
  setup <- .rate_setup(st$m_form, d, st$Tr, k, cs, levels)
  # A co-stressor the record follows through time takes its value in each interval,
  # standardised as the fit standardised its data, in place of a single level.
  for (x in cs) {
    if (is.null(intervals[[x$column]])) next
    setup$z[[paste0("z_", x$name)]] <-
      matrix(.costressor_z(x, intervals[[x$column]]), nd, nT, byrow = TRUE)
  }
  M  <- exp(.log_M_eval(setup, matrix(Ta, nd, nT, byrow = TRUE)))
  bl <- .baseline_draws(st$baseline, d, nd, st$tau_ref)
  H0 <- H0_of(st$baseline)
  growth <- H0(matrix(intervals$tau_end, nd, nT, byrow = TRUE), bl$shape, bl$scale) -
            H0(matrix(intervals$tau_start, nd, nT, byrow = TRUE), bl$shape, bl$scale)
  M * growth
}

# .record_zone() / .record_coverage() ------------------------------------------

# Whether each interval of a record is below, within or above the temperatures the
# model was fitted to. Within them the data hold the failure-rate function in
# place, so functions that fit equally well predict alike. Outside them nothing
# does, and the prediction depends on which function was chosen.
.record_zone <- function(intervals, Ta_range) {
  zone <- ifelse(intervals$Ta < Ta_range[1], "below",
          ifelse(intervals$Ta > Ta_range[2], "above", "within"))
  factor(zone, levels = c("below", "within", "above"))
}

# How much of a record lies below, within and above the temperatures the model was
# fitted to, in two senses: the share of its duration (`below`, `within`, `above`),
# and the share of the hazard accrued along it (`risk_below`, `risk_within`,
# `risk_above`, each the median over the posterior draws with any hazard, and
# missing when none has any). The two can differ greatly: a record that spends most
# of its time at cool temperatures, where little hazard accrues, may take nearly all
# of its risk from a few warm hours.
.record_coverage <- function(accrued, zone, duration, Ta_range) {
  zones <- levels(zone)
  time_share <- vapply(zones, \(z) sum(duration[zone == z]) / sum(duration), numeric(1))
  by_zone <- vapply(zones, \(z) rowSums(accrued[, zone == z, drop = FALSE]),
                    numeric(nrow(accrued)))
  by_zone <- matrix(by_zone, nrow = nrow(accrued))
  total   <- rowSums(by_zone)
  at_risk <- total > 0
  risk_share <- if (any(at_risk))
    apply(by_zone[at_risk, , drop = FALSE] / total[at_risk], 2, stats::median)
  else rep(NA_real_, length(zones))
  tibble::tibble(Ta_fitted_min = Ta_range[1], Ta_fitted_max = Ta_range[2],
                 below = time_share[["below"]], within = time_share[["within"]],
                 above = time_share[["above"]],
                 risk_below = risk_share[1], risk_within = risk_share[2],
                 risk_above = risk_share[3])
}


# ===========================================================================
# 2. Survival along a record, from a fitted model
# ===========================================================================

# .record_draws() / .pool_record() / .summarise_record() -----------------------

# Compute, for a fitted model, the posterior draws of the hazard accrued in each
# interval of a record, at the covariate values in each row of `newdata`. It checks
# its arguments, and says once when exposure has to be counted from the start of
# the record and when the record leaves the temperatures the model was fitted to.
# predict_survival() summarises these draws, and average_survival() pools them
# across models first.
.record_draws <- function(fit, record, newdata, levels, by, n_draws, grid_dT = 0.1) {
  st <- .fit_settings(fit)
  cs <- st$costressors
  cs_columns  <- vapply(cs, `[[`, character(1), "column")
  cs_names    <- vapply(cs, `[[`, character(1), "name")
  in_record   <- cs_columns %in% names(record)
  if (!is.null(by) && (!is.character(by) || length(by) != 1L || !by %in% names(record)))
    cli::cli_abort("{.arg by} must name one column of {.arg record}, the period each reading belongs to.",
                   class = "rTFT_bad_by")
  .check_record(record, c("time", "temp", cs_columns[in_record], by))
  if (is.null(newdata)) newdata <- data.frame(row.names = 1L)
  if (!is.numeric(grid_dT) || length(grid_dT) != 1L || is.na(grid_dT) || grid_dT <= 0)
    cli::cli_abort(c(
      "{.arg grid_dT} must be a single positive number, the largest change in temperature within a step along the record.",
      "i" = "Use {.code Inf} to take one step between consecutive readings."
    ), class = "rTFT_bad_grid_dT")
  intervals <- .record_intervals(record, cs_columns[in_record])
  fine      <- .refine_intervals(record, cs_columns[in_record], grid_dT)

  # A co-stressor the record follows needs no level; it is given a placeholder so
  # the levels can be checked as a whole, and its value in each interval is used.
  if (any(in_record))
    levels <- c(levels[setdiff(names(levels), cs_names[in_record])],
                stats::setNames(rep(0, sum(in_record)), cs_names[in_record]))

  pars <- c(M_catalogue[[st$m_form]]$dpars, h0_catalogue[[st$baseline]]$dpars,
            .costressor_dpars(cs))
  dr   <- .fit_draws(fit, pars, newdata, draw_ids = .thin_ids(brms::ndraws(fit), n_draws))

  # Check the co-stressor levels once, warning once about any the model does not
  # have; after that only the levels in use are passed on.
  k <- if (M_catalogue[[st$m_form]]$k_sw != 0) st$k_sw
  .rate_setup(st$m_form, purrr::map(dr, \(m) m[, 1]), st$Tr, k, cs, levels)
  levels <- if (length(cs)) levels[cs_names]

  # The hazard is accrued over the parts and added up within each interval between
  # readings, which is what survival and the shares of the risk are reported for.
  rows <- purrr::map(seq_len(nrow(newdata)), \(r) {
    accrued <- .record_hazard(st, purrr::map(dr, \(m) m[, r]), fine, levels)
    unname(t(rowsum(t(accrued), fine$parent, reorder = TRUE)))
  })
  list(newdata = newdata, time = record$time,
       record = tibble::as_tibble(record[c("time", "temp", cs_columns[in_record])]),
       period = if (!is.null(by)) record[[by]][-nrow(record)], by = by,
       zone = .record_zone(intervals, st$Ta_range),
       duration = intervals$tau_end - intervals$tau_start, Ta_range = st$Ta_range,
       baseline = st$baseline, rows = rows)
}

# Pool the draws of several models, as .record_draws() returns them, into one set in
# the same layout: for each row of newdata, the models' draws one after another.
.pool_record <- function(per_model) {
  first <- per_model[[1]]
  first$rows <- purrr::map(seq_along(first$rows), \(r)
    do.call(rbind, purrr::map(per_model, \(m) m$rows[[r]])))
  first$baseline <- unique(vapply(per_model, \(m) m$baseline, character(1)))
  first
}

# Summarise the draws from .record_draws() into the tables predict_survival()
# returns, each carrying the columns of newdata: `survival` with a row per reading of
# the record, `mortality` with one row, when the record names periods `risk` with a
# row per period, and `coverage` with one row.
.summarise_record <- function(draws, prob) {
  per_row <- purrr::map(seq_along(draws$rows), \(r) {
    accrued <- draws$rows[[r]]
    row     <- draws$newdata[r, , drop = FALSE]
    # Survival at each reading: 1 at the first, then exp(-hazard accrued so far).
    cumulative <- t(apply(accrued, 1, cumsum))
    if (ncol(accrued) == 1L) cumulative <- matrix(cumulative, ncol = 1L)
    survival   <- cbind(1, exp(-cumulative))
    out <- list(
      survival  = tibble::tibble(time = draws$time, .summarise_columns(survival, prob)),
      mortality = tibble::tibble(time = draws$time[length(draws$time)],
                                 .draw_summary(1 - survival[, ncol(survival)], prob)))
    if (!is.null(draws$by)) out$risk <- .risk_by_period(accrued, draws$period, draws$by, prob)
    out$coverage <- .record_coverage(accrued, draws$zone, draws$duration, draws$Ta_range)
    purrr::map(out, \(tbl) .with_newdata(row, tbl))
  })
  tables <- names(per_row[[1]])
  purrr::map(stats::setNames(tables, tables),
             \(tbl) dplyr::bind_rows(purrr::map(per_row, tbl)))
}

# Each period's share of the hazard accrued over the whole record. The shares of a
# draw add to one. A draw with no hazard anywhere has no shares to give, so it is
# left out, and `n_draws` counts the draws the summary rests on; with none, the
# shares are missing.
.risk_by_period <- function(accrued, period, by, prob) {
  periods <- unique(period)
  by_period <- vapply(periods, \(p) rowSums(accrued[, period == p, drop = FALSE]),
                      numeric(nrow(accrued)))
  by_period <- matrix(by_period, nrow = nrow(accrued))
  total <- rowSums(by_period)
  at_risk <- total > 0
  labels <- stats::setNames(tibble::tibble(periods), by)
  if (!any(at_risk))
    return(dplyr::bind_cols(labels, tibble::tibble(
      median = NA_real_, lower = NA_real_, upper = NA_real_, prob = prob, n_draws = 0L)))
  share <- .summarise_columns(by_period[at_risk, , drop = FALSE] / total[at_risk], prob)
  dplyr::bind_cols(labels, share[c("median", "lower", "upper", "prob", "n_draws")])
}

# .report_record() -------------------------------------------------------------

# Say, once, what a prediction along a record rests on beyond the fitted model:
# where exposure is counted from, when the baseline depends on elapsed exposure, and
# how much of the record, and of the risk predicted along it, lies outside the
# temperatures the model was fitted to. `cover` is the coverage table, a row per row
# of newdata; where the share of the risk differs among its rows, the largest is
# given.
.report_record <- function(cover, baseline) {
  first   <- cover[1, ]
  outside <- 1 - cover$risk_within
  lines <- c(
    if (any(baseline != "mem"))
      c("i" = "Elapsed exposure is counted from the first time in {.arg record}, since the baseline hazard depends on it."),
    if (first$within < 1)
      c("!" = sprintf("%.0f%% of the record, by duration, lies outside the temperatures the model was fitted to (%.4g to %.4g): %.0f%% below, %.0f%% above.",
                      100 * (1 - first$within), first$Ta_fitted_min, first$Ta_fitted_max,
                      100 * first$below, 100 * first$above),
        if (!all(is.na(outside)))
          c("!" = sprintf("About %.0f%% of the predicted risk is accrued at those temperatures.",
                          100 * max(outside, na.rm = TRUE))),
        "i" = "There the data do not constrain the failure-rate function, so the prediction depends on which function was chosen."))
  if (length(lines)) cli::cli_inform(lines, class = "rTFT_record_notes")
  invisible(NULL)
}

# predict_survival() -----------------------------------------------------------

#' Survival along a record of temperature, from a fitted failure-rate model
#'
#' Predicts, from a model fitted with [fit_hazard()], the survival of an individual
#' exposed to a measured record of temperature: the probability of surviving to
#' each time of the record and of failing by its end, each with its posterior
#' median and credible interval, at the covariate values in each row of `newdata`.
#' Random effects are set to zero, so the prediction is for a typical group.
#'
#' Along the record the cumulative hazard grows, over each step, by
#' M(T) (H0(tau_1) - H0(tau_0)): the failure rate at the temperature in the middle of
#' the step, times the growth of the cumulative baseline hazard over it. A step is
#' the interval between two readings, or a part of it where the temperature changes
#' by more than `grid_dT` between them. Survival to a time is exp(-cumulative hazard up to it). With the
#' memoryless baseline the growth of H0 is the length of the interval, so the
#' cumulative hazard is the failure rate summed over the record, and it does not
#' matter where the record begins. With the other baselines the hazard depends on
#' elapsed exposure, which is counted from the first time in the record.
#'
#' @section Outside the temperatures of the assay:
#' An assay is run at temperatures that cause failure within its duration, and a
#' field record may lie below them. Within the temperatures of the assay the data
#' hold the failure-rate function in place, so functions that fit equally well
#' predict alike. Below them nothing does, and the prediction depends on which
#' function was chosen: functions that fit the assay equally well can differ by
#' orders of magnitude there. The `coverage` table gives the share of the record
#' below, within and above the temperatures the model was fitted to, by duration
#' and by the risk accrued. The second is the one to read: a record may spend most
#' of its time at cool temperatures where little hazard accrues, and take nearly
#' all of its risk from a few warm hours that the assay did cover. A message is
#' given when any of the record lies outside. Comparing the predictions
#' of several failure-rate functions, or averaging them with [average_survival()],
#' shows how much the prediction depends on that choice.
#'
#' The zero-type failure rates (`zeroexp`, `zeropower`) are zero at the reference
#' temperature `Tr` and are taken as zero below it. The `power` function is not
#' defined below `Tr`, and a record that goes below it is refused.
#'
#' @section No recovery:
#' The hazard accrued along the record is only ever added to, so the prediction
#' assumes that an individual does not recover when the temperature returns to a
#' range in which repair is expected to exceed thermal damage. Along a record that
#' moves in and out of stressful temperatures, such as one with a daily cycle, the
#' predicted mortality is then too high if individuals recover in between, and the
#' more so the more cycles there are.
#'
#' @section Co-stressors:
#' For a model with co-stressors, a co-stressor whose column is in `record` is
#' followed through time, taking in each interval the value halfway between the two
#' readings. One that is not in `record` is held at the level given in `levels`.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param record A data frame with the columns `time` and `temp`, the record of
#'   temperature through time, in the time unit and temperature scale of the data
#'   the model was fitted to, with times increasing and no gaps. It may also hold
#'   the column of a co-stressor and the column named in `by`.
#' @param newdata A data frame with one row per set of covariate values to predict
#'   at, holding every covariate in the parameters' formulas. `NULL` (the default)
#'   suits a model with no covariates.
#' @param levels A named numeric vector giving the level of each co-stressor that
#'   `record` does not follow, on its measurement scale, e.g. `c(DO = 6)`.
#' @param by A character string, the name of a column of `record` giving the period
#'   each reading belongs to (e.g., a day). When given, the result has each period's
#'   share of the hazard accrued over the record. An interval between two readings
#'   belongs to the period of the first.
#' @param prob A single number between 0 and 1, the coverage of the credible
#'   intervals. Defaults to 0.95.
#' @param n_draws A single number, the number of posterior draws to use, evenly
#'   spaced through the posterior. Defaults to 1000. The draws are held for every
#'   step along the record, so a long record with many draws needs memory in
#'   proportion to their product.
#' @param grid_dT A single positive number, the largest change in temperature within
#'   a step along the record, in degrees. Between two readings the temperature is
#'   taken to change linearly, and the hazard is accrued in steps over which it
#'   changes by at most this much, since the failure rate at the midpoint of a step
#'   stands for the whole of it (see *The time grid* in [build_hazard_data]).
#'   Defaults to 0.1; `Inf` takes one step between consecutive readings.
#' @return A list of tibbles, each with the columns of `newdata`:
#'   * `survival`: for each time of the record, the posterior `median`, `lower` and
#'     `upper` of the probability of surviving to it (1 at the first time).
#'   * `mortality`: the probability of failing by the end of the record.
#'   * `risk` (when `by` is given): each period's share of the hazard accrued over
#'     the record. The shares of a posterior draw add to one; `n_draws` is the
#'     number of draws with any hazard, which the summary rests on.
#'   * `coverage`: the temperatures the model was fitted to; the share of the
#'     record's duration `below`, `within` and `above` them; and the share of the
#'     hazard accrued along the record at those temperatures, `risk_below`,
#'     `risk_within` and `risk_above` (the median over the posterior draws with any
#'     hazard).
#'   * `record`: the times and temperatures of the record, with any co-stressor it
#'     follows, kept so that [plot_survival()] can draw them above the survival
#'     they led to.
#' @examples
#' \dontrun{
#' # A two-day record read every hour, with a daily cycle.
#' hours  <- 0:48
#' record <- data.frame(time = hours * 60,
#'                      temp = 24 + 3 * sin(2 * pi * (hours - 9) / 24),
#'                      day  = pmin(hours %/% 24, 1))
#' pred <- predict_survival(fit, record, newdata = data.frame(mass = 10), by = "day")
#' pred$mortality
#' pred$risk
#' }
#' @export
predict_survival <- function(fit, record, newdata = NULL, levels = NULL, by = NULL,
                             prob = 0.95, n_draws = 1000, grid_dT = 0.1) {
  .check_prob(prob)
  draws <- .record_draws(fit, record, newdata, levels, by, n_draws, grid_dT)
  out   <- .summarise_record(draws, prob)
  .report_record(out$coverage, draws$baseline)
  c(out, list(record = draws$record))
}


# ===========================================================================
# 3. Survival along a record, averaged over models
# ===========================================================================

# average_survival() -----------------------------------------------------------

#' Survival along a record of temperature, averaged over several fitted models
#'
#' Predicts survival along a record as [predict_survival()] does, averaged over
#' several models fitted to the same data. Each model's posterior draws are pooled
#' in proportion to its weight and the pool is summarised, so the credible intervals
#' reflect both the uncertainty within each model and the uncertainty about which
#' model holds. This matters most where the record lies outside the temperatures of
#' the assay, where failure-rate functions that fit the data equally well can
#' predict very differently (see [predict_survival()]).
#'
#' The weights say how well each model predicts failures within the assay. They do
#' not say which model extrapolates better below it, which the data cannot tell.
#'
#' @inheritParams average_tdt_ttl
#' @inheritParams predict_survival
#' @return The list [predict_survival()] returns, with `weights`, the table of
#'   weights with the number of draws each model contributed. In `coverage`, the
#'   shares of the risk are those of the pooled draws.
#' @examples
#' \dontrun{
#' fits <- list("exp" = fit_exp, "zeroexp" = fit_zeroexp)
#' w    <- hazard_weights(lapply(fits, loo_hazard))
#' average_survival(fits, w, record, newdata = data.frame(mass = 10))$mortality
#' }
#' @export
average_survival <- function(fits, weights, record, newdata = NULL, levels = NULL,
                             by = NULL, prob = 0.95, n_draws = 1000, grid_dT = 0.1) {
  .check_prob(prob)
  setup <- .average_setup(fits, weights, n_draws)
  per_model <- .each_model(setup, \(fit, n)
    .record_draws(fit, record, newdata, levels, by, n, grid_dT))
  draws <- .pool_record(per_model)
  out   <- .summarise_record(draws, prob)
  .report_record(out$coverage, draws$baseline)
  c(out, list(record = draws$record, weights = setup$weights))
}
