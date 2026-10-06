# derive-tdt-ttl.R

# Functions to derive thermal tolerance from a fitted failure-rate model, with the
# posterior uncertainty carried through. Nothing here is fitted: survival to an
# elapsed exposure tau at a constant temperature T is S(tau, T) = exp(-M(T) H0(tau)),
# and every quantity reported is read off that expression, one posterior draw at a
# time. The failure rate M(T) is evaluated from the same log M(T) expression the
# model was fitted with, and the cumulative baseline hazard H0 and its inverse have
# closed forms. The quantities are the model's parameters on their natural scale,
# the static critical thermal maximum sCTmax (the constant temperature at which a
# given fraction fail within a reference time), the thermal death time (TDT) curves
# (the time until a given fraction have failed, against temperature) and the
# thermal tolerance landscape (survival over temperature and elapsed time).
#
# The path is derive_tdt_ttl(), which calls .derive_draws() to compute the posterior
# draws of every quantity and .summarise_derived() to reduce them to a median and a
# credible interval. The sections below build up to it: the baseline cumulative
# hazard, the failure rate for each draw, the summaries, sCTmax by bisection, and
# the reading of draws out of a fitted model at chosen covariate values. A draw
# that never reaches a quantity counts as infinite and is kept in the summaries.
# M_draws() returns the fitted failure rate itself, and H0_of() and H0inv_of() the
# baseline cumulative hazard and its inverse.
#
# "Thermal death time" is the literature's name for that curve, kept here so the
# package connects to it, whatever the failure endpoint: loss of equilibrium,
# knockdown or death. Everywhere else the package speaks of failure.

# ===========================================================================
# 1. Baseline cumulative hazard H0(tau) and its inverse
# ===========================================================================

# The closed-form cumulative baseline hazard for each function of h0_catalogue, and
# its inverse:
#   mem       h0 = 1                    H0(tau) = tau                 H0^-1(x) = x
#   weibull   h0 = (k/s)(tau/s)^(k-1)   H0(tau) = (tau/s)^k           H0^-1(x) = s x^(1/k)
#   llogis    log-logistic hazard       H0(tau) = log(1 + (tau/s)^k)  H0^-1(x) = s (e^x - 1)^(1/k)
#   gompertz  h0 = e^(g tau)            H0(tau) = (e^(g tau) - 1)/g    H0^-1(x) = log(1 + g x)/g
#
# Shape (k) and scale (s) arrive as per-draw vectors; tau or x may be a vector or a
# matrix whose rows are draws, and recycle column-wise over draws. The Weibull scale
# is not estimated: it is the reference time tau_ref, the same for every draw, so
# its H0 is 1 at tau_ref. The Gompertz shape
# holds its rate constant gamma in the first slot and ignores the second, because
# its amplitude is loga in M(T) rather than a parameter of its own.
#
# Both Gompertz expressions are singular at gamma = 0, where they are 0/0 rather than
# undefined: the limits are H0(tau) = tau and H0^-1(x) = x, the memoryless case. gamma
# passes through zero freely during sampling, so near gamma = 0 the functions below
# switch to the first term of the series (tau * (1 + u/2) and x * (1 - u/2)) rather
# than let the division produce NaN at the one value of the parameter the shape is
# meant to generalise. .gz_tol sets how close to zero counts as near.
.gz_tol <- 1e-8

# Make sure each baseline returns one value per posterior draw. Shape and scale come
# with one value per draw, and the Weibull, log-logistic and Gompertz baselines
# combine them with elapsed time, so their results already have one value per draw.
# The memoryless baseline, H0(tau) = tau, uses neither shape nor scale, so a single
# elapsed time would come back as a single value. That value is right for every draw,
# but later calculations read the result one draw at a time, taking draw i's value
# from position i, and would find nothing after the first. This helper repeats the
# shorter input until it is as long as the longest one. Elapsed time given as a
# matrix, with one row per draw, is already the longest input and comes back
# unchanged.
.recycle_to <- function(x, ...) {
  n <- max(length(x), vapply(list(...), length, integer(1)))
  if (length(x) == n) x else rep_len(x, n)
}


# H0_of() ----------------------------------------------------------------------

#' The cumulative baseline hazard of a given shape
#'
#' Returns the function H0(tau) for the chosen baseline. Under the memoryless shape
#' it is elapsed exposure itself; the Weibull, log-logistic and Gompertz shapes let
#' the hazard rise or fall with how long the individual has already been exposed. Shape
#' and scale arrive as per-draw vectors, and elapsed time may be a vector or a matrix
#' whose rows are draws.
#'
#' For the Gompertz shape the first argument is the rate constant gamma and the
#' second is ignored. A negative gamma is a valid value and means a hazard that
#' falls the longer the individual has been exposed. The cumulative hazard then levels
#' off: it rises towards 1/|gamma| but never reaches it, instead of growing without
#' limit. Survival therefore never drops to zero, and the model predicts that some
#' individuals never fail (a cure fraction, in survival-analysis terms). [H0inv_of()]
#' gives the time by which a chosen fraction of individuals has failed; when that
#' fraction is never reached it returns `NaN` for the time, rather than a large
#' number standing in for "never".
#'
#' @param baseline A character string, a key of [h0_catalogue]: `"mem"`, `"weibull"`,
#'   `"llogis"` or `"gompertz"`.
#' @return A function of `(tau, shape, scale)` giving the cumulative baseline hazard.
#'   What `shape` and `scale` mean depends on the baseline:
#'   * `"mem"`: both are ignored.
#'   * `"weibull"` and `"llogis"`: `shape` is k and `scale` is s.
#'   * `"gompertz"`: `shape` is the rate constant gamma and `scale` is ignored.
#' @examples
#' # The memoryless baseline is elapsed exposure itself.
#' H0_of("mem")(tau = c(10, 30, 60), shape = 1, scale = 1)
#'
#' # A Weibull baseline with shape 2 and scale 60: H0(tau) = (tau / 60)^2.
#' H0_of("weibull")(tau = 30, shape = 2, scale = 60)
#' @export
H0_of <- function(baseline) {
  .baseline_spec(baseline)   # errors rTFT_unknown_baseline on a bad key
  switch(baseline,
    mem      = function(tau, shape, scale) .recycle_to(tau, shape, scale),
    weibull  = function(tau, shape, scale) (tau / scale)^shape,
    llogis   = function(tau, shape, scale) log1p((tau / scale)^shape),
    gompertz = function(tau, shape, scale) {
      u <- shape * tau
      ifelse(abs(u) < .gz_tol, tau * (1 + u / 2), expm1(u) / shape)
    })
}


# H0inv_of() -------------------------------------------------------------------

#' The inverse of the cumulative baseline hazard
#'
#' Returns the function that gives the elapsed time at which the cumulative baseline
#' hazard reaches a given value, which is what turns a failure probability into a
#' thermal death time.
#'
#' Under the Gompertz shape with a negative rate constant the cumulative hazard
#' levels off below 1/|gamma|, so a target above that level is never reached. The
#' time returned there is `NaN`, rather than a large number standing in for
#' "never". A target exactly at the level is reached only after infinite time, and
#' returns `Inf`.
#'
#' @param baseline A character string, a key of [h0_catalogue]: `"mem"`, `"weibull"`,
#'   `"llogis"` or `"gompertz"`.
#' @return A function of `(x, shape, scale)` giving the elapsed time at which the
#'   cumulative baseline hazard equals `x`. `shape` and `scale` mean the same as in
#'   [H0_of()]:
#'   * `"mem"`: both are ignored.
#'   * `"weibull"` and `"llogis"`: `shape` is k and `scale` is s.
#'   * `"gompertz"`: `shape` is the rate constant gamma and `scale` is ignored.
#' @examples
#' # Invert a Weibull baseline (shape 2, scale 60): the time at which H0 reaches 0.25.
#' H0inv_of("weibull")(x = 0.25, shape = 2, scale = 60)   # 30
#' @export
H0inv_of <- function(baseline) {
  .baseline_spec(baseline)   # errors rTFT_unknown_baseline on a bad key
  switch(baseline,
    mem      = function(x, shape, scale) .recycle_to(x, shape, scale),
    weibull  = function(x, shape, scale) scale * x^(1 / shape),
    llogis   = function(x, shape, scale) scale * expm1(x)^(1 / shape),
    gompertz = function(x, shape, scale) {
      u <- shape * x
      # With a negative rate constant the cumulative hazard levels off at 1/|gamma|.
      # Because gamma is negative, u = gamma * x runs opposite to the target x: a
      # target below the level has u between -1 and 0, a target at the level has
      # u = -1, and a target above it has u < -1. Below the level the target is
      # reached at a finite time, which grows without limit as the target nears the
      # level; at the level, log1p(-1) / gamma gives Inf; above it the target is
      # never reached and the time is NaN. ifelse() computes both of its
      # branches for every target, so log1p() would still meet the targets above
      # the level and warn "NaNs produced", which reads as a numerical fault. The
      # inner ifelse() hands it 0 for those targets instead, and the outer one puts
      # NaN in their place.
      val <- ifelse(u < -1, NaN, log1p(ifelse(u < -1, 0, u)) / shape)
      ifelse(abs(u) < .gz_tol, x * (1 - u / 2), val)
    })
}


# ===========================================================================
# 2. The failure rate M(T) for each posterior draw
# ===========================================================================

# .rate_setup() / .log_M_eval() ------------------------------------------------

# Check and prepare everything the failure rate needs apart from the temperatures:
# the function, the draws, the reference temperature, the softening steepness and
# the co-stressors at their chosen levels. Returns the pieces .log_M_eval() needs.
# .M_draws() and .sctmax_draws() share it, so both refuse the same bad inputs the
# same way.
.rate_setup <- function(m_form, draws, Tr, k_sw = NULL, costressors = list(),
                        levels = NULL) {
  m_spec <- .m_form_spec(m_form)
  .check_costressors(costressors)
  k_sw_used <- .resolve_k_sw(m_form, k_sw)

  if (!is.numeric(Tr) || length(Tr) != 1L || !is.finite(Tr))
    cli::cli_abort(c(
      "{.arg Tr} must be a single number, the reference temperature of the fit.",
      "x" = "You supplied {.val {Tr}}."
    ), class = "rTFT_bad_Tr")

  # Every parameter the rate needs must have draws, and the same number of them. A
  # missing one would otherwise fail later with an unhelpful "object not found", and
  # a short one would be silently reused across the draws, shrinking the spread of
  # the result.
  pars    <- c(m_spec$dpars, .costressor_dpars(costressors))
  draws   <- as.list(draws)
  missing <- setdiff(pars, names(draws))
  if (length(missing))
    cli::cli_abort(c(
      "{.arg draws} has no draws for {.val {missing}}.",
      "i" = "{.val {m_form}} with these co-stressors needs {.val {pars}}."
    ), class = "rTFT_missing_draws")
  n_draws <- lengths(draws[pars])
  if (length(unique(n_draws)) != 1L)
    cli::cli_abort(c(
      "Every parameter in {.arg draws} must have the same number of draws.",
      "x" = "Lengths: {paste(names(n_draws), n_draws, sep = ' = ', collapse = ', ')}."
    ), class = "rTFT_bad_draws")

  # Each co-stressor enters at its chosen level, standardised exactly as the fit
  # standardised its data.
  cs_names <- vapply(costressors, `[[`, character(1), "name")
  if (length(costressors)) {
    if (!is.numeric(levels) || is.null(names(levels)) || anyNA(levels))
      cli::cli_abort(c(
        "{.arg levels} must be a named numeric vector, one level per co-stressor.",
        "i" = "For example {.code levels = c({cs_names[1]} = 6)}."
      ), class = "rTFT_bad_levels")
    no_level <- setdiff(cs_names, names(levels))
    if (length(no_level))
      cli::cli_abort(c(
        "{.arg levels} gives no level for {.val {no_level}}.",
        "i" = "Give one for each co-stressor: {.val {cs_names}}."
      ), class = "rTFT_bad_levels")
  }
  # A level for something the model does not have is ignored, with a warning, since
  # it usually means a misspelled name or the wrong model.
  if (!is.null(levels)) {
    if (!length(costressors))
      cli::cli_warn(c(
        "!" = "Ignoring {.arg levels}: the model has no co-stressors."
      ), class = "rTFT_unused_levels")
    else {
      extra <- setdiff(names(levels), cs_names)
      if (length(extra))
        cli::cli_warn(c(
          "!" = "Ignoring {.arg levels} with no matching co-stressor: {.val {extra}}.",
          "i" = "The model's co-stressors are {.val {cs_names}}."
        ), class = "rTFT_unused_levels")
    }
  }
  z <- purrr::map(costressors, \(cs) .costressor_z(cs, levels[[cs$name]]))
  names(z) <- if (length(z)) paste0("z_", cs_names)

  list(draws = draws[pars], nd = n_draws[[1]], Tr = Tr, z = z,
       log_M = str2lang(.log_M_string(m_form, k_sw_used, costressors)))
}

# Evaluate log M(T) for every draw at the temperatures in `Ta_mat`, a matrix with
# one row per draw: row i holds the temperatures at which draw i is evaluated. Each
# parameter is laid out the same shape, one value per draw repeated along its row,
# so the expression gives log M for every cell at once. The expression sees only
# these inputs and base R.
.log_M_eval <- function(setup, Ta_mat) {
  nd <- setup$nd
  nT <- ncol(Ta_mat)
  env <- c(purrr::map(setup$draws, \(v) matrix(v, nd, nT)),
           list(Ta = Ta_mat, Tr = setup$Tr),
           setup$z)
  eval(setup$log_M, env, enclos = baseenv())
}


# .M_draws() -------------------------------------------------------------------

# Evaluate the failure rate M(T) at each temperature in `Ta` for each posterior
# draw, returning a matrix with one row per draw and one column per temperature.
# It evaluates the same log M(T) expression the model was fitted with, including
# the softening steepness and any co-stressor terms, so the rate reported is the
# rate that was fitted. With co-stressors the rate is M(T, Z), evaluated at a chosen
# level of each.
#
# The settings (`m_form`, `Tr`, `k_sw`, `costressors`) must be the ones the model
# was fitted with, so this is kept internal: the reporting functions read them from
# the fit and pass them on (a user never re-enters them to avoid errors). The inputs:
#   m_form       a key of M_catalogue
#   draws        a named list or data frame, one vector of draws per parameter; it
#                must hold every parameter of m_form and of the co-stressors
#   Ta           the temperatures at which to evaluate the rate
#   Tr           the reference temperature of the fit
#   k_sw         the softening steepness set at fitting, or NULL for the catalogue's
#   costressors  the list of costressor() objects the model was fitted with
#   levels       a named vector giving each co-stressor's level on its measurement
#                scale, e.g. c(DO = 6) for 6 mg/L of dissolved oxygen
.M_draws <- function(m_form, draws, Ta, Tr, k_sw = NULL, costressors = list(),
                     levels = NULL) {
  setup <- .rate_setup(m_form, draws, Tr, k_sw, costressors, levels)

  if (!is.numeric(Ta) || !length(Ta) || anyNA(Ta))
    cli::cli_abort("{.arg Ta} must be a numeric vector of temperatures with no missing values.",
                   class = "rTFT_bad_Ta")
  # zeroexp, zeropower and power take a logarithm or a fractional power of Ta - Tr,
  # so they have no rate below the reference, and a grid reaching below it is
  # refused, as the fit refuses data that do. Ta = Tr itself is allowed: the zero
  # functions give a rate of exactly zero there, which is how they are defined.
  if (m_form %in% c("zeroexp", "zeropower", "power") && any(Ta < Tr))
    cli::cli_abort(c(
      "The {.val {m_form}} failure rate is not defined below the reference {.arg Tr}.",
      "x" = "min(Ta) = {.val {min(Ta)}}, Tr = {.val {Tr}}.",
      "i" = "Start the temperature grid at or above {.arg Tr}."
    ), class = "rTFT_Ta_below_Tr")

  # Every draw at every temperature: each temperature is repeated down the draws.
  exp(.log_M_eval(setup, matrix(Ta, setup$nd, length(Ta), byrow = TRUE)))
}


# ===========================================================================
# 3. Posterior summaries
# ===========================================================================

# .check_prob() ----------------------------------------------------------------

# Check an interval's coverage: a single number strictly between 0 and 1.
.check_prob <- function(prob) {
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1)
    cli::cli_abort(c(
      "{.arg prob} must be a single number between 0 and 1, the interval's coverage.",
      "x" = "You supplied {.val {prob}}."
    ), class = "rTFT_bad_prob")
  invisible(prob)
}


# .summarise_columns() ---------------------------------------------------------

# The posterior summary described for .draw_summary() below, applied to every column
# of a draws-by-columns matrix in one pass, one row per column: the thermal death
# time curves and the tolerance landscape have one column per temperature or per
# cell. .draw_summary() is its one-column case, so both summarise alike.
.summarise_columns <- function(mat, prob = 0.95) {
  .check_prob(prob)
  if (!is.numeric(mat) || !length(mat))
    cli::cli_abort("The draws must be a non-empty numeric vector or matrix.",
                   class = "rTFT_bad_summary_draws")
  mat <- as.matrix(mat)
  if (anyNA(mat))
    cli::cli_abort(c(
      "{sum(is.na(mat))} of {length(mat)} draws {?is/are} missing (NA or NaN).",
      "i" = "A quantity that is never reached should be passed as {.code Inf} or {.code -Inf}."
    ), class = "rTFT_missing_draw_values")

  probs <- c(0.5, (1 - prob) / 2, 1 - (1 - prob) / 2)
  q <- apply(mat, 2, stats::quantile, probs = probs, names = FALSE)
  # Interpolating between a -Inf draw and a +Inf draw gives NaN. Only there, take
  # the nearest draw instead.
  for (j in which(colSums(is.na(q)) > 0)) {
    bad <- is.na(q[, j])
    q[bad, j] <- stats::quantile(mat[, j], probs[bad], names = FALSE, type = 1)
  }

  tibble::tibble(median = q[1, ], lower = q[2, ], upper = q[3, ], prob = prob,
                 n_draws = nrow(mat),
                 n_unreached = as.integer(colSums(is.infinite(mat))))
}


# .draw_summary() --------------------------------------------------------------

# Summarise one quantity's posterior draws as a median and an equal-tailed credible
# interval of coverage `prob` (95% by default), with the number of draws and how
# many of them never reach the quantity. Every quantity the package reports is
# summarised here, so all are summarised alike.
#
# A draw for which the quantity is never reached (a thermal death time under a
# cure fraction, a CTmax the rate cannot attain) arrives as +Inf or -Inf and is
# kept. Dropping it would summarise only the draws where the quantity was reached,
# pulling the result towards earlier times and lower temperatures than the
# posterior implies. Kept, it enters the quantiles honestly: while fewer than
# (1 - prob) / 2 of the draws never reach the quantity (2.5% for a 95% interval),
# both bounds stay finite; above that share the outer bound is infinite, and from
# half onwards the median is. `n_unreached` counts these draws.
#
# Since a quantity that is not reached must arrive as +Inf or -Inf, a missing value
# (NA or NaN) can only mean a fault, and stops here rather than being dropped.
.draw_summary <- function(x, prob = 0.95) {
  .summarise_columns(matrix(x, ncol = 1L), prob)
}


# ===========================================================================
# 4. The static critical thermal maximum for each posterior draw
# ===========================================================================

# .sctmax_draws() --------------------------------------------------------------

# Find, for each posterior draw, the static critical thermal maximum sCTmax: the
# constant temperature at which a fraction p of individuals fail within a reference
# exposure time t_ref (Jorgensen et al. 2021 report it for one hour, sCTmax(1h)).
# Under a constant temperature survival is exp(-M(T) * H0(t_ref)), so sCTmax solves
# M(T) = target, with target = -log(1 - p) / H0(t_ref). The caller computes the
# target, one per draw.
#
# The equation is solved on the shared log M(T) expression, by bisection run on
# every draw at once, rather than read off a grid of temperatures. The answer then
# does not depend on a grid, and there is no second, hand-written copy of each
# function to drift from the fit. Every function of M_catalogue rises with
# temperature over its domain (or, for `const`, stays flat), so a draw has at most
# one crossing.
#
# The search starts at the lower end of the function's domain: Tr for zeroexp,
# zeropower and power, which are undefined below it; Tc for zeroexp_threshold and
# zeropower_threshold, whose softened switch makes a small spurious rise just below
# Tc, so that the crossing found is on the branch the function means; otherwise
# just below Tr. It widens in doubling steps until the crossing is bracketed, or
# until it reaches .sctmax_span either side of Tr. A draw with no crossing is
# returned as +Inf when its rate stays below the target throughout (a constant rate
# below it, a threshold rate whose ceiling a lies below it, a crossing beyond the
# span) and as -Inf when its rate is above the target throughout its domain (a
# constant rate above it, a power rate whose floor c lies above it).
#
# The one way a rate can fall with temperature is through a co-stressor's
# temperature interaction, when it outweighs the rate's own rise at the chosen
# level. sCTmax is not defined there, and a search that assumes a rising rate
# would report such a draw as never reaching the target, which is wrong. So a draw
# whose rate is lower at the top of its search range than at the bottom stops the
# calculation with an error; such a level usually lies outside the range the data
# support.
.sctmax_span <- 1000
.sctmax_tol  <- 1e-8

.sctmax_draws <- function(m_form, draws, target, Tr, k_sw = NULL,
                          costressors = list(), levels = NULL) {
  setup <- .rate_setup(m_form, draws, Tr, k_sw, costressors, levels)
  nd <- setup$nd
  if (!is.numeric(target) || length(target) != nd || !all(is.finite(target)) ||
      any(target <= 0))
    cli::cli_abort(c(
      "{.arg target} must hold one positive, finite failure rate per draw.",
      "x" = "There are {nd} draws and {length(target)} target{?s}, of which {sum(!is.finite(target) | target <= 0)} {?is/are} not positive and finite."
    ), class = "rTFT_bad_target")

  log_target <- log(target)
  gap <- function(Tv) as.vector(.log_M_eval(setup, matrix(Tv, nd, 1L))) - log_target

  # The lower end of each draw's domain, or NA where the function has none.
  domain_lo <- switch(m_form,
    zeroexp = , zeropower = , power           = rep(Tr, nd),
    zeroexp_threshold = , zeropower_threshold = setup$draws$Tc,
    rep(NA_real_, nd))
  bounded <- !is.na(domain_lo)
  lo <- ifelse(bounded, domain_lo, Tr - 1)
  hi <- pmax(lo, Tr) + 1

  # Widen until the rate is below the target at `lo` and above it at `hi`.
  gap_lo <- gap(lo)
  gap_hi <- gap(hi)
  step <- 1
  repeat {
    raise <- gap_hi < 0 & hi < Tr + .sctmax_span
    lower <- gap_lo > 0 & !bounded & lo > Tr - .sctmax_span
    if (!any(raise | lower)) break
    step <- 2 * step
    hi[raise] <- pmin(hi[raise] + step, Tr + .sctmax_span)
    lo[lower] <- pmax(lo[lower] - step, Tr - .sctmax_span)
    gap_lo <- gap(lo)
    gap_hi <- gap(hi)
  }
  falling <- gap_hi < gap_lo
  if (any(falling))
    cli::cli_abort(c(
      "sCTmax is not defined where the failure rate falls with temperature.",
      "x" = paste0("The rate falls with temperature for {sum(falling)} of {nd} draws",
                   if (length(costressors)) " at {.val {paste(names(levels), levels, sep = ' = ')}}", "."),
      "i" = "A temperature interaction strong enough to reverse the rate's rise usually means the co-stressor level lies outside the range the data support."
    ), class = "rTFT_falling_rate")
  out <- rep(NA_real_, nd)
  out[gap_hi < 0] <- Inf                                 # never reaches the target
  out[gap_lo > 0] <- -Inf                                # above it throughout

  # Halve each bracket until it is narrower than .sctmax_tol.
  open <- is.na(out)
  while (any(open)) {
    mid <- (lo + hi) / 2
    below <- gap(mid) < 0
    lo[open & below]  <- mid[open & below]
    hi[open & !below] <- mid[open & !below]
    open <- open & (hi - lo) > .sctmax_tol
  }
  ifelse(is.na(out), (lo + hi) / 2, out)
}


# ===========================================================================
# 5. Posterior draws from a fitted model
# ===========================================================================

# .thin_ids() / .complete_newdata() --------------------------------------------

# Pick `n_keep` evenly spaced draws out of `n_total`, always including the first and
# the last, or all of them when `n_keep` is NULL or not smaller. Evenly spaced rather
# than random, so the same call gives the same draws without a seed.
.thin_ids <- function(n_total, n_keep = NULL) {
  if (is.null(n_keep) || n_keep >= n_total) return(seq_len(n_total))
  if (!is.numeric(n_keep) || length(n_keep) != 1L || !is.finite(n_keep) || n_keep < 1)
    cli::cli_abort(c(
      "The number of draws to keep must be a single positive number.",
      "x" = "You supplied {.val {n_keep}}."
    ), class = "rTFT_bad_thin")
  unique(round(seq(1, n_total, length.out = n_keep)))
}

# Make `newdata` acceptable to brms, which checks it against the whole model and so
# asks for every column the model's formula uses, including those (the temperature,
# the reference, the log exposure, the failure indicator) that no parameter's own
# formula uses. Those are filled from the first row of the fitted data; their values
# cannot reach the parameters' draws. The covariates in the parameters' formulas are
# another matter: they set the values the draws are reported at, so a missing one is
# an error rather than a silent fill with some individual's value.
.complete_newdata <- function(fit, newdata, pars) {
  if (!is.data.frame(newdata) || !nrow(newdata))
    cli::cli_abort("{.arg newdata} must be a data frame with at least one row.",
                   class = "rTFT_bad_newdata")
  covariates <- unique(unlist(purrr::map(pars, \(p) {
    terms <- .rhs_terms(fit$formula$pforms[[p]])
    unlist(purrr::map(terms[!vapply(terms, .is_bar, logical(1))], all.vars))
  })))
  missing <- setdiff(covariates, names(newdata))
  if (length(missing))
    cli::cli_abort(c(
      "{.arg newdata} has no column for {.val {missing}}.",
      "i" = "The fitted model uses {cli::qty(length(missing))}{?it/them} in its parameters' formulas; give the value{?s} to report at."
    ), class = "rTFT_missing_covariates")
  for (col in setdiff(names(fit$data), names(newdata)))
    newdata[[col]] <- fit$data[[col]][1]
  newdata
}


# .fit_draws() -----------------------------------------------------------------

# Take the posterior draws of the model's parameters (loga, logb, Tc, the baseline
# and co-stressor parameters) out of a fitted model, at the covariate values in
# each row of `newdata`, with random effects at zero: the population-level answer.
# Returns a named list, one matrix per parameter, with a row per draw and a column
# per row of `newdata`.
#
# This is the one place the reporting functions reach into a fitted model; every
# later step works on plain draws. Draw i of one parameter and draw i of another
# are one joint posterior sample, and their pairing is what carries the correlation
# between parameters into every derived quantity. So every parameter is read at the
# same `draw_ids`.
.fit_draws <- function(fit, pars, newdata, draw_ids = NULL) {
  if (!inherits(fit, "brmsfit"))
    cli::cli_abort("{.arg fit} must be a fitted model ({.cls brmsfit}).",
                   class = "rTFT_not_brmsfit")
  unknown <- setdiff(pars, names(fit$formula$pforms))
  if (length(unknown))
    cli::cli_abort(c(
      "The fitted model has no parameter {.val {unknown}}.",
      "i" = "Its parameters are {.val {names(fit$formula$pforms)}}."
    ), class = "rTFT_missing_nlpar")
  newdata <- .complete_newdata(fit, newdata, pars)
  if (is.null(draw_ids)) draw_ids <- seq_len(brms::ndraws(fit))

  out <- purrr::map(pars, \(p)
    brms::posterior_linpred(fit, newdata = newdata, nlpar = p, re_formula = NA,
                            draw_ids = draw_ids))
  stats::setNames(out, pars)
}


# ===========================================================================
# 6. The fitted failure rate, from a fitted model
# ===========================================================================

# M_draws() --------------------------------------------------------------------

#' The fitted failure rate over temperature, for each posterior draw
#'
#' Evaluates the failure rate M(T) of a model fitted with [fit_hazard()] at each
#' temperature in `Ta`, for every posterior draw, at the covariate values in each row
#' of `newdata`. The rate is computed from the same expression the model was fitted
#' with, and the failure-rate function, reference temperature, softening steepness
#' and co-stressors are read from the fit, so none of them is entered by hand. Random
#' effects are set to zero, so the rate is that of a typical group.
#'
#' The draws are returned as they are, not summarised, for plotting the rate with
#' its uncertainty or for further calculation.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param Ta A numeric vector of temperatures at which to evaluate the rate.
#' @param newdata A data frame of covariate values, one row per setting at which to
#'   evaluate the rate. It must hold every covariate in the model's parameter
#'   formulas; other columns the model needs are filled in. May be `NULL` for a model
#'   without covariates.
#' @param levels A named numeric vector giving the level of each co-stressor, on its
#'   measurement scale, e.g. `c(DO = 6)` for 6 mg/L of dissolved oxygen. Required
#'   when the model has co-stressors.
#' @param n_draws An optional single number, the number of posterior draws to keep,
#'   evenly spaced through the posterior. `NULL` (the default) keeps them all.
#' @return A tibble with one row per draw, temperature and row of `newdata`: the
#'   columns of `newdata`, `.draw` (the draw's number in the fit, so draws from
#'   different calls can be matched), `Ta` and the failure rate `M`.
#' @examples
#' \dontrun{
#' fit  <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#' rate <- M_draws(fit, Ta = seq(25, 40, by = 0.5), newdata = data.frame(mass = c(8, 12)))
#' }
#' @export
M_draws <- function(fit, Ta, newdata = NULL, levels = NULL, n_draws = NULL) {
  st   <- .fit_settings(fit)
  pars <- c(M_catalogue[[st$m_form]]$dpars, .costressor_dpars(st$costressors))
  if (is.null(newdata)) newdata <- data.frame(row.names = 1L)

  ids <- .thin_ids(brms::ndraws(fit), n_draws)
  dr  <- .fit_draws(fit, pars, newdata, draw_ids = ids)

  # The rows of newdata are stacked as further draws, so the rate for all of them is
  # evaluated in one call: the draws of row 1, then those of row 2, and so on. The
  # fit records the steepness it used; for a function without a threshold that is
  # 0, which is passed on as no override at all.
  nd  <- length(ids)
  nr  <- nrow(newdata)
  k   <- if (M_catalogue[[st$m_form]]$k_sw != 0) st$k_sw
  M   <- .M_draws(st$m_form, purrr::map(dr, as.vector), Ta, st$Tr, k,
                  st$costressors, levels)

  # One row per draw, temperature and newdata row, ordered by newdata row, then
  # draw, then temperature.
  row_of  <- rep(rep(seq_len(nr), each = nd), times = length(Ta))
  draw_of <- rep(rep(ids, times = nr), times = length(Ta))
  Ta_of   <- rep(Ta, each = nd * nr)
  ord     <- order(row_of, draw_of, Ta_of)
  dplyr::bind_cols(
    tibble::as_tibble(newdata[row_of[ord], , drop = FALSE]),
    tibble::tibble(.draw = draw_of[ord], Ta = Ta_of[ord], M = as.vector(M)[ord]))
}


# ===========================================================================
# 7. Thermal tolerance from a fitted model
# ===========================================================================

# .baseline_draws() / .param_draws() / .with_newdata() / .report_unreached() --

# Each draw's baseline shape and scale, in the (shape, scale) form H0_of() and
# H0inv_of() take: 1 and 1 for the memoryless baseline, the rate constant gamma and
# an unused 1 for the Gompertz, exp(logshape) and the reference time `tau_ref` the
# scale is fixed at for the Weibull, and exp(logshape), exp(logscale) for the
# log-logistic.
.baseline_draws <- function(baseline, d, nd, tau_ref = NULL) {
  if (baseline == "weibull" && is.null(tau_ref))
    cli::cli_abort(c(
      "This model has a {.val weibull} baseline but no recorded reference time.",
      "i" = "It was fitted before the Weibull scale was fixed at {.field tau_ref}; refit it with {.fn fit_hazard}."
    ), class = "rTFT_no_tau_ref")
  switch(baseline,
    mem      = list(shape = rep(1, nd),       scale = rep(1, nd)),
    gompertz = list(shape = d$gam,            scale = rep(1, nd)),
    weibull  = list(shape = exp(d$logshape),  scale = rep(tau_ref, nd)),
    llogis   = list(shape = exp(d$logshape),  scale = exp(d$logscale)))
}

# The draws of each reported parameter, on its natural scale: a, b and c from their
# logarithms, Tc as it is, the baseline's shape, its scale where that is estimated
# (the log-logistic) or its rate constant, and each co-stressor's coefficients. Which parameters are reported is decided by the
# failure-rate function, the baseline and the co-stressors, never by which draws
# happen to be in `d`. Where a co-stressor interacts with temperature,
# the thermal sensitivity at the chosen levels, b + sum of dT z, is added as
# `b_at_levels`. z is added only where the thermal death time curve is straight in
# log time against temperature, the exponential rate with the memoryless or Weibull
# baseline, where the time to a given fraction failing scales as M^(-1/k), so
# log10 t falls by b / (k log 10) per degree and z = k log(10) / b.
.param_draws <- function(m_form, baseline, d, costressors = list(), z_at = list(),
                         tau_ref = NULL) {
  nd     <- length(d$loga)
  m_pars <- M_catalogue[[m_form]]$dpars
  out    <- list(a = exp(d$loga))
  if ("logb" %in% m_pars) out$b  <- exp(d$logb)
  if ("logc" %in% m_pars) out$c  <- exp(d$logc)
  if ("Tc"   %in% m_pars) out$Tc <- d$Tc
  bl <- .baseline_draws(baseline, d, nd, tau_ref)
  if (baseline %in% c("weibull", "llogis")) out$shape <- bl$shape
  if (baseline == "llogis")                 out$scale <- bl$scale
  if (baseline == "gompertz") out$gamma <- bl$shape
  for (cs in costressors) {
    out[[paste0("c", cs$name)]] <- d[[paste0("c", cs$name)]]
    if (cs$interaction) out[[paste0("dT", cs$name)]] <- d[[paste0("dT", cs$name)]]
  }
  if (!is.null(out$b)) {
    interacting <- purrr::keep(costressors, \(cs) cs$interaction)
    b_at <- out$b
    for (cs in interacting) b_at <- b_at + d[[paste0("dT", cs$name)]] * z_at[[cs$name]]
    if (length(interacting)) out$b_at_levels <- b_at
    if (m_form == "exp" && baseline %in% c("mem", "weibull"))
      out$z <- bl$shape * log(10) / b_at
  }
  out
}

# Put one row of newdata in front of every row of a summary table.
.with_newdata <- function(row, tbl) {
  dplyr::bind_cols(tibble::as_tibble(row[rep(1L, nrow(tbl)), , drop = FALSE]), tbl)
}

# One message listing the quantities with draws that never reach them, and the
# largest share of such draws across the rows of newdata, or nothing when there are
# none. A message rather than a warning: an unreached draw is part of the answer,
# not a fault.
.report_unreached <- function(tables) {
  hits <- purrr::imap(tables, \(tbl, nm) {
    hit <- tbl[tbl$n_unreached > 0, , drop = FALSE]
    if (!nrow(hit)) return(NULL)
    label <- switch(nm,
      sctmax      = sprintf("sCTmax (p = %g, t_ref = %g)", hit$p, hit$t_ref),
      tdt_curves  = sprintf("time until p = %g fail, at some temperatures", hit$p),
      ttl_surface = "survival",
      hit$quantity)
    worst <- tapply(hit$n_unreached / hit$n_draws, label, max)
    sprintf("%s: up to %.1f%% of draws", names(worst), 100 * worst)
  })
  lines <- unlist(hits, use.names = FALSE)
  if (!length(lines)) return(invisible(NULL))
  cli::cli_inform(c(
    "i" = "Some posterior draws never reach these quantities, and count as infinite:",
    stats::setNames(lines, rep("*", length(lines))),
    "i" = "See the section on such draws in {.fn derive_tdt_ttl}'s help for how to read them."
  ), class = "rTFT_unreached")
}


# .tdt_draws() -----------------------------------------------------------------

# The time until a fraction p has failed under a constant temperature, for each
# draw (a row of M) at each temperature (a column): M(T) H0(t_p) = -log(1 - p), so
# t_p = H0^-1(-log(1 - p) / M(T)). Two cases give a time that is never reached,
# returned as Inf: a cumulative hazard that levels off below the target (a Gompertz
# baseline with a negative rate constant, where H0^-1 returns NaN), and a rate of
# zero (a zero function at Tr), where the target is infinite.
.tdt_draws <- function(M, p, H0inv, shape, scale) {
  x  <- -log1p(-p) / M
  tp <- H0inv(x, shape, scale)
  tp[is.nan(tp) | is.infinite(x)] <- Inf
  tp
}


# .stop_no_t_ref() / .derive_draws() / .summarise_derived() --------------------

# Stop when no reference exposure time was given. `t_ref` has no default: its unit
# is that of the data the model was fitted to, which the fit does not record, and a
# default of 60 would silently mean 60 hours for a model fitted in hours.
.stop_no_t_ref <- function() {
  cli::cli_abort(c(
    "Give {.arg t_ref}, the reference exposure time for sCTmax.",
    "i" = "It is in the time unit of the data the model was fitted to, e.g. {.code t_ref = 60} for one hour when times are in minutes."
  ), class = "rTFT_bad_t_ref")
}

# Compute, for a fitted model, the posterior draws of everything derive_tdt_ttl()
# reports, before any summary: for each row of `newdata`, the draws of the
# parameters, of sCTmax for each p, of the time to failure for each p at each
# temperature, and of survival at each temperature and elapsed time. It checks its
# arguments and sets the default grids. derive_tdt_ttl() summarises these draws, and
# average_tdt_ttl() pools them across models first, so both rest on one calculation.
.derive_draws <- function(fit, newdata, p, t_ref, Ta, tau, levels, n_draws) {
  st <- .fit_settings(fit)
  if (!is.numeric(p) || !length(p) || anyNA(p) || any(p <= 0 | p >= 1))
    cli::cli_abort(c(
      "{.arg p} must hold fractions of individuals failing, each strictly between 0 and 1.",
      "x" = "You supplied {.val {p}}."
    ), class = "rTFT_bad_p")
  if (!is.numeric(t_ref) || length(t_ref) != 1L || !is.finite(t_ref) || t_ref <= 0)
    cli::cli_abort(c(
      "{.arg t_ref} must be a single positive exposure time, in the time unit of the data the model was fitted to.",
      "x" = "You supplied {.val {t_ref}}."
    ), class = "rTFT_bad_t_ref")
  if (is.null(newdata)) newdata <- data.frame(row.names = 1L)

  # The grids default to the temperatures and elapsed times the data cover: 50
  # temperatures evenly spaced, and 50 times evenly spaced on a log scale, since
  # times to failure and survival change over orders of magnitude in time.
  if (is.null(Ta))  Ta  <- seq(st$Ta_range[1], st$Ta_range[2], length.out = 50)
  if (is.null(tau)) tau <- exp(seq(log(st$tau_range[1]), log(st$tau_range[2]),
                                   length.out = 50))
  if (!is.numeric(Ta) || !length(Ta) || !all(is.finite(Ta)))
    cli::cli_abort("{.arg Ta} must be a numeric vector of finite temperatures.",
                   class = "rTFT_bad_Ta")
  if (!is.numeric(tau) || !length(tau) || !all(is.finite(tau)) || any(tau < 0))
    cli::cli_abort(c(
      "{.arg tau} must be a numeric vector of elapsed exposure times, none negative.",
      "i" = "They are in the time unit of the data the model was fitted to."
    ), class = "rTFT_bad_tau")

  cs   <- st$costressors
  pars <- c(M_catalogue[[st$m_form]]$dpars, h0_catalogue[[st$baseline]]$dpars,
            .costressor_dpars(cs))
  ids  <- .thin_ids(brms::ndraws(fit), n_draws)
  nd   <- length(ids)
  dr   <- .fit_draws(fit, pars, newdata, draw_ids = ids)
  k    <- if (M_catalogue[[st$m_form]]$k_sw != 0) st$k_sw

  # Check the co-stressor levels once, warning once about any the model does not
  # have; after that only the levels in use are passed on, so the calculations for
  # each row and each p do not repeat the warning.
  .rate_setup(st$m_form, purrr::map(dr, \(m) m[, 1]), st$Tr, k, cs, levels)
  cs_names <- vapply(cs, `[[`, character(1), "name")
  levels   <- if (length(cs)) levels[cs_names]
  z_at     <- stats::setNames(purrr::map(cs, \(x) .costressor_z(x, levels[[x$name]])),
                              cs_names)
  H0    <- H0_of(st$baseline)
  H0inv <- H0inv_of(st$baseline)

  rows <- purrr::map(seq_len(nrow(newdata)), \(r) {
    d  <- purrr::map(dr, \(m) m[, r])
    bl <- .baseline_draws(st$baseline, d, nd, st$tau_ref)
    sctmax <- purrr::map(p, \(pp) {
      target <- -log1p(-pp) / H0(t_ref, bl$shape, bl$scale)
      .sctmax_draws(st$m_form, d, target, st$Tr, k, cs, levels)
    })
    # The thermal death time curves: for each p, the time until that fraction has
    # failed at each temperature of the grid.
    M   <- .M_draws(st$m_form, d, Ta, st$Tr, k, cs, levels)
    tdt <- purrr::map(p, \(pp) .tdt_draws(M, pp, H0inv, bl$shape, bl$scale))
    # The thermal tolerance landscape: survival exp(-M(T) H0(tau)) over the grid of
    # temperatures and elapsed times. H0(tau) is laid out one row per draw, and
    # multiplying by the column of M at one temperature scales each draw's row by
    # its own rate.
    H0_tau  <- H0(matrix(tau, nd, length(tau), byrow = TRUE), bl$shape, bl$scale)
    surface <- purrr::map(seq_along(Ta), \(j) exp(-M[, j] * H0_tau))
    list(params = .param_draws(st$m_form, st$baseline, d, cs, z_at, st$tau_ref),
         sctmax = sctmax, tdt = tdt, surface = surface)
  })
  list(newdata = newdata, p = p, t_ref = t_ref, Ta = Ta, tau = tau, rows = rows)
}

# Summarise the draws from .derive_draws() into the tables derive_tdt_ttl() returns,
# each carrying the columns of newdata: `sctmax` with a row per p, `tdt_curves` with
# a row per p and temperature, `ttl_surface` with a row per temperature and elapsed
# time, and, unless `params` is FALSE, `params` with a row per parameter.
.summarise_derived <- function(draws, prob, params = TRUE) {
  p   <- draws$p
  Ta  <- draws$Ta
  tau <- draws$tau
  per_row <- purrr::map(seq_along(draws$rows), \(r) {
    d   <- draws$rows[[r]]
    row <- draws$newdata[r, , drop = FALSE]
    out <- list(
      sctmax = dplyr::bind_rows(purrr::map(seq_along(p), \(i)
        tibble::tibble(p = p[i], t_ref = draws$t_ref, .draw_summary(d$sctmax[[i]], prob)))),
      tdt_curves = dplyr::bind_rows(purrr::map(seq_along(p), \(i)
        tibble::tibble(p = p[i], Ta = Ta, .summarise_columns(d$tdt[[i]], prob)))),
      ttl_surface = dplyr::bind_rows(purrr::map(seq_along(Ta), \(j)
        tibble::tibble(Ta = Ta[j], tau = tau, .summarise_columns(d$surface[[j]], prob)))))
    if (params)
      out <- c(list(params = dplyr::bind_rows(purrr::imap(d$params, \(x, nm)
        tibble::tibble(quantity = nm, .draw_summary(x, prob))))), out)
    purrr::map(out, \(tbl) .with_newdata(row, tbl))
  })
  tables <- names(per_row[[1]])
  purrr::map(stats::setNames(tables, tables),
             \(tbl) dplyr::bind_rows(purrr::map(per_row, tbl)))
}


# derive_tdt_ttl() -------------------------------------------------------------

#' Thermal tolerance derived from a fitted failure-rate model
#'
#' Derives, from a model fitted with [fit_hazard()], the failure-rate parameters on
#' their natural scale, the static critical thermal maximum sCTmax, the thermal death
#' time curves and the thermal tolerance landscape, each with its posterior median
#' and credible interval, at the covariate values in each row of `newdata`. The model's 
#' settings (failure-rate and baseline functions, reference temperature, softening 
#' steepness, co-stressors) are read from the fit. Random effects are set to zero, so 
#' the results are those of a typical group.
#'
#' sCTmax is the constant temperature at which a fraction `p` of individuals fail within
#' the reference exposure time `t_ref`, the static critical thermal maximum of
#' Jorgensen et al. (2021), who report it for one hour. Under a constant
#' temperature, survival is exp(-M(T) H0(t_ref)), so for each posterior draw sCTmax
#' solves M(T) = -log(1 - p) / H0(t_ref). It is not the critical thermal maximum a
#' ramping experiment records (their dynamic CTmax).
#'
#' The parameters reported are a, and where the failure-rate function has them b,
#' c and Tc; the baseline's shape (Weibull, log-logistic), its scale (log-logistic;
#' the Weibull's is fixed at the reference time, see [h0_catalogue]) or its rate
#' constant gamma (Gompertz); and each co-stressor's coefficients. With a
#' co-stressor that interacts with temperature, `b_at_levels` is the thermal
#' sensitivity b + dT z at the chosen levels. The thermal sensitivity z, the
#' temperature change that shortens the time to failure tenfold, is reported only
#' where the thermal death time curve is a straight line in log time against
#' temperature, so that one number describes it: the exponential failure rate with
#' the memoryless or Weibull baseline, where z = k log(10) / b (k = 1 for the
#' memoryless baseline, the Weibull shape otherwise). Elsewhere the curve bends and
#' no single z stands for it.
#'
#' The thermal death time curve for a fraction `p` gives, at each constant
#' temperature in `Ta`, the time until that fraction has failed:
#' t_p = H0^-1(-log(1 - p) / M(T)). The thermal tolerance landscape gives survival
#' exp(-M(T) H0(tau)) over the grid of temperatures `Ta` and elapsed exposure times
#' `tau`. "Thermal death time" is the literature's name for the curve, whatever the
#' failure endpoint.
#'
#' @section Draws that never reach the quantity:
#' A posterior draw for which a quantity is never reached (an sCTmax the rate
#' cannot attain within 1000 degrees of the reference temperature, or at all; a time
#' to failure that never comes, under a Gompertz baseline whose cumulative hazard
#' levels off) counts as infinite rather than being dropped, so the summaries
#' describe the whole posterior. `n_unreached` counts these draws, and a message
#' lists the quantities that have any. With coverage `prob`, let u = (1 - prob) / 2
#' (0.025 for a 95% interval):
#'
#' * fewer than a share u of the draws unreached: both bounds are finite; report
#'   the summary as it stands.
#' * between u and half: the outer bound is infinite; report it as unbounded (for
#'   example "above 40 degC", "not reached"), and consider whether the question lies
#'   beyond the temperatures or exposure times the data cover.
#' * half or more: the median is infinite; the quantity is not estimable there. Ask
#'   at a lower `p` or a shorter `t_ref`.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param newdata A data frame of covariate values, one row per setting at which to
#'   derive. It must hold every covariate in the model's parameter formulas. May be
#'   `NULL` for a model without covariates.
#' @param p A numeric vector of the fractions of individuals failing at which sCTmax and
#'   the thermal death time curves are derived, each strictly between 0 and 1.
#' @param t_ref A single positive number, the reference exposure time for sCTmax, in
#'   the time unit of the data the model was fitted to (for example `60` for one hour
#'   when times are in minutes). It has no default, because the fit does not record
#'   that unit.
#' @param Ta An optional numeric vector of the constant temperatures for the thermal
#'   death time curves and the tolerance landscape. `NULL` (the default) takes 50
#'   evenly spaced across the temperatures the data cover.
#' @param tau An optional numeric vector of the elapsed exposure times for the
#'   tolerance landscape, in the time unit of the data, none negative. `NULL` (the
#'   default) takes 50 spaced evenly on a log scale across the times the data cover.
#'   A time of zero is allowed, and gives a survival of 1, so that a survival curve
#'   can start there.
#' @param levels A named numeric vector giving the level of each co-stressor, on its
#'   measurement scale, e.g. `c(DO = 6)`. Required when the model has co-stressors.
#' @param prob A single number strictly between 0 and 1, the coverage of the
#'   credible intervals.
#' @param n_draws An optional single number, the number of posterior draws to use,
#'   evenly spaced through the posterior. `NULL` uses them all.
#' @return A list of four tibbles, each starting with the columns of `newdata` and
#'   ending with the posterior summary `median`, `lower`, `upper`, `prob`, `n_draws`
#'   and `n_unreached`:
#'   * `params`: one row per parameter, named in `quantity`;
#'   * `sctmax`: one row per value of `p`, with `p` and `t_ref`;
#'   * `tdt_curves`: one row per value of `p` and temperature, with `p` and `Ta`; the
#'     summary is of the time until that fraction has failed;
#'   * `ttl_surface`: one row per temperature and elapsed time, with `Ta` and `tau`;
#'     the summary is of survival.
#' @references
#' Jorgensen, L. B., Malte, H., Orsted, M., Klahn, N. A. and Overgaard, J. (2021).
#' A unifying model to estimate thermal tolerance limits in ectotherms across
#' static, dynamic and fluctuating exposures to thermal stress. *Scientific
#' Reports* 11, 12840.
#' @examples
#' \dontrun{
#' fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#' tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), p = 0.5,
#'                       t_ref = 60)
#' tol$sctmax
#' }
#' @export
derive_tdt_ttl <- function(fit, newdata = NULL, p = c(0.1, 0.5, 0.9), t_ref,
                           Ta = NULL, tau = NULL, levels = NULL, prob = 0.95,
                           n_draws = 1000) {
  if (missing(t_ref)) .stop_no_t_ref()
  .check_prob(prob)
  draws <- .derive_draws(fit, newdata, p, t_ref, Ta, tau, levels, n_draws)
  out   <- .summarise_derived(draws, prob)
  .report_unreached(out)
  out
}
