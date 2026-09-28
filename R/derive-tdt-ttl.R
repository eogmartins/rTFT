# derive-tdt-ttl.R

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
# matrix whose rows are draws, and recycle column-wise over draws. The Gompertz shape
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
#' the hazard rise or fall with how long the animal has already been exposed. Shape
#' and scale arrive as per-draw vectors, and elapsed time may be a vector or a matrix
#' whose rows are draws.
#'
#' For the Gompertz shape the first argument is the rate constant gamma and the
#' second is ignored. A negative gamma is a valid value and means a hazard that
#' falls the longer the animal has been exposed. The cumulative hazard then levels
#' off: it rises towards 1/|gamma| but never reaches it, instead of growing without
#' limit. Survival therefore never drops to zero, and the model predicts that some
#' animals never fail (a cure fraction, in survival-analysis terms). [H0inv_of()]
#' gives the time by which a chosen fraction of animals has failed; when that
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
