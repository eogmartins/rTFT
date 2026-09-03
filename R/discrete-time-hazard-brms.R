# discrete-time-hazard-brms.R

# ===========================================================================
# 0a. Numerically stable cloglog Bernoulli family (interval mode)
# ===========================================================================

# The built-in bernoulli(link = "cloglog") forms p = inv_cloglog(eta), which
# underflows to 1 for eta > ~3.6 and returns -Inf log-lik for survived
# intervals. This family computes the log-lik straight from eta.

# cloglog_bern_family() --------------------------------------------------

#' Complementary log-log Bernoulli response, on the identity link
#'
#' The response distribution for interval-censored failures. brms's own
#' `bernoulli(link = "cloglog")` forms the failure probability first, and that
#' probability rounds to exactly one once the linear predictor passes about 3.6,
#' after which every survived interval returns minus infinity. This is not a
#' rare edge case in thermal failure data: near CTmax the failure rate grows
#' very large, so the cumulative hazard over an interval (and with it the linear
#' predictor) routinely passes that threshold, and the fit becomes unusable.
#' This family keeps the predictor on the identity link and forms the log
#' density from it directly, which is numerically stable across the whole range.
#'
#' @return A `customfamily` for use with the rTFT package.
#' @export
cloglog_bern_family <- function() {
  brms::custom_family(
    "cloglog_bern",
    dpars = "mu", links = "identity", type = "int", lb = 0L, ub = 1L
  )
}

## .cloglog_bern_stanvars: Stan cloglog_bern_lpmf -------------------------------

# This returns log P(y | mu): the log-probability of the outcome recorded in the
# row, one branch per outcome. The survived branch reads as a bare -lam because
# the survival probability is exp(-lam), so its log is just -lam (the log and the
# exp cancel). The failed branch is log(1 - exp(-lam)). When lam is small,
# exp(-lam) is just below 1, so computing 1 - exp(-lam) directly loses most of its
# precision and underflows to log(0) = -Inf once lam is tiny enough. log1m_exp
# avoids these issues by never subtracting exp(-lam) from 1.
.cloglog_bern_stanvars <- brms::stanvar(
  scode = "
  real cloglog_bern_lpmf(int y, real mu) {
    real lam = exp(mu);            // mu is eta (identity link)
    if (y == 1) return log1m_exp(-lam);
    else        return -lam;
  }
  ",
  block = "functions"
)

## log_lik_cloglog_bern() -----------------------------------------------------

#' Pointwise log-likelihood of a `cloglog_bern` observation
#'
#' The R counterpart of the `cloglog_bern_lpmf` in Stan. It returns the same log
#' density in R. brms needs this function after fitting to compute the pointwise
#' log-likelihood (used by LOO and WAIC), as the Stan function cannot be called
#' from R.
#'
#' @param i Observation index.
#' @param prep A `brmsprep` object.
#' @return For observation `i`, a vector of log-likelihood values, one per draw;
#'   brms binds these columns into the draws-by-observations matrix.
#' @keywords internal
#' @export
log_lik_cloglog_bern <- function(i, prep) {
  mu  <- brms::get_dpar(prep, "mu", i = i)
  y   <- prep$data$Y[i]
  lam <- exp(mu)
  y * log(-expm1(-lam)) + (1 - y) * (-lam)
}

## posterior_predict_cloglog_bern() -------------------------------------------

#' Draw from the posterior predictive of a `cloglog_bern` observation
#'
#' @param i Observation index.
#' @param prep A `brmsprep` object.
#' @param ... Unused.
#' @return For observation `i`, a vector of zeros and ones, one per draw; brms
#'   binds these columns into the draws-by-observations matrix.
#' @keywords internal
#' @export
posterior_predict_cloglog_bern <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p  <- pmin(1 - exp(-exp(mu)), 1 - 1e-15)
  stats::rbinom(length(p), 1L, p)
}

## posterior_epred_cloglog_bern() ---------------------------------------------

#' Expected value of a `cloglog_bern` observation
#'
#' The probability that the failure falls within the interval,
#' 1 - exp(-cumulative hazard).
#'
#' @param prep A `brmsprep` object.
#' @return A draws-by-observations matrix.
#' @keywords internal
#' @export
posterior_epred_cloglog_bern <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  pmin(1 - exp(-exp(mu)), 1 - 1e-15)
}


# ===========================================================================
# 0b. Unified family for mixed interval-censored + known-time failures
# ===========================================================================

# hazard_mixed_family() --------------------------------------------------

#' Response distribution for a table mixing interval-censored and known-time failures
#'
#' Some failures are known only to have fallen within an interval, others at a
#' known time. This family fits both kinds in one model. A row where the animal
#' survived its interval contributes the same amount either way, so only the row
#' recording the failure differs: an interval-censored failure contributes the
#' probability of failing somewhere inside the interval, whereas a failure seen at a
#' known time contributes the density at that recorded time. A per-row flag, passed
#' to brms through `vint()`, tells the model which of the two to use, so a single
#' fit can contain animals whose failures were recorded as either kind.
#'
#' @return A `customfamily` for use with the rTFT package.
#' @export
hazard_mixed_family <- function() {
  brms::custom_family(
    "hazard_mixed",
    dpars = "mu", links = "identity", type = "int", lb = 0L, ub = 1L,
    vars = "vint1[n]"
  )
}

## .hazard_mixed_stanvars: Stan hazard_mixed_lpmf -----------------------------

.hazard_mixed_stanvars <- brms::stanvar(
  scode = "
  real hazard_mixed_lpmf(int y, real mu, int vint1) {
    real lam = exp(mu);
    if (y == 0)     return -lam;              // survived interval (both kinds)
    if (vint1 == 0) return log1m_exp(-lam);   // interval-censored failure (cloglog)
    return mu - lam;                          // known-time failure (Poisson piecewise-exp)
  }
  ",
  block = "functions"
)

## log_lik_hazard_mixed() -----------------------------------------------------

#' Pointwise log-likelihood of a `hazard_mixed` observation
#'
#' The R counterpart of the `hazard_mixed_lpmf` in Stan, reading the per-row
#' known-time flag from `prep$data$vint1`. It returns the same log density in R.
#' brms needs this function after fitting to compute the pointwise log-likelihood
#' (used by LOO and WAIC), as the Stan function cannot be called from R.
#'
#' @param i Observation index.
#' @param prep A `brmsprep` object.
#' @return For observation `i`, a vector of log-likelihood values, one per draw;
#'   brms binds these columns into the draws-by-observations matrix.
#' @keywords internal
#' @export
log_lik_hazard_mixed <- function(i, prep) {
  mu  <- brms::get_dpar(prep, "mu", i = i)
  y   <- prep$data$Y[i]
  kt  <- prep$data$vint1[i]
  lam <- exp(mu)
  if (y == 0) return(-lam)
  if (kt == 0) return(log(-expm1(-lam)))      # interval-censored failure
  mu - lam                                    # known-time failure
}

## posterior_predict_hazard_mixed() -------------------------------------------

# This and posterior_epred model the 0/1 failure indicator through the interval
# failure probability 1 - exp(-lam), and ignore vint1. A 0/1 response can only be
# predicted as a Bernoulli, and 1 - exp(-lam) is the probability of failing in the
# row's window either way, so a pp_check on failures is unaffected. The known failure
# time enters the likelihood (log_lik) but not the indicator, so only LOO and the
# predictive differ in how they weight a known-time failure.

#' Draw from the posterior predictive of a `hazard_mixed` observation
#'
#' @param i Observation index.
#' @param prep A `brmsprep` object.
#' @param ... Unused.
#' @return For observation `i`, a vector of zeros and ones, one per draw; brms
#'   binds these columns into the draws-by-observations matrix.
#' @keywords internal
#' @export
posterior_predict_hazard_mixed <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p  <- pmin(1 - exp(-exp(mu)), 1 - 1e-15)
  stats::rbinom(length(p), 1L, p)
}

## posterior_epred_hazard_mixed() ---------------------------------------------

#' Expected value of a `hazard_mixed` observation
#'
#' The probability that the failure falls within the interval,
#' 1 - exp(-cumulative hazard).
#'
#' @param prep A `brmsprep` object.
#' @return A draws-by-observations matrix.
#' @keywords internal
#' @export
posterior_epred_hazard_mixed <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  pmin(1 - exp(-exp(mu)), 1 - 1e-15)
}
