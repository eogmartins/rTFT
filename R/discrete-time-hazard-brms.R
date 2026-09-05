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


# ===========================================================================
# 1. Catalogues and log M(T) / log h0(tau) expressions
# ===========================================================================

# M_catalogue -----------------------------------------------------------------

#' Failure rate functions M(T)
#'
#' This catalogue contains every function the failure rate may take: the nine
#' failure rate functions M(T) of Kingsolver and Umbanhowar (2018, Table 1), each
#' catalogue key named for its model there. Without a threshold: `const` (a
#' constant rate), `exp` (exponential), `zeroexp` (exponential minus its value at
#' the reference temperature `Tr`, so the rate is zero there), `power` (a power
#' of temperature above a baseline rate) and `zeropower` (the same power with no
#' baseline, zero at `Tr`). The threshold variants `const_threshold`,
#' `exp_threshold`, `zeroexp_threshold` and `zeropower_threshold` hold the rate at
#' zero below a fitted critical temperature Tc; above Tc it is constant for
#' `const_threshold` and accelerates for the other three. The exponential
#' functions give a straight thermal death time curve, whereas the power functions
#' give a curved one.
#'
#' `Tr` (`T0` in Kingsolver and Umbanhowar's notation) is a single reference
#' temperature set for the whole model, not the per-trial starting temperature. The
#' exponential is anchored at `Tr` and the zero forms vanish there, but because
#' the failure rate must depend only on the current temperature, `Tr` cannot be
#' each trial's own starting temperature. Trials that begin at different
#' temperatures (acclimation groups, field conditions) enter through their
#' temperature trajectories, not through `Tr`.
#'
#' To keep fitting, simulation and reporting from diverging, all three evaluate
#' the same expression strings, which are held in [brms_log_M].
#'
#' Each entry gives the integer the Stan code dispatches on (`stan_id`), the
#' parameters the function needs (`dpars`), and the steepness of the threshold
#' switch (`k_sw`). The rate parameters a, b and c are represented on the log
#' scale (`loga`, `logb`, `logc`), so the linear predictor is unconstrained and
#' their strictly positive values are recovered by exponentiation; the critical
#' temperature Tc is left on the natural scale, since it is not sign-constrained.
#' `k_sw` is roughly the inverse of the width of the transition in temperature, so
#' a larger value is a sharper switch over a narrower band. The switch is smoothed
#' rather than hard because Stan fits the model with Hamiltonian Monte Carlo,
#' which explores the posterior by following the gradient of the log-density with
#' respect to each parameter; the gradient of a hard step with respect to the
#' critical temperature is zero on either side of it and undefined at the step
#' itself, so the sampler gets no slope to move that temperature along, whereas a
#' smooth switch supplies one.
#'
#' @references
#' Kingsolver, J. G. and Umbanhowar, J. (2018). The analysis and interpretation
#' of critical temperatures. Journal of Experimental Biology 221, jeb167858.
#' \doi{10.1242/jeb.167858}
#'
#' @format A named list, one entry per failure rate function M(T).
#' @export
M_catalogue <- list(
  const               = list(stan_id = 1L, dpars = "loga",                    k_sw = 0),
  const_threshold     = list(stan_id = 2L, dpars = c("loga", "Tc"),           k_sw = 50),
  exp                 = list(stan_id = 3L, dpars = c("loga", "logb"),         k_sw = 0),
  exp_threshold       = list(stan_id = 4L, dpars = c("loga", "logb", "Tc"),   k_sw = 50),
  zeroexp             = list(stan_id = 5L, dpars = c("loga", "logb"),         k_sw = 0),
  zeroexp_threshold   = list(stan_id = 6L, dpars = c("loga", "logb", "Tc"),   k_sw = 50),
  power               = list(stan_id = 7L, dpars = c("loga", "logb", "logc"), k_sw = 0),
  zeropower           = list(stan_id = 8L, dpars = c("loga", "logb"),         k_sw = 0),
  zeropower_threshold = list(stan_id = 9L, dpars = c("loga", "logb", "Tc"),   k_sw = 50)
)


# h0_catalogue -----------------------------------------------------------------

#' Baseline hazard functions h0(tau)
#'
#' This catalogue contains the four baseline hazard functions h0(tau), the factor
#' in elapsed exposure time that multiplies the failure rate M(T). `mem` is the
#' memoryless case, in which the hazard does not depend on how long the animal
#' has already been exposed; `weibull`, `llogis` and `gompertz` allow
#' it to rise or fall with elapsed exposure.
#'
#' `mem` is the baseline canonical to Kingsolver and Umbanhowar (2018), and is
#' the default throughout. Their framework has no baseline-hazard-in-time
#' component at all: the only hazard is M(T) evaluated along the temperature
#' trajectory T(t). The other three are extensions of their framework.
#'
#' `gompertz`, h0(tau) = exp(gamma * tau), is the only one of the four baselines
#' whose hazard increases exponentially with elapsed exposure time, which is also
#' the shape an exponential M(T) exhibits as temperature increases through a ramp.
#' Because the two share that shape, a single ramp rate cannot separate them, as
#' it only sees the combined acceleration, not how much comes from the increasing
#' temperature and how much from time under exposure. Two or more ramp rates can
#' separate them, because the temperature part steepens with the ramp rate while the
#' elapsed time part does not. The same confound occurs between `weibull`, whose
#' hazard is a power of elapsed time, and a power M(T) (`zeropower` or
#' `zeropower_threshold`), since a power of temperature along a ramp is also a power of
#' elapsed time. So the two are confounded under a single ramp and, as for
#' `gompertz`, are separated by two or more ramp rates. `llogis`, whose hazard
#' rises then falls, matches none of the failure-rate functions of temperature,
#' which all rise monotonically, so it is not confounded with any of them.
#'
#' A Gompertz baseline hazard is a different object from the Gompertz
#' distribution of failure times that Kingsolver and Umbanhowar's Eqn 6 derives.
#' The distribution is an outcome of the model: an exponential M(T) traversed
#' along a linear ramp, under a memoryless baseline, produces Gompertz-distributed
#' failure times. The baseline hazard presented here is an assumption about how
#' risk depends on elapsed exposure time, independent of temperature.
#'
#' Each entry gives the integer the Stan code dispatches on (`stan_id`) and the
#' parameters of each baseline function (`dpars`), as in [M_catalogue]. `mem` has
#' no parameters. `weibull` and `llogis` have shape and scale parameters on the
#' log scale (`logshape`, `logscale`). `gompertz` has the rate constant gamma (`gam`)
#' parameter in reciprocal time, left on the natural scale and unconstrained in
#' sign (a positive value is a hazard rising with elapsed exposure time and a negative
#' one a hazard falling). The amplitude of the Gompertz hazard is not included,
#' because it is the same quantity as `loga` in M(T), and `mem` omits a parameter
#' for the same reason.
#'
#' @references
#' Kingsolver, J. G. and Umbanhowar, J. (2018). The analysis and interpretation
#' of critical temperatures. Journal of Experimental Biology 221, jeb167858.
#' \doi{10.1242/jeb.167858}
#'
#' @format A named list, one entry per baseline function.
#' @export
h0_catalogue <- list(
  mem      = list(stan_id = 1L, dpars = character(0)),
  weibull  = list(stan_id = 2L, dpars = c("logshape", "logscale")),
  llogis   = list(stan_id = 3L, dpars = c("logshape", "logscale")),
  gompertz = list(stan_id = 4L, dpars = "gam")
)


# brms_log_M -------------------------------------------------------------------

#' log M(T) for each failure rate function, as a brms expression string
#'
#' Each entry is the expression string of log M(T) for one function of [M_catalogue],
#' which brms uses to build the model formula. The strings are written in terms
#' of the parameter names and the data columns `Ta` (the row's temperature: the
#' interval midpoint, or the exact temperature for a known-time failure) and `Tr`
#' (the starting temperature). The threshold functions add a small constant
#' (`1e-12`) inside a logarithm so it stays finite at the critical temperature.
#' The temperature column is named `Ta` rather than
#' Kingsolver and Umbanhowar's `T`, because `T` is an alias for `TRUE` in R and a
#' column of that name would silently shadow it in the data-handling code. The
#' threshold functions also reference `k_sw`, the softening steepness, taken from
#' the `k_sw` field of [M_catalogue] when the model formula is built by brms. The
#' strings are the single shared definition that fitting, simulation and reporting
#' all evaluate.
#'
#' @format A named list, one expression string per failure rate function of
#'   [M_catalogue].
#' @export
brms_log_M <- list(
  const               = "loga",
  const_threshold     = "loga - log1p(exp(-k_sw * (Ta - Tc)))",
  exp                 = "loga + exp(logb) * (Ta - Tr)",
  exp_threshold       = "loga + exp(logb) * (Ta - Tc) - log1p(exp(-k_sw * (Ta - Tc)))",
  zeroexp             = "loga + log(exp(exp(logb) * (Ta - Tr)) - 1)",
  zeroexp_threshold   = "loga + log(abs(exp(exp(logb) * (Ta - Tc)) - 1) + 1e-12) - log1p(exp(-k_sw * (Ta - Tc)))",
  power               = "log(exp(logc) + exp(loga) * (Ta - Tr)^exp(logb))",
  zeropower           = "loga + exp(logb) * log(Ta - Tr)",
  zeropower_threshold = "loga + exp(logb) * log(abs(Ta - Tc) + 1e-12) - log1p(exp(-k_sw * (Ta - Tc)))"
)


# brms_log_h0 ------------------------------------------------------------------

#' log h0(tau) for each baseline function, as a brms expression string
#'
#' Each entry is the expression string of log h0(tau) for one function of
#' [h0_catalogue], which brms uses to build the model formula. The strings are
#' written in terms of the parameter names and the data columns `log_tau` (the log
#' of the row's elapsed exposure time) and, for `gompertz`, `tau` itself (the
#' elapsed exposure time: the midpoint for an interval row, or the exact time for
#' a known-time failure). `mem` is the empty string, because a memoryless
#' baseline contributes nothing: h0(tau) = 1 and the hazard is M(T) alone.
#'
#' Each string is the baseline **hazard** on the log scale, not the density. For
#' `llogis` the two differ by one `log1p` term: the hazard has a single
#' `- log1p(...)`, the density would have two.
#'
#' @format A named list, one expression string per baseline function of
#'   [h0_catalogue].
#' @export
brms_log_h0 <- list(
  mem      = "",
  weibull  = "logshape - logscale + (exp(logshape) - 1) * (log_tau - logscale)",
  llogis   = paste0("logshape - logscale + (exp(logshape) - 1) * (log_tau - logscale)",
                    " - log1p(exp(exp(logshape) * (log_tau - logscale)))"),
  gompertz = "gam * tau"
)
