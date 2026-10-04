# discrete-time-hazard-brms.R

# Functions to fit thermal failure times (i.e., the time at which an individual fails 
# under a thermal exposure) as discrete-time complementary log-log survival models 
# in brms, following the failure-rate framework of Kingsolver and Umbanhowar (2018). Each
# at-risk interval of a trial contributes a Bernoulli failure indicator whose hazard
# is h0(tau) * M(T): a failure rate M(T) set by the current temperature, and a
# baseline hazard h0(tau) in elapsed exposure. The failure rate takes one of a
# catalogue of forms (M_catalogue) and the baseline one of a catalogue of shapes
# (h0_catalogue). Both are written once as log-scale expression strings (brms_log_M,
# brms_log_h0) that fitting, simulation and reporting all read, so the model has a
# single definition.
#
# The path is build_hazard_data() -> build_hazard_spec() -> build_hazard_prior() ->
# fit_hazard(): reshape one row per individual into one row per at-risk interval,
# assemble the nonlinear formula and its response family, supply weakly-informative
# default priors, and hand it to brms. Failure times may be interval-censored,
# observed at a known time, or a mix of the two; each has its own response
# family (the numerically stable cloglog Bernoulli family here, brms's Poisson, or
# the unified mixed family). A non-temperature co-stressor Z (e.g., dissolved oxygen, pH)
# enters through costressor(), turning the failure rate into a surface
# M(T, Z). ramp_path() and assess_ramp_linearity() help when the exposure is a
# controlled linear temperature ramp.

# ===========================================================================
# 0a. Numerically stable cloglog Bernoulli family (interval mode)
# ===========================================================================

# The built-in bernoulli(link = "cloglog") forms p = inv_cloglog(eta), which
# underflows to 1 for eta > ~3.6 and returns -Inf log-lik for survived
# intervals. This family computes the log-lik straight from eta. (Continuous
# mode does not need this -- brms's built-in poisson() is stable.)

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
#' @return A `customfamily` for use with [fit_hazard()].
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
#' known time. This family fits both kinds in one model. A row where the individual
#' survived its interval contributes the same amount either way, so only the row
#' recording the failure differs: an interval-censored failure contributes the
#' probability of failing somewhere inside the interval, whereas a failure seen at a
#' known time contributes the density at that recorded time. A per-row flag, passed
#' to brms through `vint()`, tells the model which of the two to use, so a single
#' fit can contain individuals whose failures were recorded as either kind.
#'
#' @return A `customfamily` for use with [fit_hazard()].
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
#' memoryless case, in which the hazard does not depend on how long the individual
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


# ===========================================================================
# Internal helpers: argument validation
# ===========================================================================

# Shared by the functions below so the checks (and their cli condition classes) are
# defined once. Each validates a user argument against the catalogues and, on a bad
# value, raises the same cli error the rest of the package uses.

# Resolve the observation mode, keeping match.arg()'s partial matching and
# defaulting but rethrowing its plain stop() as a cli condition.
.match_observation <- function(observation) {
  choices <- c("interval", "continuous", "mixed")
  tryCatch(
    match.arg(observation, choices),
    error = function(e) {
      show <- cli::cli_vec(choices, list("vec-last" = ", or "))
      cli::cli_abort(c(
        "{.arg observation} must be one of {.val {show}}.",
        "x" = "You supplied {.val {observation}}."
      ), class = "rTFT_bad_observation")
    }
  )
}

# The M_catalogue entry for a failure-rate function, or a cli error naming the keys.
.m_form_spec <- function(m_form) {
  spec <- M_catalogue[[m_form]]
  if (is.null(spec))
    cli::cli_abort(c(
      "Unknown {.arg m_form}: {.val {m_form}}.",
      "i" = "It must be a key of {.code M_catalogue}: {.val {names(M_catalogue)}}."
    ), class = "rTFT_unknown_m_form")
  spec
}

# The h0_catalogue entry for a baseline function, or a cli error naming the keys.
.baseline_spec <- function(baseline) {
  spec <- h0_catalogue[[baseline]]
  if (is.null(spec))
    cli::cli_abort(c(
      "Unknown {.arg baseline}: {.val {baseline}}.",
      "i" = "It must be a key of {.code h0_catalogue}: {.val {names(h0_catalogue)}}."
    ), class = "rTFT_unknown_baseline")
  spec
}

# Choose the softening steepness a failure-rate function is fitted and reported
# with: the user's `k_sw` when one is given, otherwise the catalogue value. Only the
# threshold functions have a gate to soften, so a `k_sw` given for any other function
# is ignored with a warning. Fitting and reporting both call this, so a model fitted
# with a custom steepness is also reported with it.
.resolve_k_sw <- function(m_form, k_sw = NULL) {
  m_spec <- .m_form_spec(m_form)
  is_threshold <- m_spec$k_sw != 0
  if (!is.null(k_sw)) {
    if (!is.numeric(k_sw) || length(k_sw) != 1L || !is.finite(k_sw) || k_sw <= 0)
      cli::cli_abort(c(
        "{.arg k_sw} must be a single positive number.",
        "x" = "You supplied {.val {k_sw}}."
      ), class = "rTFT_bad_k_sw")
    if (!is_threshold)
      cli::cli_warn(c(
        "!" = "{.arg k_sw} is ignored for {.val {m_form}}, which has no threshold to soften.",
        "i" = "It applies only to the threshold functions: {.val {names(M_catalogue)[vapply(M_catalogue, function(x) x$k_sw != 0, logical(1))]}}."
      ), class = "rTFT_k_sw_ignored")
  }
  if (!is.null(k_sw) && is_threshold) k_sw else m_spec$k_sw
}

# Build log M(T) as one expression string: the catalogue string for the failure-rate
# function, with the softening steepness written in as a number (so it stays a
# value, not a symbol the sampler would try to estimate), and each co-stressor's
# main effect and temperature interaction appended, which turns the rate into a
# surface M(T, Z). The model formula and every rate reported from a fit are built
# from this one string.
.log_M_string <- function(m_form, k_sw_used, costressors = list()) {
  log_M <- brms_log_M[[m_form]]
  if (.m_form_spec(m_form)$k_sw != 0)
    log_M <- gsub("k_sw", format(k_sw_used, scientific = FALSE), log_M, fixed = TRUE)
  cs_add <- .costressor_logM(costressors)
  if (nzchar(cs_add)) log_M <- paste(log_M, "+", cs_add)
  log_M
}


# ===========================================================================
# Co-stressor declarations
# ===========================================================================

# costressor() -----------------------------------------------------------------

#' Declare a co-stressor for a failure-rate hazard model
#'
#' Describes one non-temperature environmental variable (Z) the individual is exposed to
#' and that co-determines the failure rate: dissolved oxygen, salinity, pH, a
#' toxicant concentration, or, for a terrestrial ectotherm, humidity. The failure
#' rate then becomes a surface M(T, Z). The co-stressor enters log M(T) as its own
#' term plus, by default, an interaction with temperature, so a crossed
#' temperature-by-stressor design (a ramp run at two or more stressor levels) can
#' estimate how the stressor changes the thermal response. Build one for each
#' co-stressor and pass them together, as a list, to [build_hazard_data],
#' [build_hazard_spec], [build_hazard_prior] and [fit_hazard], so all four read the
#' same declaration.
#'
#' The raw measurement `X` becomes a standardised column
#' z = direction * (X - ref) / scale, which is zero at the benign reference. There
#' the failure rate reduces to the pure thermal M(T): the amplitude `loga` and
#' critical temperature `Tc` keep their unstressed meaning, and a model with the
#' co-stressor held at its reference is the same as one without it. This follows
#' from centring the column at a chosen reference rather than at the sample mean.
#' 
#' The three pieces the standardisation needs (a reference, a scale, and a
#' direction) are all required, as none has a safe default. A domain-chosen scale
#' (per mg/L of oxygen, per unit concentration) is preferred to a data standard
#' deviation, which is unstable with only two or three treatment levels. As an 
#' example, for dissolved oxygen with a reference at air-saturation (ref = 8.5 mg/L), 
#' a scale of 1 mg/L (scale = 1), and less oxygen meaning more stress (direction = -1),
#' z = -1 * (X - 8.5) / 1 = 8.5 - X, the oxygen deficit is: z = 0 at
#' saturation and z = 2.5 at 6 mg/L (below saturation). The direction re-orients every
#' stressor onto a common axis on which a larger z always means more stress, so the
#' main effect `c<name>` has the same sign meaning for every stressor: a positive value
#' is a hazard that rises under stress.
#'
#' This signed-linear column assumes the stressor harm runs one way. For a stressor that is
#' harmful in either direction (salinity or pH, where a departure above or below the
#' optimum both add stress), supply a `transform` function so z measures distance from the
#' reference rather than a signed deficit: a squared standardised deviation,
#' transform = function(u) u^2, makes z = ((X - ref) / scale)^2, zero at the optimum
#' and rising on both sides. The `direction` then no longer applies, since z is
#' symmetric.
#'
#' The co-stressor contributes two parameters, named from `name`. The main effect
#' `c<name>` is the change in the log failure rate per `scale` of stressor beyond
#' the reference, and the interaction `dT<name>` (kept unless switched off) is how
#' that change varies with temperature. These are the names brms estimates and that
#' [build_hazard_prior] sets priors on. `c<name>` is a continuous dose term: with
#' three or more stressor levels it is a per-`scale` slope through them, but with
#' only two it cannot be told apart from a per-treatment shift of the amplitude
#' `loga`. A main effect on its own is in fact just that shift: it is equivalent to
#' a covariate on the amplitude, `loga ~ 1 + z` in [build_hazard_spec]'s
#' `dpar_formulas`. The additive `dT<name>` term is not: no `dpar_formulas` entry
#' reproduces it, so it is the part `costressor()` supplies. A temperature
#' interaction of a different shape can still be had through `dpar_formulas` (the
#' multiplicative form in \sQuote{Two kinds of temperature interaction} below); what
#' is unique to `costressor()` is this particular additive term.
#'
#' Because `c<name>` and `dT<name>` are ordinary parameters, `dpar_formulas` in
#' [build_hazard_spec] can put covariates on either, so the stressor's effect can
#' itself vary: `cDO = ~ 1 + acclim` lets the oxygen sensitivity depend on
#' acclimation (a stressor-by-covariate interaction), and the same on `dT<name>`
#' lets the temperature interaction vary. These are higher-order terms, so they need
#' a design that varies temperature, the stressor and the covariate together enough
#' to identify them; the formula is built either way, but a thin design shows up as
#' poor mixing rather than a clean fit.
#'
#' @section Categorical stressor treatments:
#' `costressor()` expects a stressor measured on a numeric scale, so `ref`, `scale`
#' and `direction` mean something. Quantitative treatment levels belong here even
#' though they are "treatments": dissolved oxygen dosed at 8.5, 5 and 2 mg/L is
#' numeric, so pass the measured level as `column`. A genuinely unordered
#' categorical stressor, with no scale on which a reference and a per-unit step make
#' sense, does not. Model it through `dpar_formulas` in [build_hazard_spec] instead:
#' a factor on the amplitude `loga` gives a per-treatment level shift, the same as a
#' stressor main effect. To let the treatment effect vary with temperature, put the
#' factor on the parameter that governs the temperature response (`logb` for the
#' exponential forms, `Tc` for the threshold forms). That is a genuine
#' temperature-by-treatment interaction, though a different parameterisation from the
#' co-stressor's `dT<name>` term, as it scales the temperature slope rather than adding
#' to it (see \sQuote{Two kinds of temperature interaction} below). A two-level 
#' treatment can instead be encoded 0/1 and passed to `costressor()`, which then gives
#' the co-stressor's own main and interaction terms.
#'
#' @section Two kinds of temperature interaction:
#' The default `dT<name>` interaction is additive on the temperature slope: it adds
#' `dT<name> * (Ta - Tr) * z` to log M, so the per-degree effect shifts by a constant
#' amount per unit of stressor. This is the standard temperature-by-stressor
#' interaction, it works the same on every M(T) function, and it lets the slope
#' change sign, so a stressor can in principle flatten or reverse the thermal
#' effect. The sign of `dT<name>` then reads as the synergy or antagonism of
#' temperature and the stressor (on the proportional-hazards scale): positive
#' steepens the thermal slope, negative flattens it, and a large enough negative
#' value reverses it (an antagonistic crossover where warming becomes protective
#' under severe harm caused by the co-stressor). The alternative is to let the 
#' stressor rescale thermal sensitivity by a proportion, which keeps it positive: 
#' set `interaction = FALSE` (so the co-stressor adds only its main effect and the
#' `z_<name>` column) and put that column on the sensitivity parameter through 
#' `dpar_formulas` instead: `logb ~ 1 + z_<name>` for the exponential forms (for the 
#' threshold forms a factor on `Tc` shifts the threshold location, a different effect).
#' Choose the additive form for a conventional interaction and the freedom for the
#' slope to change sign (e.g., an antagonistic crossover), the multiplicative form 
#' when the biology is a proportional change in sensitivity that should stay positive
#' (it can attenuate the thermal slope but not reverse it).
#'
#' @param name A character string, a short name for the co-stressor, used to name
#'   its column (`z_<name>`) and its parameters (`c<name>`, `dT<name>`). Letters and
#'   digits only, starting with a letter, because brms allows no dot or underscore
#'   in a parameter name. For example `"DO"` gives the column `z_DO` and the
#'   parameters `cDO` and `dTDO`.
#' @param column A character string naming the column that holds the raw measurement
#'   `X`: a column of `ind` (one value per individual) for a stressor held fixed
#'   within a trial, or a column of `trajectory` for one that changes over the trial.
#' @param ref The benign, no-additional-stress reference level of the raw
#'   measurement (air-saturation for dissolved oxygen, zero for a toxicant, the
#'   control level for salinity or pH). The failure rate reduces to the pure M(T) here.
#' @param scale A single positive number, the unit a coefficient is expressed per
#'   (per mg/L or per kPa for oxygen, per unit concentration for a toxicant).
#' @param direction Either `1` or `-1`, set so that a positive `z` always means more
#'   stress: `-1` for dissolved oxygen (less oxygen is more stress), `1` for a
#'   toxicant (more is more stress).
#' @param interaction A single logical value. When `TRUE` (the default) the
#'   co-stressor also enters through a temperature interaction `dT<name>`, so it can
#'   change the temperature slope of log M; when `FALSE` it enters through its main
#'   effect alone. For a multiplicative alternative that rescales thermal sensitivity
#'   instead, see \sQuote{Two kinds of temperature interaction} below.
#' @param trajectory An optional data frame giving the co-stressor over time, with
#'   columns `path_id`, `time` and the `column` named above, in the shape of the
#'   `path` argument of [build_hazard_data]. Supply it for a stressor that changes
#'   during a trial; it is interpolated to the same interval midpoints as temperature.
#'   Leave `NULL` for a stressor held fixed within a trial, which is then read from `ind`.
#' @param transform An optional function of the standardised deviation
#'   (X - ref) / scale that returns the column `z`. Use it for a stressor whose harm
#'   is not signed-linear, for example a squared deviation for a U-shaped response to
#'   salinity or pH, where a departure either way is harmful. When `NULL` (the
#'   default) the signed-linear deficit direction * (X - ref) / scale is used.
#' @return A `costressor` object: a named list holding the declaration, for the
#'   `costressors` argument of [build_hazard_data], [build_hazard_spec],
#'   [build_hazard_prior] and [fit_hazard].
#' @examples
#' # Dissolved oxygen held at a fixed level per trial: benign at air-saturation
#' # (about 8.5 mg/L), expressed per mg/L, with less oxygen meaning more stress.
#' costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
#'
#' # A toxicant that changes over the trial, supplied as its own trajectory, and
#' # entering through its main effect alone (no temperature interaction).
#' tox_path <- data.frame(path_id = "A", time = c(0, 60), conc = c(0, 4))
#' costressor("Tox", column = "conc", ref = 0, scale = 1, direction = 1,
#'            interaction = FALSE, trajectory = tox_path)
#'
#' # Salinity, harmful above or below an optimum: a squared standardised deviation
#' # makes z measure distance from the optimum (about 15 PSU, per 5 PSU) in either
#' # direction. With a transform supplied, `direction` no longer applies.
#' costressor("Sal", column = "sal", ref = 15, scale = 5, direction = 1,
#'            transform = function(u) u^2)
#' @export
costressor <- function(name, column, ref, scale, direction,
                       interaction = TRUE, trajectory = NULL, transform = NULL) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !grepl("^[A-Za-z][A-Za-z0-9]*$", name))
    cli::cli_abort(c(
      "{.arg name} must be a single name of letters and digits, starting with a letter.",
      "x" = "You supplied {.val {name}}.",
      "i" = "It names the parameters {.code c<name>} and {.code dT<name>}, and brms allows no dot or underscore there."
    ), class = "rTFT_bad_costressor_name")

  if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column))
    cli::cli_abort(c(
      "{.arg column} must be a single column name.",
      "x" = "You supplied {.val {column}}."
    ), class = "rTFT_bad_costressor_column")

  # A single finite number, optionally strictly positive.
  chk_num1 <- function(x, arg, cls, positive = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && (!positive || x > 0)
    if (!ok)
      cli::cli_abort(c(
        "{.arg {arg}} must be a single {if (positive) 'positive ' else ''}number.",
        "x" = "You supplied {.val {x}}."
      ), class = cls)
  }
  chk_num1(ref,   "ref",   "rTFT_bad_costressor_ref")
  chk_num1(scale, "scale", "rTFT_bad_costressor_scale", positive = TRUE)

  if (!is.numeric(direction) || length(direction) != 1L || !direction %in% c(-1, 1))
    cli::cli_abort(c(
      "{.arg direction} must be either {.val {1}} or {.val {-1}}.",
      "x" = "You supplied {.val {direction}}.",
      "i" = "Set it so a positive {.code z} always means more stress."
    ), class = "rTFT_bad_costressor_direction")

  if (!is.logical(interaction) || length(interaction) != 1L || is.na(interaction))
    cli::cli_abort("{.arg interaction} must be a single {.code TRUE} or {.code FALSE}.",
                   class = "rTFT_bad_costressor_interaction")

  if (!is.null(trajectory) && !is.data.frame(trajectory))
    cli::cli_abort("{.arg trajectory} must be a data frame or {.code NULL}.",
                   class = "rTFT_bad_costressor_trajectory")

  if (!is.null(transform) && !is.function(transform))
    cli::cli_abort("{.arg transform} must be a function or {.code NULL}.",
                   class = "rTFT_bad_costressor_transform")

  structure(
    list(name = name, column = column, ref = ref, scale = scale,
         direction = direction, interaction = interaction,
         trajectory = trajectory, transform = transform),
    class = "rTFT_costressor"
  )
}


# .check_costressors() / .costressor_dpars() / .costressor_logM() --------------

# Validate the `costressors` list argument shared by the fitting path: a list of
# costressor() objects with distinct names. Returns invisibly so it can guard a
# function's entry.
.check_costressors <- function(costressors) {
  if (!length(costressors)) return(invisible(costressors))
  if (!is.list(costressors) ||
      !all(vapply(costressors, inherits, logical(1), "rTFT_costressor")))
    cli::cli_abort(c(
      "{.arg costressors} must be a list of {.fn costressor} objects.",
      "i" = "Build each one with {.fn costressor}."
    ), class = "rTFT_bad_costressors")
  nms <- vapply(costressors, `[[`, character(1), "name")
  dup <- unique(nms[duplicated(nms)])
  if (length(dup))
    cli::cli_abort(c(
      "Each co-stressor needs a distinct {.arg name}.",
      "x" = "Repeated: {.val {dup}}."
    ), class = "rTFT_duplicate_costressor")
  invisible(costressors)
}

# The nonlinear parameters a set of co-stressors adds: the main effect `c<name>`
# for each, plus the interaction `dT<name>` for those keeping the temperature
# interaction. build_hazard_spec adds these to the parameters needing a sub-formula
# and build_hazard_prior to those needing an intercept prior, so both read one source.
.costressor_dpars <- function(costressors) {
  if (!length(costressors)) return(character(0))
  unlist(purrr::map(costressors, \(cs)
    c(paste0("c", cs$name), if (cs$interaction) paste0("dT", cs$name))),
    use.names = FALSE)
}

# The terms a set of co-stressors appends to log M(T), as one expression string
# (empty when there are none): the main effect `c<name> * z_<name>` and, when the
# interaction is kept, `dT<name> * (Ta - Tr) * z_<name>`.
.costressor_logM <- function(costressors) {
  if (!length(costressors)) return("")
  terms <- unlist(purrr::map(costressors, \(cs) {
    z    <- paste0("z_", cs$name)
    main <- paste0("c", cs$name, " * ", z)
    if (cs$interaction) c(main, paste0("dT", cs$name, " * (Ta - Tr) * ", z)) else main
  }), use.names = FALSE)
  paste(terms, collapse = " + ")
}

# Standardise raw co-stressor measurements X to z: the signed-linear deficit
# direction * (X - ref) / scale, or the supplied transform of (X - ref) / scale.
# build_hazard_data() uses it to build the z_<name> column the model is fitted to,
# and the reporting functions to place a chosen level on the same scale.
.costressor_z <- function(cs, X) {
  u <- (X - cs$ref) / cs$scale
  if (is.null(cs$transform)) cs$direction * u else cs$transform(u)
}


# ===========================================================================
# 2. Build the row-per-interval long table
# ===========================================================================

# build_hazard_data() ----------------------------------------------------------

#' Build the row-per-interval table a failure-rate model is fitted to
#'
#' Reshapes one row per individual, together with a temperature trajectory, into
#' one row per individual and at-risk interval. An interval-censored failure
#' gives the failure its whole censoring interval; a known-time one truncates
#' the final at-risk interval at the recorded time. The cumulative hazard is
#' integrated numerically over the at-risk intervals, with the failure-check times
#' serving as the integration nodes. For a known-time failure, the fit integrates
#' on the checks (its final interval ending exactly at the failure time) and no
#' finer sub-grid is required but can be provided.
#'
#' @param ind A data frame with one row per individual. Required: `id`, `path_id` (key
#'   into `path`), `event` (1 = failed, 0 = right-censored), plus the covariate
#'   columns referenced in the model. Timing depends on how the failure was
#'   observed: interval-censored individuals have `t_upper` (the check at which
#'   failure was found, or a survivor's last observed check); known-time ones
#'   have `t_event` (the exact failure or censoring time). For
#'   `observation = "mixed"` the data also needs an `obs_type` column per
#'   individual (`"interval"` or `"continuous"`), and each type supplies its own
#'   timing column.
#' @param path A data frame giving the temperature trajectory for each path.
#'   Required columns: `path_id`, `time`, `temp`. The function interpolates these
#'   points linearly in time, so a minimal two-point path (start and end) encodes a
#'   linear ramp. Interpolation uses the full series, so a temperature data logger
#'   whose records are finer than the checks still supplies the true temperature at
#'   each interval midpoint (never a check-to-check straight line); the checks control
#'   only the integration grid, not the temperature. Temperatures outside the recorded
#'   time range are held flat at the nearest endpoint rather than returned as `NA`,
#'   because the at-risk grid (the checks, the exact failure time, and the `t = 0`
#'   start) need not sit strictly inside the logged times: a check or failure time
#'   can fall just past the last reading, or the first interval can begin before the
#'   first reading, and an `NA` temperature there would propagate into `Ta` and the
#'   linear predictor.
#' @param observation How each failure time was observed. One of `"interval"`
#'   (every failure known only to fall within a censoring interval), `"continuous"`
#'   (every failure time known exactly), or `"mixed"` (both kinds in one table,
#'   each individual's kind given by its `obs_type` column). This determines which
#'   response family the table is fitted with, and whether the `exact` flag is
#'   added; the response column is `failed` in every mode.
#' @param covariates A character vector of column names in `ind` to carry through
#'   onto every row of the output table, so they are available in the model formula
#'   as fixed-effect predictors or as grouping factors for random effects. A trial
#'   identifier is the common case: include it here, then place a trial-level random
#'   intercept on a parameter in [build_hazard_spec], e.g.
#'   `dpar_formulas = list(loga = ~ 1 + (1 | trial))`. The individual identifier `id`
#'   is always kept, so an individual-level grouping term needs no entry here.
#' @param Tr The reference temperature, written to a constant `Tr` column so the
#'   functions can form `Ta - Tr`. It is a single model reference, not the
#'   per-trial starting temperature (see [M_catalogue]); the per-trial starts live
#'   in the temperature trajectories (`path`). For the forms that take a logarithm or
#'   fractional power of `Ta - Tr` (`zeroexp`, `power`, `zeropower`) it must lie at
#'   or below the minimum temperature in the data, so it defaults (`NULL`) to the
#'   coldest temperature in `path`. It must be set to a common value across datasets 
#'   when comparing `loga`, whose meaning is anchored at `Tr`.
#' @param checks An optional data frame with `path_id`, `time` giving the
#'   failure-check schedule per path (i.e., the times at which each individual was
#'   inspected for failure). Defaults to the distinct `path` times (logger reading =
#'   check). As noted above, these times double as the integration nodes for the
#'   cumulative hazard, so their spacing sets the midpoint-rule resolution. For a
#'   known-time failure, a regular grid can be used instead (see `grid_dt`).
#' @param grid_dt An optional numeric giving a fixed time step (in the units of
#'   `path$time`). If set, a known-time failure is integrated on a regular grid of
#'   this step, from `t = 0` to the failure time, instead of on `checks`, providing 
#'   a finer cumulative hazard at the cost of more rows. Defaults to `NULL` (integrate
#'   on the checks).
#' @param costressors An optional list of [costressor] objects declaring
#'   non-temperature stressors (dissolved oxygen and the like). Each adds a
#'   standardised column `z_<name>`, read from its own trajectory (interpolated to
#'   the interval midpoints) or, for a stressor held fixed within a trial, from a
#'   column of `ind`. Defaults to none.
#' @return A long-format data frame, one row per individual and at-risk interval,
#'   holding the identifiers `id` and `path_id`; the model columns `Ta`, `Tr`,
#'   `tau`, `log_tau`, `exposure`, `log_exposure`; the carried covariates; one
#'   `z_<name>` column per co-stressor; the 0/1 response column `failed`; and, for
#'   mixed observation, the integer `exact` flag the family branches on.
#' @examples
#' # A single linear ramp encoded as two points: 20 -> 35 degrees over 60 minutes,
#' # with failure checked for every 10 minutes.
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))
#'
#' # Interval-censored: one failure found at the 30-minute check, one survivor.
#' ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
#' build_hazard_data(ind, path, "interval", checks = checks)
#'
#' # Known-time: the exact failure time, integrated on the same checks.
#' ind_c <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
#' build_hazard_data(ind_c, path, "continuous", checks = checks)
#'
#' # A finer, regular integration grid instead of the checks (known-time only);
#' # the last interval still closes exactly on the failure time.
#' build_hazard_data(ind_c, path, "continuous", grid_dt = 5)
#'
#' # Mixed: each individual carries its kind in obs_type; the unused timing
#' # column is left NA.
#' ind_m <- data.frame(id = 1:2, path_id = "A", event = 1,
#'                     obs_type = c("interval", "continuous"),
#'                     t_upper = c(30, NA), t_event = c(NA, 42))
#' build_hazard_data(ind_m, path, "mixed", checks = checks)
#'
#' # A covariate carried through to every row, and a shared reference Tr.
#' ind$mass <- c(12.1, 9.8)
#' build_hazard_data(ind, path, "interval", covariates = "mass", Tr = 15, checks = checks)
#' @export
build_hazard_data <- function(ind, path, observation = c("interval", "continuous", "mixed"),
                              covariates = character(0), Tr = NULL, checks = NULL,
                              grid_dt = NULL, costressors = list()) {
  observation <- .match_observation(observation)
  .check_costressors(costressors)

  # Fail on a missing input column, rather than letting a NULL slip into an
  # approxfun() call: approxfun(time, NULL) does not error, it silently
  # reinterprets its arguments and returns garbage temperatures.
  need_cols <- function(df, cols, arg) {
    miss <- setdiff(cols, names(df))
    if (length(miss))
      cli::cli_abort("{.arg {arg}} is missing column{?s} {.field {miss}}.",
                     class = "rTFT_missing_column")
  }
  need_cols(ind, c("id", "path_id", "event"), "ind")
  need_cols(path, c("path_id", "time", "temp"), "path")
  if (!is.null(checks)) need_cols(checks, c("path_id", "time"), "checks")

  # In the two non-mixed modes, the observation mode fixes a single timing column
  # for the whole table: `t_upper` for interval, `t_event` for continuous. Guard it
  # the same way as the input columns above: a missing column reads back as NULL, and
  # findInterval() on it then returns a single at-risk interval rather than erroring. Mixed
  # mode carries a column per level and is guarded separately below.
  if (observation != "mixed") {
    timing <- c(interval = "t_upper", continuous = "t_event")[[observation]]
    need_cols(ind, timing, "ind")
    if (!all(is.finite(ind[[timing]])))
      cli::cli_abort(
        "{.arg ind} column {.field {timing}} is not finite for {sum(!is.finite(ind[[timing]]))} of {nrow(ind)} individual{?s}.",
        class = "rTFT_nonfinite_timing")
  }

  # Reference temperature defaults to the coldest point in the trajectories, which
  # keeps `Ta - Tr` positive for the log/power forms. A form aware check that a
  # user-supplied Tr is not above the data (needed only for zeroexp/power/zeropower)
  # belongs in the formula assembly step, which knows the function used.
  if (is.null(Tr)) Tr <- min(path$temp)

  # In mixed mode, each individual has its own scheme in `obs_type`. Guard the
  # levels explicitly, because build_one() reads that string twice and the two
  # reads disagree on anything but the exact levels "interval" and "continuous".
  # It picks the integration grid with `type == "interval"` (so any other value
  # takes the continuous grid, built from t_event) but sets the likelihood's
  # `exact` flag with `type == "continuous"` (so any other value comes out 0, i.e.
  # interval-censored). A typo or NA level would thus get a continuous grid while
  # being flagged interval-censored (two observation schemes for one individual, with
  # nothing to signal it). Rejecting unknown levels here keeps them aligned.
  if (observation == "mixed") {
    if (!"obs_type" %in% names(ind))
      cli::cli_abort(c(
        "{.arg observation} = {.val mixed} needs an {.field obs_type} column on {.arg ind}.",
        "i" = "It holds {.val interval} or {.val continuous} for each individual."
      ), class = "rTFT_missing_obs_type")
    ot  <- as.character(ind[["obs_type"]])
    ok  <- c("interval", "continuous")
    bad <- setdiff(unique(ot), ok)
    if (length(bad))
      cli::cli_abort(c(
        "Unrecognised {.field obs_type}: {.val {bad}}.",
        "i" = "The accepted levels are {.val interval} and {.val continuous}."
      ), class = "rTFT_bad_obs_type")
    needs <- c(interval = "t_upper", continuous = "t_event")
    for (lv in intersect(ok, unique(ot))) {
      col <- needs[[lv]]
      if (!col %in% names(ind))
        cli::cli_abort(
          "{.field obs_type} == {.val {lv}} reads column {.field {col}}, which {.arg ind} lacks.",
          class = "rTFT_missing_timing_col")
      have <- ind[[col]][ot == lv]
      if (!all(is.finite(have)))
        cli::cli_abort(
          "{.field obs_type} == {.val {lv}} reads {.field {col}}, but it is not finite for {sum(!is.finite(have))} of {length(have)} such individual{?s}.",
          class = "rTFT_nonfinite_timing")
    }
  }

  interp_of <- split(path, path$path_id) |>
    purrr::map(\(d) stats::approxfun(d$time, d$temp, rule = 2))
  chk_src   <- if (is.null(checks)) path else checks
  checks_of <- split(chk_src, chk_src$path_id) |>
    purrr::map(\(d) sort(unique(d$time)))

  # One individual's rows. `type` is "interval" or "continuous" (per individual in
  # mixed mode via obs_type; otherwise the global observation mode).
  build_one <- function(r, type) {
    pid <- as.character(r$path_id)
    f   <- interp_of[[pid]]
    ck  <- checks_of[[pid]]
    if (type == "interval") {
      kD <- max(findInterval(r$t_upper - 1e-9, ck), 1)   # failure in [ck[kD], ck[kD+1]]
      ts <- ck[seq_len(kD)]; te <- ck[seq_len(kD) + 1]   # failure row = full interval
    } else {
      td    <- r$t_event
      nodes <- if (is.null(grid_dt)) sort(unique(c(0, ck[ck > 0 & ck < td], td)))
               else { g <- seq(0, td, by = grid_dt); sort(unique(c(g, td))) }
      ts <- nodes[-length(nodes)]; te <- nodes[-1]       # failure row ends exactly at td
    }
    n   <- length(ts)
    # `tau` and `Ta` are taken at the midpoint of every at-risk interval, the
    # failure interval included: they are the node of the midpoint rule that
    # integrates the cumulative hazard over the interval, not the event time. A
    # known-time failure is exact in the interval boundary (te = t_event, hence
    # `exposure`), not in where the hazard is sampled.
    mid <- (ts + te) / 2
    tibble::tibble(id = r$id, path_id = r$path_id, tau = mid, exposure = te - ts,
                   Ta = f(mid), failed = as.integer(seq_len(n) == n & r$event == 1),
                   exact = as.integer(type == "continuous"))
  }

  type_of <- function(r) if (observation == "mixed") as.character(r$obs_type) else observation
  long <- purrr::map(seq_len(nrow(ind)),
                     \(i) { r <- ind[i, , drop = FALSE]; build_one(r, type_of(r)) }) |>
    purrr::list_rbind()

  cov_tbl <- dplyr::select(ind, id, dplyr::all_of(covariates))
  out <- long |>
    dplyr::left_join(cov_tbl, by = "id") |>
    dplyr::mutate(Tr = Tr, log_tau = log(tau), log_exposure = log(exposure))

  # One standardised column per co-stressor, z = direction * (X - ref) / scale (or
  # the supplied transform). The raw measurement X is read from the co-stressor's
  # own trajectory, interpolated to the interval midpoints exactly as temperature is
  # (the midpoint is already carried as `tau`), or, when no trajectory is given,
  # from a column of `ind` held fixed within a trial and broadcast to its rows.
  for (cs in costressors) {
    if (!is.null(cs$trajectory)) {
      tr <- cs$trajectory
      need_cols(tr, c("path_id", "time", cs$column),
                paste0("costressor(\"", cs$name, "\") trajectory"))
      f_of <- split(tr, tr$path_id) |>
        purrr::map(\(d) stats::approxfun(d$time, d[[cs$column]], rule = 2))
      X <- rep(NA_real_, nrow(out))
      for (pid in unique(out$path_id)) {
        fp <- f_of[[as.character(pid)]]
        if (is.null(fp))
          cli::cli_abort(
            "Co-stressor {.val {cs$name}} has no trajectory for path {.val {pid}}.",
            class = "rTFT_missing_costressor_path")
        at <- out$path_id == pid
        X[at] <- fp(out$tau[at])
      }
    } else {
      need_cols(ind, cs$column, paste0("ind (co-stressor ", cs$name, ")"))
      X <- ind[[cs$column]][match(out$id, ind$id)]
    }
    out[[paste0("z_", cs$name)]] <- .costressor_z(cs, X)
  }

  # The response is the same 0/1 failure indicator in every mode, named `failed`
  # throughout; only the family (and, for mixed, the `exact` flag) differs.
  if (observation != "mixed") out$exact <- NULL
  dplyr::relocate(out, id, dplyr::all_of(covariates))
}


# ===========================================================================
# 3. Formula and family builder
# ===========================================================================

# .rhs_str() -------------------------------------------------------------------

# The right-hand side of a one-sided formula, as a string:
#   ~ 1 + species + (1 | trial)  ->  "1 + species + (1 | trial)"
# Used to attach a user's covariate formula to a parameter's sub-formula.
.rhs_str <- function(f) paste(deparse(f[[length(f)]]), collapse = " ")


# build_hazard_spec() --------------------------------------------------------

#' Assemble the brms formula and family for a failure-rate hazard model
#'
#' Builds a nonlinear brms formula whose linear predictor is the log failure rate
#' log M(T), plus (unless the baseline is memoryless) the log baseline hazard
#' log h0(tau), plus the log exposure of each at-risk interval. The failure rate
#' and baseline are named by their keys in [M_catalogue] and [h0_catalogue], and
#' their log-scale expression strings ([brms_log_M], [brms_log_h0]) are pasted
#' straight into the predictor, so the formula fitted here is the same expression
#' [build_hazard_data] assumes and that reporting later evaluates.
#'
#' Each parameter the chosen functions need (`loga`, `logb`, `Tc`, and so on) gets
#' its own sub-formula, an intercept only (`~ 1`) by default. Pass a one-sided
#' formula in `dpar_formulas` to put covariates on a parameter; a formula named for
#' a parameter the chosen functions do not use is ignored with a warning.
#'
#' @param m_form A character string denoting a key of [M_catalogue] (e.g. `"exp"` or
#'   `"exp_threshold"`).
#' @param baseline A character string denoting a key of [h0_catalogue]: `"mem"`
#'   (default), `"weibull"`, `"llogis"` or `"gompertz"`.
#' @param dpar_formulas A named list of one-sided formulas, one per parameter you
#'   want covariates on. A parameter the functions need but you omit defaults to
#'   `~ 1` (intercept only); a name not among the parameters is ignored with a warning.
#'   For example, `list(loga = ~ 1 + acclim + (1 | trial), logb = ~ 1 + (1 | trial))`.
#' @param observation A character string denoting how each failure time was observed,
#'   as in [build_hazard_data]: `"interval"` (the [cloglog_bern_family]), `"continuous"`
#'   (a Poisson piecewise-exponential), or `"mixed"` (the [hazard_mixed_family], which
#'   reads the per-row `exact` flag through `vint()`).
#' @param k_sw An optional positive numeric representing the softening steepness of a
#'   threshold function's gate, in inverse temperature units (the transition band is
#'   roughly its reciprocal). This is an advanced override: leave `NULL` (the default)
#'   to use the value in [M_catalogue] (50), which suits most fits. A larger value is
#'   a sharper switch over a narrower band, but too large returns the sampler to the
#'   hard step geometry, where the critical temperature `Tc` mixes poorly (low
#'   effective sample size, sometimes divergences). This is a sampling-geometry
#'   effect, separate from identifiability: `Tc` is poorly identified whenever the
#'   observed temperatures do not bracket it, whatever `k_sw` is. Ignored, with a
#'   warning, for a function that has no threshold.
#' @param costressors An optional list of [costressor] objects. Each appends its
#'   terms to log M(T) (a main effect `c<name>`, and unless switched off a
#'   temperature interaction `dT<name>`), so the failure rate becomes a surface
#'   M(T, Z). Their coefficients join the parameters getting a sub-formula, so a
#'   `dpar_formulas` entry can put covariates on one. Defaults to none.
#' @return A named list with the brms `formula`, the `family`, and the `observation`
#'   mode.
#' @examples
#' # An exponential failure rate and memoryless baseline, with interval-censored
#' # observation (the default).
#' build_hazard_spec("exp")
#'
#' # An exponential-threshold failure rate with a Weibull baseline hazard in elapsed
#' # exposure time.
#' build_hazard_spec("exp_threshold", baseline = "weibull")
#'
#' # A different structure per parameter: acclim as a fixed covariate on the
#' # amplitude loga, and a trial-level random intercept on both loga and the
#' # temperature sensitivity logb (`acclim` and `trial` are carried through by
#' # build_hazard_data()).
#' build_hazard_spec("exp",
#'   dpar_formulas = list(loga = ~ 1 + acclim + (1 | trial), logb = ~ 1 + (1 | trial)))
#'
#' # Mixed interval-censored and known-time failures: the response holds the per-row
#' # `exact` flag through vint().
#' build_hazard_spec("exp", observation = "mixed")
#'
#' # A crossed temperature-by-oxygen design: dissolved oxygen enters log M as a main
#' # effect (cDO) and a temperature interaction (dTDO).
#' do <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
#' build_hazard_spec("exp", costressors = list(do))
#'
#' # Advanced: a gentler threshold switch than the catalogue default.
#' build_hazard_spec("exp_threshold", k_sw = 20)
#' @export
build_hazard_spec <- function(m_form, baseline = "mem", dpar_formulas = list(),
                                observation = c("interval", "continuous", "mixed"),
                                k_sw = NULL, costressors = list()) {
  observation <- .match_observation(observation)
  m_spec  <- .m_form_spec(m_form)
  h0_spec <- .baseline_spec(baseline)
  .check_costressors(costressors)

  needed <- c(m_spec$dpars, h0_spec$dpars, .costressor_dpars(costressors))

  # Paste the stored log-scale strings into the linear predictor: log M(T) (with
  # the softening steepness written in and any co-stressor terms appended), the log
  # baseline hazard (unless memoryless, whose string is empty), and the log exposure
  # of the interval. build_hazard_data() supplies Ta, Tr, tau, log_tau and
  # log_exposure; the parameters (loga, logb, Tc, ...) are what the model estimates.
  k_sw_used <- .resolve_k_sw(m_form, k_sw)   # validates k_sw, even where it is unused
  log_M     <- .log_M_string(m_form, k_sw_used, costressors)
  log_h0    <- brms_log_h0[[baseline]]

  eta  <- if (nzchar(log_h0)) paste(log_h0, "+", log_M, "+ log_exposure")
          else                paste(log_M, "+ log_exposure")
  lhs  <- if (observation == "mixed") "failed | vint(exact)" else "failed"
  main <- stats::as.formula(paste(lhs, "~", eta))

  extra <- setdiff(names(dpar_formulas), needed)
  if (length(extra))
    cli::cli_warn(c(
      "!" = "Ignoring {.arg dpar_formulas} not used by {.val {m_form}} + {.val {baseline}}: {.val {extra}}.",
      "i" = "The parameters in play are {.val {needed}}."
    ), class = "rTFT_unused_dpar_formulas")

  # One sub-formula per parameter: the user's covariate formula if given, else an
  # intercept. do.call(brms::bf, ...) assembles the nonlinear formula (nl = TRUE),
  # in which the main predictor is an expression in these parameters and the data.
  dp_fmls <- purrr::map(needed, \(dp) {
    rhs   <- dpar_formulas[[dp]]
    rhs_s <- if (is.null(rhs)) "1" else .rhs_str(rhs)
    stats::as.formula(paste(dp, "~", rhs_s))
  })

  bf_obj <- do.call(brms::bf, c(list(main, nl = TRUE), dp_fmls))
  fam <- switch(observation,
                interval   = cloglog_bern_family(),
                continuous = brms::brmsfamily("poisson", link = "log"),
                mixed      = hazard_mixed_family())

  list(formula = bf_obj, family = fam, observation = observation)
}


# ===========================================================================
# 4. Prior builder
# ===========================================================================

# .rhs_terms() / .is_bar() / .has_re() / .has_slopes() -------------------------

# Whether a parameter's one-sided formula carries fixed slopes and/or a random
# effect decides which extra priors build_hazard_prior() adds. Both read the
# formula's parse tree rather than its text, so a `(a | b)` group with inner
# parentheses (e.g., a brms grouping helper such as `(1 | gr(site, cor = TRUE))`)
# is still recognised as a single random-effect term.

# The additive terms of a one-sided formula's RHS, split on top-level `+` only, so a
# `(a | b)` group stays whole.
.rhs_terms <- function(f) {
  split_plus <- function(e)
    if (is.call(e) && identical(e[[1]], quote(`+`)))
      c(split_plus(e[[2]]), split_plus(e[[3]]))
    else list(e)
  split_plus(f[[length(f)]])
}

# Checks if a term is a random effect group `(... | ...)`, or `(... || ...)` for
# uncorrelated effects: a group parses as a parenthesis call wrapping a `|` or `||`
# call.
.is_bar <- function(term)
  is.call(term) && identical(term[[1]], quote(`(`)) && is.call(term[[2]]) &&
    (identical(term[[2]][[1]], quote(`|`)) || identical(term[[2]][[1]], quote(`||`)))

# Checks if a one-sided formula contain a random effect term
.has_re <- function(f)
  !is.null(f) && any(vapply(.rhs_terms(f), .is_bar, logical(1)))

# Checks if the formula has a fixed predictor beyond the intercept. Random effect terms
# are dropped, and the intercept indicators (the numeric 0 / 1) do not count, so only a
# named covariate or an interaction reads as a slope.
.has_slopes <- function(f) {
  if (is.null(f)) return(FALSE)
  terms <- .rhs_terms(f)
  fixed <- terms[!vapply(terms, .is_bar, logical(1))]
  any(!vapply(fixed, is.numeric, logical(1)))
}


# build_hazard_prior() ---------------------------------------------------------

#' Weakly-informative default priors for a failure-rate hazard model
#'
#' Supplies a sensible starting prior for the model [build_hazard_spec] assembles,
#' keyed to each parameter's intercept and robust to whatever covariates the
#' `dpar_formulas` add. For every parameter the chosen functions need, it sets a
#' weakly-informative prior on the intercept (the parameter's value at the reference
#' temperature and start of exposure), a generic Normal(0, 1) on any covariate
#' slopes, and an Exponential(2) on the standard deviation of any random effect. The
#' slope and random-effect priors are added only when the parameter actually has
#' those terms, so brms raises no "no parameter" warning. Override the centre of any
#' parameter through `intercept_priors`, or the shared slope and random-effect
#' standard-deviation priors through `slope_prior` and `sd_prior`; or replace the
#' whole set with your own `prior` when fitting with [fit_hazard()].
#'
#' The default intercept priors, all weakly-informative and meant to be changed to
#' suit your study, are:
#'
#' * `loga ~ normal(-5, 2)` -- log amplitude (the failure rate at the reference `Tr`)
#' * `logb ~ normal(-0.7, 1)` -- log temperature sensitivity
#' * `logc ~ normal(-5, 2)` -- log baseline term of the power form
#' * `Tc ~ normal(30, 3)` -- critical temperature
#' * `logshape ~ normal(0, 1)` -- log Weibull / log-logistic shape
#' * `logscale ~ normal(log(60), 1)` -- log Weibull / log-logistic scale
#' * `gam ~ normal(0, 1)` -- Gompertz rate
#'
#' with `normal(0, 1)` on any covariate slope and `exponential(2)` on any
#' random-effect standard deviation. Each co-stressor in `costressors` adds a
#' `normal(0, 1)` on its main effect `c<name>` and a tighter `normal(0, 0.5)` on its
#' temperature interaction `dT<name>`, both overridable through `intercept_priors`.
#'
#' Prior values for `loga` and `logscale` assume exposure time is measured in minutes
#' (a small failure rate per minute, and a trial about an hour long), so rescale them
#' for other time units. `Tc` assumes a temperature near many ectotherms' critical range.
#' `gam` is centred at zero, the memoryless case where the baseline hazard is flat in 
#' elapsed exposure time, so unless the data show a trend the fit stays near that flat 
#' baseline rather than assuming the hazard rises or falls over time.
#'
#' @param m_form A character string denoting a key of [M_catalogue] (e.g. `"exp"` or
#'   `"exp_threshold"`).
#' @param baseline A character string denoting a key of [h0_catalogue]: `"mem"`
#'   (default), `"weibull"`, `"llogis"` or `"gompertz"`.
#' @param dpar_formulas A named list of one-sided formulas, one per parameter, as in
#'   [build_hazard_spec]. A parameter's fixed and random terms decide which slope and
#'   standard-deviation priors are added.
#' @param intercept_priors A named list of prior strings, one per parameter, that
#'   replace the default intercept (centre) for those parameters, e.g.
#'   `list(Tc = "normal(26, 4)")`. A parameter left out keeps its default; a name
#'   that is not a parameter of the chosen functions is ignored with a warning.
#' @param slope_prior A single prior string for every covariate slope (a
#'   population-level effect beyond the intercept), applied to each parameter that
#'   has one. Defaults to `"normal(0, 1)"`.
#' @param sd_prior A single prior string for every random-effect standard deviation,
#'   applied to each parameter that has one. Defaults to `"exponential(2)"`.
#' @param costressors An optional list of [costressor] objects, as in
#'   [build_hazard_spec]. Each adds an intercept prior for its main effect and, when
#'   the temperature interaction is kept, its interaction coefficient. Defaults to none.
#' @return A `brmsprior`.
#' @examples
#' # Default weakly-informative priors for an exponential-threshold + Weibull model.
#' build_hazard_prior("exp_threshold", baseline = "weibull")
#'
#' # Override a couple of priors (a cooler species, a tighter amplitude) and narrow
#' # the covariate-slope prior; every other parameter keeps its default.
#' build_hazard_prior("exp_threshold", baseline = "weibull",
#'   dpar_formulas = list(loga = ~ 1 + acclim),
#'   intercept_priors = list(Tc = "normal(26, 4)", loga = "normal(-6, 1.5)"),
#'   slope_prior = "normal(0, 0.5)")
#' @export
build_hazard_prior <- function(m_form, baseline = "mem", dpar_formulas = list(),
                               intercept_priors = list(),
                               slope_prior = "normal(0, 1)",
                               sd_prior = "exponential(2)",
                               costressors = list()) {
  m_spec  <- .m_form_spec(m_form)
  h0_spec <- .baseline_spec(baseline)
  .check_costressors(costressors)
  needed  <- c(m_spec$dpars, h0_spec$dpars, .costressor_dpars(costressors))

  # Intercept prior per parameter, on the parameter's own scale (the log-scale
  # parameters carry priors in log units).
  icpt <- c(loga = "normal(-5, 2)", logb = "normal(-0.7, 1)", Tc = "normal(30, 3)",
            logc = "normal(-5, 2)", logshape = "normal(0, 1)",
            logscale = paste0("normal(", round(log(60), 3), ", 1)"),
            gam = "normal(0, 1)")

  # Each co-stressor's coefficients get a weakly-informative prior: the main effect
  # `c<name>` a wider Normal, the temperature interaction `dT<name>` a tighter one,
  # since the interaction is the part harder to identify. Both are overridable
  # through `intercept_priors`.
  for (cs in costressors) {
    icpt[[paste0("c", cs$name)]] <- "normal(0, 1)"
    if (cs$interaction) icpt[[paste0("dT", cs$name)]] <- "normal(0, 0.5)"
  }

  # A user's intercept overrides replace entries in that lookup, so anything left
  # unnamed keeps its default. Warn (rather than silently ignore) an override for a
  # parameter the chosen functions do not use.
  if (length(intercept_priors)) {
    if (is.null(names(intercept_priors)) || any(names(intercept_priors) == ""))
      cli::cli_abort("{.arg intercept_priors} must be a named list, one prior per parameter.",
                     class = "rTFT_bad_intercept_priors")
    unused <- setdiff(names(intercept_priors), needed)
    if (length(unused))
      cli::cli_warn(c(
        "!" = "Ignoring {.arg intercept_priors} not used by {.val {m_form}} + {.val {baseline}}: {.val {unused}}.",
        "i" = "The parameters in play are {.val {needed}}."
      ), class = "rTFT_unused_intercept_prior")
    icpt[names(intercept_priors)] <- unlist(intercept_priors)
  }

  # Accumulate in a list and combine once with do.call(c, ...), which uses the
  # brmsprior method for c() and preserves the class. The slope and
  # standard-deviation priors are the same string across parameters (each
  # overridable through slope_prior / sd_prior).
  plist <- list()
  for (dp in needed) {
    plist <- c(plist, list(brms::set_prior(icpt[[dp]], class = "b",
                                           coef = "Intercept", nlpar = dp)))
    if (.has_slopes(dpar_formulas[[dp]]))
      plist <- c(plist, list(brms::set_prior(slope_prior, class = "b", nlpar = dp)))
    if (.has_re(dpar_formulas[[dp]]))
      plist <- c(plist, list(brms::set_prior(sd_prior, class = "sd", nlpar = dp)))
  }
  do.call(c, plist)
}


# ===========================================================================
# 5. Model fitting
# ===========================================================================

# .check_Tr() ------------------------------------------------------------------

# Fail before fitting if the temperature reaches the reference Tr for a failure
# rate that takes a logarithm or fractional power of (Ta - Tr). `zeroexp` and
# `zeropower` are -Inf exactly at Ta == Tr (log of zero), so they need min(Ta)
# strictly above Tr; `power` is fine at Ta == Tr (0^b = 0) but NaN below it, so it
# needs min(Ta) at or above Tr. The other forms have no such constraint.
.check_Tr <- function(long, m_form) {
  needs_strict <- c("zeroexp", "zeropower")
  needs_ge     <- c("power")
  if (!m_form %in% c(needs_strict, needs_ge)) return(invisible(NULL))
  if (!all(c("Ta", "Tr") %in% names(long)))   return(invisible(NULL))  # the formula will complain
  Tr <- long$Tr[1]
  min_Ta <- min(long$Ta)
  strict <- m_form %in% needs_strict
  too_low <- if (strict) min_Ta <= Tr else min_Ta < Tr
  if (too_low) {
    rel <- if (strict) "above" else "at or above"
    cli::cli_abort(c(
      "The {.val {m_form}} failure rate requires the temperature to stay {rel} the reference {.arg Tr}.",
      "x" = "min(Ta) = {.val {min_Ta}}, Tr = {.val {Tr}}.",
      "i" = "Set a lower {.arg Tr} in {.fn build_hazard_data}, or use a form without a logarithm or fractional power of (Ta - Tr)."
    ), class = "rTFT_Tr_below_data")
  }
}

# .hazard_brm_args() -----------------------------------------------------------

# Assemble the argument list for brms::brm(), everything fit_hazard() does except
# the fit itself, so the wiring (formula, family, prior, stanvars, the Tr guard) can
# be checked without compiling a model.
.hazard_brm_args <- function(long, m_form, baseline = "mem", dpar_formulas = list(),
                             observation = "interval", prior = NULL,
                             intercept_priors = list(), slope_prior = NULL,
                             sd_prior = NULL, costressors = list(), k_sw = NULL) {
  observation <- .match_observation(observation)
  .check_Tr(long, m_form)
  spec <- build_hazard_spec(m_form, baseline, dpar_formulas, observation,
                            k_sw = k_sw, costressors = costressors)

  # Build the default prior only when none is supplied, forwarding the override
  # arguments; slope_prior / sd_prior are passed on only when set, so their default
  # strings stay defined in one place (build_hazard_prior).
  if (is.null(prior)) {
    pa <- list(m_form = m_form, baseline = baseline, dpar_formulas = dpar_formulas,
               intercept_priors = intercept_priors, costressors = costressors)
    if (!is.null(slope_prior)) pa$slope_prior <- slope_prior
    if (!is.null(sd_prior))    pa$sd_prior    <- sd_prior
    prior <- do.call(build_hazard_prior, pa)
  }

  args <- list(formula = spec$formula, data = long, family = spec$family, prior = prior)
  # The custom families need their Stan log-mass functions; continuous mode uses the
  # built-in Poisson and needs none.
  if (observation == "interval") args$stanvars <- .cloglog_bern_stanvars
  if (observation == "mixed")    args$stanvars <- .hazard_mixed_stanvars
  args
}

# .hazard_settings() / .fit_settings() -----------------------------------------

# Record the settings a failure-rate model is fitted with, which the reporting
# functions need to evaluate its failure rate exactly as it was fitted: the
# failure-rate and baseline functions, the observation mode, the reference
# temperature (the constant `Tr` column of `long`), the softening steepness actually
# used, and the co-stressors. It also records the ranges of temperature and of
# elapsed exposure time the data cover, which the reporting functions use as their
# default grids, since a model with a memoryless baseline does not keep elapsed
# time in its data, and the individual each row of `long` belongs to, which the
# fitted data do not keep and which cross-validation needs to leave whole
# individuals out. fit_hazard() attaches the record to the fit it returns. The
# steepness is resolved quietly here, since build_hazard_spec() has already warned
# about a `k_sw` it ignores.
.hazard_settings <- function(long, m_form, baseline, observation, k_sw, costressors) {
  list(m_form      = m_form,
       baseline    = baseline,
       observation = .match_observation(observation),
       Tr          = long$Tr[1],
       k_sw        = suppressWarnings(.resolve_k_sw(m_form, k_sw)),
       costressors = costressors,
       Ta_range    = range(long$Ta),
       tau_range   = range(long$tau),
       id          = long$id)
}

# Read the settings fit_hazard() attached to a fit, or stop when there are none: a
# model fitted some other way, whose settings the reporting functions cannot know. A
# model changed afterwards with update() keeps its original record, so refit with
# fit_hazard() instead.
.fit_settings <- function(fit) {
  if (!inherits(fit, "brmsfit"))
    cli::cli_abort("{.arg fit} must be a fitted model ({.cls brmsfit}).",
                   class = "rTFT_not_brmsfit")
  if (is.null(fit$rtft))
    cli::cli_abort(c(
      "{.arg fit} carries no record of the settings it was fitted with.",
      "i" = "Fit the model with {.fn fit_hazard}, which records the settings the reporting functions need."
    ), class = "rTFT_no_settings")
  fit$rtft
}


# fit_hazard() -----------------------------------------------------------------

#' Fit a failure-rate hazard model with brms
#'
#' Assembles the formula and family with [build_hazard_spec], the default prior with
#' [build_hazard_prior] (unless you supply your own), selects the Stan code the custom
#' families need, and hands it all to [brms::brm]. Anything in `...` is forwarded to
#' `brm` (chains, iterations, cores, backend, initial values, and so on).
#'
#' Before fitting, the temperature is checked against the reference `Tr` for the
#' failure rates that take a logarithm or fractional power of `Ta - Tr` (`zeroexp`,
#' `power`, `zeropower`): these are undefined where the temperature reaches `Tr`
#' (`zeroexp`, `zeropower`) or drops below it (`power`), so a run that would hit
#' that is stopped with a clear error rather than a cryptic sampler failure.
#'
#' @param long A data frame, the row-per-interval table from [build_hazard_data].
#' @param m_form,baseline,dpar_formulas,observation As in [build_hazard_spec].
#' @param prior An optional `brmsprior` for full manual control. When `NULL` (the
#'   default) a weakly-informative prior is built with [build_hazard_prior]. A
#'   supplied prior is used as given, in place of that default rather than merged
#'   with it, so any coefficient it leaves out falls back to brms's own default, not
#'   the weakly-informative one. To adjust only some priors and keep the rest, use
#'   `intercept_priors` / `slope_prior` / `sd_prior` instead of this.
#' @param intercept_priors,slope_prior,sd_prior As in [build_hazard_prior],
#'   forwarded when the default prior is built (ignored if you supply `prior`). Each
#'   changes only the priors it names and leaves the rest at their defaults, so a
#'   user who only calls `fit_hazard` can retune without building a prior by hand.
#' @param costressors An optional list of [costressor] objects, as in
#'   [build_hazard_spec], forwarded to both the formula and the default prior. The
#'   same list must have been passed to [build_hazard_data], so `long` carries the
#'   `z_<name>` columns the formula reads. Defaults to none.
#' @param k_sw An optional softening steepness for a threshold function's gate, as
#'   in [build_hazard_spec]. An advanced override: leave `NULL` (the default) to use
#'   the value in [M_catalogue]. A fit made with a custom value must be reported with
#'   the same one.
#' @param ... Further arguments forwarded to [brms::brm]. Any argument it accepts
#'   (e.g. `chains`, `iter`, `cores`, `backend`, `init`), except `formula`, `data`,
#'   `family`, `prior` and `stanvars`, which `fit_hazard` sets itself. See
#'   [brms::brm] for the full set.
#' @return A `brmsfit`, carrying in `$rtft` the settings it was fitted with: the
#'   failure-rate and baseline functions, the observation mode, the reference
#'   temperature `Tr`, the softening steepness, the co-stressors, the ranges of
#'   temperature and elapsed exposure time the data cover, and the individual each
#'   row of `long` belongs to. The reporting functions read them from there, so they
#'   are never re-entered by hand.
#' @examples
#' \dontrun{
#' path <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' ind  <- data.frame(id = 1:20, path_id = "A", event = 1,
#'                    t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long <- build_hazard_data(ind, path, "interval")
#' fit  <- fit_hazard(long, "exp", chains = 4, iter = 2000, cores = 4)
#' }
#' @export
fit_hazard <- function(long, m_form, baseline = "mem", dpar_formulas = list(),
                       observation = c("interval", "continuous", "mixed"), prior = NULL,
                       intercept_priors = list(), slope_prior = NULL, sd_prior = NULL,
                       costressors = list(), k_sw = NULL, ...) {
  args <- .hazard_brm_args(long, m_form, baseline, dpar_formulas, observation, prior,
                           intercept_priors, slope_prior, sd_prior, costressors, k_sw)
  fit <- do.call(brms::brm, c(args, list(...)))
  fit$rtft <- .hazard_settings(long, m_form, baseline, observation, k_sw, costressors)
  fit
}


# ===========================================================================
# 6. Fixed-ramp helpers
# ===========================================================================

# ramp_path() ------------------------------------------------------------------

#' Build a linear ramp temperature trajectory
#'
#' A convenience function for the `path` argument of [build_hazard_data], useful when
#' there is no logged temperature series: when temperature follows a controlled linear
#' ramp, `temp(t) = T_start + rate * t`, the nominal ramp stands in for a logged
#' trace. It is vectorised over its arguments: give several ramps at once and it
#' stacks them into one table keyed by `path_id`.
#'
#' @param path_id,T_start,rate,t_max The ramp specification, one value per ramp:
#'   `path_id` the identifier, `T_start` the start temperature, `rate` in degrees
#'   per time unit, and `t_max` the duration. The four are recycled together, so
#'   equal-length vectors describe several ramps and a length-one argument is reused
#'   across all of them.
#' @param by An optional numeric giving the sampling step in time. When `NULL` (the
#'   default) the path is just its two endpoints, which the linear interpolation in
#'   [build_hazard_data] reproduces exactly. A value such as `by = 5` creates a row
#'   every `by` time units, which then serves as the check schedule in interval mode
#'   or a ready-made integration grid in continuous mode.
#' @return A data frame with `path_id`, `time` and `temp`, ready to pass to
#'   [build_hazard_data] as its `path`.
#' @examples
#' # One ramp: 20 -> 35 degrees over 60 minutes, as a two-point path.
#' ramp_path("A", T_start = 20, rate = 0.25, t_max = 60)
#'
#' # Two ramps at different rates, each sampled every 10 minutes (T_start recycled).
#' ramp_path(c("slow", "fast"), T_start = 20, rate = c(0.2, 0.4),
#'           t_max = c(75, 40), by = 10)
#' @export
ramp_path <- function(path_id, T_start, rate, t_max, by = NULL) {
  purrr::pmap(list(path_id, T_start, rate, t_max), function(pid, Ts, rr, tm) {
    t <- if (is.null(by)) c(0, tm) else sort(unique(c(seq(0, tm, by = by), tm)))
    tibble::tibble(path_id = pid, time = t, temp = Ts + rr * t)
  }) |> purrr::list_rbind()
}


# assess_ramp_linearity() ------------------------------------------------------

#' Assess how linear a logged temperature trace is
#'
#' A pre-fitting diagnostic with two uses. As quality control on the ramp, it fits
#' `temp ~ time` per path and reports the realised start and rate, the fit
#' (`r_squared`), and, given `nominal`, the `rate_bias` against the setpoint, so
#' you can see whether the controller held a linear ramp at the intended rate. And
#' as a check on whether a dense or noisy log can be simplified for fitting: a finely
#' sampled temperature series turns every reading into an at-risk-interval boundary in
#' [build_hazard_data], so a linear one can be replaced by a clean two-point ramp
#' from [ramp_path] (exact under the linear interpolation, free of sensor jitter,
#' and far fewer rows). A faithful series is otherwise best fitted directly; this only
#' judges that faithfulness.
#'
#' As a rule of thumb, a series is linear enough to simplify when its deviations (and
#' especially `max_dev_hot`, the worst deviation in the upper temperature quartile,
#' where the failure rate M(T) is steepest) are small relative to `log(10) / b`, the
#' number of degrees for a tenfold change in an exponential M(T): half that distance
#' already multiplies M(T) by `sqrt(10)`, about threefold.
#'
#' @param path A data frame with `path_id`, `time` and `temp` -- the logged
#'   temperature trace(s).
#' @param nominal An optional named numeric of setpoint rates (names matching
#'   `path_id`). When given, adds `rate_nominal` and `rate_bias` (realised minus
#'   nominal), the controller's rate error.
#' @return A tibble, one row per `path_id`, with `n` (the number of readings), the
#'   realised `T_start` and `rate`, `r_squared`, `rmse`, `max_abs_dev`, and
#'   `max_dev_hot` (the worst absolute deviation in the hot quartile); plus
#'   `rate_nominal` and `rate_bias` when `nominal` is supplied.
#' @examples
#' # A logged temperature series close to a 0.3 degree/min ramp, with a small mid-run wobble.
#' path <- data.frame(path_id = "r1", time = 0:20,
#'                    temp = 20 + 0.3 * (0:20) + c(rep(0, 8), 0.4, 0.5, 0.3, rep(0, 10)))
#' assess_ramp_linearity(path)
#'
#' # With the setpoint rate, to see the controller's rate bias.
#' assess_ramp_linearity(path, nominal = c(r1 = 0.3))
#' @export
assess_ramp_linearity <- function(path, nominal = NULL) {
  miss <- setdiff(c("path_id", "time", "temp"), names(path))
  if (length(miss))
    cli::cli_abort("{.arg path} is missing column{?s} {.field {miss}}.",
                   class = "rTFT_missing_column")

  purrr::imap(split(path, path$path_id), function(d, pid) {
    d   <- d[order(d$time), ]
    m   <- stats::lm(temp ~ time, data = d)
    dev  <- stats::residuals(m)
    ss_t <- sum((d$temp - mean(d$temp))^2)           # total sum of squares
    hot  <- d$temp >= stats::quantile(d$temp, 0.75)  # the upper temperature quartile
    out <- tibble::tibble(
      path_id     = pid,
      n           = nrow(d),
      T_start     = unname(stats::coef(m)[1]),
      rate        = unname(stats::coef(m)[2]),
      # R^2 computed directly rather than via summary.lm(), which warns "essentially
      # perfect fit" on an exactly linear trace (e.g. a ramp_path() output). NA when
      # temp is constant (no variation to explain).
      r_squared   = if (ss_t > 0) 1 - sum(dev^2) / ss_t else NA_real_,
      rmse        = sqrt(mean(dev^2)),
      max_abs_dev = max(abs(dev)),
      max_dev_hot = max(abs(dev[hot]))
    )
    if (!is.null(nominal) && pid %in% names(nominal)) {
      out$rate_nominal <- nominal[[pid]]
      out$rate_bias    <- out$rate - nominal[[pid]]
    }
    out
  }) |> purrr::list_rbind()
}
