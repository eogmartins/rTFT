# Tests for R/discrete-time-hazard-brms.R.
#
# The sections below match the blocks of the source.
#
# Shared across sections: stub_prep() builds the smallest object
# brms::get_dpar() will accept (i.e., a list containing the linear predictor mu
# as a draws-by-observations matrix, tagged with the brmsprep class), so the
# post-processing methods can be tested without fitting a model. Every assertion
# was written after running the function and reading the output, and each
# expected value is derived a second way (from the mathematical definition, not
# from the function under test), so the test checks correctness rather than
# re-recording what the function returned.

stub_prep <- function(mu, y, vint1 = NULL) {
  if (is.null(dim(mu))) mu <- matrix(mu, ncol = 1)
  data <- list(Y = y)
  if (!is.null(vint1)) data$vint1 <- vint1     # per-row known-time flag, for hazard_mixed
  structure(
    list(dpars = list(mu = mu), data = data,
         ndraws = nrow(mu), nobs = ncol(mu)),
    class = "brmsprep"
  )
}

# ===========================================================================
# 0a. Numerically stable cloglog Bernoulli family (interval mode)
# ===========================================================================

# Interval-censored cloglog Bernoulli response ---------------------------------
#
# The family declaration (cloglog_bern_family), its Stan log-mass function
# (.cloglog_bern_stanvars), and the three R methods brms calls after fitting.

# The Stan log-mass definition, written out in R as this family's reference: a
# death in the interval contributes log(1 - exp(-Lambda)), a survived interval
# -Lambda, with Lambda = exp(mu).
ref_loglik <- function(mu, y) {
  lam <- exp(mu)
  if (y == 1) log(1 - exp(-lam)) else -lam
}


test_that("cloglog_bern_family() is the family brms is given", {
  fam <- cloglog_bern_family()
  expect_s3_class(fam, "brmsfamily")
  expect_equal(fam$name, "cloglog_bern")
  expect_equal(fam$dpars, "mu")
  # The identity link is the whole point: the family uses eta directly and
  # forms the log density from it, rather than a probability that saturates to 1.
  expect_equal(fam$link, "identity")
  expect_equal(fam$type, "int")
  # A binary response
  expect_equal(fam$lb$mu, "0")
  expect_equal(fam$ub$mu, "1")
})


test_that(".cloglog_bern_stanvars adds cloglog_bern_lpmf to Stan's functions block", {
  sv <- .cloglog_bern_stanvars
  expect_s3_class(sv, "stanvars")
  expect_equal(sv[[1]]$block, "functions")
  expect_match(sv[[1]]$scode, "real cloglog_bern_lpmf\\(int y, real mu\\)")
  # The survived branch returns -lam directly, never log(1 - p). That is the line
  # that keeps the density finite when the cumulative hazard is large.
  expect_match(sv[[1]]$scode, "return -lam")
})


test_that("log_lik_cloglog_bern reproduces the Stan log-mass and is never positive", {
  set.seed(1)
  mu <- matrix(rnorm(200 * 3), nrow = 200)   # 200 draws, 3 observations
  y  <- c(1, 0, 1)
  prep <- stub_prep(mu, y)

  for (i in seq_len(3)) {
    ll <- log_lik_cloglog_bern(i, prep)
    expect_equal(ll, ref_loglik(mu[, i], y[i]))
    expect_length(ll, 200)
    # A log-probability is never above zero, so this checks the sign directly
    # because a sign error would still return numbers of the right length.
    expect_true(all(ll <= 0))
  }
})


test_that("posterior_epred_cloglog_bern returns the interval failure probability", {
  set.seed(2)
  mu <- matrix(rnorm(50 * 4), nrow = 50)     # 50 draws, 4 observations
  prep <- stub_prep(mu, y = rep(0, 4))       # Y is unused by the expected value

  ep <- posterior_epred_cloglog_bern(prep)
  expect_equal(dim(ep), c(50L, 4L))          # a draws-by-observations matrix
  expect_equal(ep, 1 - exp(-exp(mu)))        # the expected value is p itself
  expect_true(all(ep >= 0 & ep <= 1))
})


test_that("posterior_predict_cloglog_bern draws Bernoulli(p), not the probability", {
  set.seed(3)
  # One observation, many draws sharing the same mu, so the draws should average
  # to p. mu = 0 gives p = 1 - exp(-exp(0)) ~ 0.632.
  ndraw <- 20000
  prep  <- stub_prep(matrix(0, nrow = ndraw, ncol = 1), y = 1)

  pp <- posterior_predict_cloglog_bern(1L, prep)
  expect_length(pp, ndraw)
  expect_true(all(pp %in% c(0, 1)))          # draws, not probabilities
  # And they are genuinely Bernoulli(p): the sample mean tracks p.
  expect_equal(mean(pp), 1 - exp(-exp(0)), tolerance = 0.02)
})


test_that("the family stays finite where the naive inverse-cloglog underflows", {
  # This is the reason the custom family exists. Here, we test a survived interval
  # whose linear predictor is large (routine near CTmax), so the cumulative hazard
  # is huge. It should return finite with the custom family and underflow with
  # the naive inverse-cloglog.
  prep <- stub_prep(matrix(40, nrow = 1, ncol = 1), y = 0)   # lam = exp(40)

  ll <- log_lik_cloglog_bern(1L, prep)
  expect_true(is.finite(ll))                 # the custom family returns -lam
  expect_equal(ll, -exp(40))

  # The route brms's bernoulli(cloglog) would take: form p, then log(1 - p).
  p_naive  <- 1 - exp(-exp(40))              # rounds to exactly 1
  expect_true(is.infinite(log(1 - p_naive))) # log(0) = -Inf, the failure
})


# ===========================================================================
# 0b. Unified family for mixed interval-censored + known-time failures
# ===========================================================================

# Mixed interval-censored and known-time response ------------------------------
#
# The family declaration (hazard_mixed_family), its Stan log-mass and density
# function (.hazard_mixed_stanvars), and the three R methods brms calls after
# fitting.

# The Stan log-mass and density definition, written out in R as this family's
# reference. Three branches, keyed on the survived indicator y and the known-time
# flag kt: a survived interval contributes -Lambda, an interval-censored failure
# log(1 - exp(-Lambda)), and a known-time failure the log density mu - Lambda,
# with Lambda = exp(mu).
ref_loglik_mixed <- function(mu, y, kt) {
  lam <- exp(mu)
  if (y == 0)  return(-lam)
  if (kt == 0) return(log(1 - exp(-lam)))
  mu - lam
}


test_that("hazard_mixed_family() is the family brms is given", {
  fam <- hazard_mixed_family()
  expect_s3_class(fam, "brmsfamily")
  expect_equal(fam$name, "hazard_mixed")
  expect_equal(fam$dpars, "mu")
  expect_equal(fam$link, "identity")
  expect_equal(fam$type, "int")
  expect_equal(fam$lb$mu, "0")
  expect_equal(fam$ub$mu, "1")
  # The per-row known-time flag travels as an integer vint, which is what lets one
  # fit contain both kinds of failure.
  expect_true(any(grepl("vint1", fam$vars)))
})


test_that(".hazard_mixed_stanvars adds hazard_mixed_lpmf to Stan's functions block", {
  sv <- .hazard_mixed_stanvars
  expect_s3_class(sv, "stanvars")
  expect_equal(sv[[1]]$block, "functions")
  expect_match(sv[[1]]$scode, "real hazard_mixed_lpmf\\(int y, real mu, int vint1\\)")
  # the three branches: survived, interval-censored, known-time
  expect_match(sv[[1]]$scode, "return -lam")
  expect_match(sv[[1]]$scode, "log1m_exp\\(-lam\\)")
  expect_match(sv[[1]]$scode, "return mu - lam")
})


test_that("log_lik_hazard_mixed reproduces the survived, interval, and known-time branches", {
  set.seed(4)
  mu <- matrix(rnorm(200 * 3), nrow = 200)   # 200 draws, 3 observations
  y  <- c(0, 1, 1)                           # survived, failure, failure
  kt <- c(0, 0, 1)                           # kt is ignored for the survived row
  prep <- stub_prep(mu, y, vint1 = kt)

  for (i in seq_len(3)) {
    ll <- log_lik_hazard_mixed(i, prep)
    expect_equal(ll, ref_loglik_mixed(mu[, i], y[i], kt[i]))
    expect_length(ll, 200)
  }
})


test_that("vint1 selects the failure term, and is ignored when the individual survived", {
  set.seed(5)
  mu  <- matrix(rnorm(50), nrow = 50, ncol = 1)   # 50 draws, 1 observation
  lam <- exp(mu[, 1])

  # On a failure row (y = 1) the flag picks interval vs known-time, and they differ.
  ll_interval <- log_lik_hazard_mixed(1L, stub_prep(mu, y = 1, vint1 = 0))
  ll_known    <- log_lik_hazard_mixed(1L, stub_prep(mu, y = 1, vint1 = 1))
  expect_equal(ll_interval, log(1 - exp(-lam)))    # interval-censored term
  expect_equal(ll_known,    mu[, 1] - lam)         # known-time term
  expect_false(isTRUE(all.equal(ll_interval, ll_known)))
  # The interval-censored term is exactly the cloglog family's failure term.
  expect_equal(ll_interval, log_lik_cloglog_bern(1L, stub_prep(mu, y = 1)))

  # On a survived row (y = 0) the flag makes no difference: both give -lam.
  ll_surv0 <- log_lik_hazard_mixed(1L, stub_prep(mu, y = 0, vint1 = 0))
  ll_surv1 <- log_lik_hazard_mixed(1L, stub_prep(mu, y = 0, vint1 = 1))
  expect_equal(ll_surv0, ll_surv1)
  expect_equal(ll_surv0, -lam)
})


test_that("posterior_epred_hazard_mixed returns the interval failure probability", {
  set.seed(6)
  mu <- matrix(rnorm(50 * 4), nrow = 50)     # 50 draws, 4 observations
  prep <- stub_prep(mu, y = rep(0, 4))       # epred reads neither Y nor vint1

  ep <- posterior_epred_hazard_mixed(prep)
  expect_equal(dim(ep), c(50L, 4L))
  expect_equal(ep, 1 - exp(-exp(mu)))
  expect_true(all(ep >= 0 & ep <= 1))
})


test_that("posterior_predict_hazard_mixed draws Bernoulli(p)", {
  set.seed(7)
  prep <- stub_prep(matrix(0, nrow = 20000, ncol = 1), y = 1)
  pp <- posterior_predict_hazard_mixed(1L, prep)
  expect_length(pp, 20000)
  expect_true(all(pp %in% c(0, 1)))
  expect_equal(mean(pp), 1 - exp(-exp(0)), tolerance = 0.02)
})


# ===========================================================================
# 1. Catalogues and log M(T) / log h0(tau) expressions
# ===========================================================================

# The catalogues (M_catalogue, h0_catalogue) and the expression-string objects
# (brms_log_M, brms_log_h0). The correctness checks evaluate each stored string in
# R and compare it to the closed form written out independently from the source of
# the model, not from the string under test: Kingsolver & Umbanhowar (2018) Table 1
# for M(T), and the standard survival-analysis baseline hazards for h0(tau). This
# checks the algebra the strings encode; it is not an end-to-end check that brms
# and Stan translate them identically (they use the same functions exp, log,
# log1p, abs, ^m, so R evaluation is a faithful proxy for the math).

# Keys and structure -----------------------------------------------------------

test_that("brms_log_M keys and stan_ids line up with M_catalogue", {
  expect_identical(names(brms_log_M), names(M_catalogue))
  expect_length(M_catalogue, 9L)
  # stan_id is the integer Stan dispatches on: the nine distinct codes 1..9
  expect_identical(unname(vapply(M_catalogue, `[[`, integer(1), "stan_id")), 1:9)
})


test_that("brms_log_h0 keys and stan_ids line up with h0_catalogue", {
  expect_identical(names(brms_log_h0), names(h0_catalogue))
  expect_length(h0_catalogue, 4L)
  expect_identical(unname(vapply(h0_catalogue, `[[`, integer(1), "stan_id")), 1:4)
})


test_that("every parameter a function declares appears in its expression string", {
  # M_catalogue$dpars are the parameters the string must use; each must occur in it.
  for (k in names(M_catalogue)) {
    dp <- M_catalogue[[k]]$dpars
    expect_true(all(vapply(dp, grepl, logical(1), x = brms_log_M[[k]], fixed = TRUE)),
                info = k)
  }
  for (k in names(h0_catalogue)) {
    dp <- h0_catalogue[[k]]$dpars                 # mem has none: all() of nothing is TRUE
    expect_true(all(vapply(dp, grepl, logical(1), x = brms_log_h0[[k]], fixed = TRUE)),
                info = k)
  }
})


test_that("k_sw is nonzero exactly for the functions with a softening gate", {
  # The threshold functions are those whose string carries a log1p gate, and only
  # those carry a nonzero softening steepness. This ties the k_sw field to the string.
  k_sw     <- vapply(M_catalogue, `[[`, numeric(1), "k_sw")
  has_gate <- grepl("log1p", unlist(brms_log_M[names(M_catalogue)]))
  expect_identical(unname(k_sw > 0), unname(has_gate))
})


# brms_log_M against Kingsolver & Umbanhowar Table 1 ---------------------------

test_that("brms_log_M reproduces the Kingsolver-Umbanhowar Table 1 forms above Tc", {
  a <- 2.82e-7; b <- 6.13; cc <- 1e-4; Tc <- 36; Tr <- 30; k_sw <- 50; Ta <- 42
  env <- list(loga = log(a), logb = log(b), logc = log(cc),
              Tc = Tc, Tr = Tr, k_sw = k_sw, Ta = Ta)  # above Tc and Tr
  ev <- function(s) exp(eval(parse(text = s), env))

  # KU Table 1 closed forms; threshold entries take the above-Tc branch
  u <- Ta - Tr; v <- Ta - Tc
  ku <- list(
    const               = a,
    const_threshold     = a,
    exp                 = a * exp(b * u),
    exp_threshold       = a * exp(b * v),
    zeroexp             = a * exp(b * u) - a,
    zeroexp_threshold   = a * exp(b * v) - a,
    power               = cc + a * u^b,
    zeropower           = a * u^b,
    zeropower_threshold = a * v^b
  )
  for (k in names(brms_log_M)) {
    expect_equal(ev(brms_log_M[[k]]), ku[[k]], info = k)
  }
})


test_that("no failure-rate or baseline function puts a number in scientific notation into the Stan program", {
  # brms spaces out the sign of an exponent when it writes the Stan program, so
  # `1e-12` reaches Stan as `1e - 12`, which does not compile. The expression
  # strings must therefore hold no such number,
  strings <- c(unlist(brms_log_M), unlist(brms_log_h0))
  expect_identical(names(strings)[grepl("[0-9][eE][-+]?[0-9]", strings)], character(0))
  # and none may appear, broken, in the program of any function with any baseline.
  ramp <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
  ind  <- data.frame(id = 1:6, path_id = "A", event = 1, t_event = c(30, 35, 40, 42, 45, 50))
  programs <- purrr::map(names(M_catalogue), \(m) purrr::map_chr(names(h0_catalogue), \(b) {
    Tr   <- if (m %in% c("zeroexp", "power", "zeropower")) 19
    long <- suppressMessages(build_hazard_data(ind, ramp, "continuous", Tr = Tr))
    a    <- .hazard_brm_args(long, m, b, observation = "continuous")
    brms::make_stancode(a$formula, data = a$data, family = a$family, prior = a$prior,
                        stanvars = a$stanvars)
  }))
  programs <- unlist(programs)
  expect_length(programs, length(M_catalogue) * length(h0_catalogue))   # all 36 written
  expect_false(any(grepl("[0-9][eE] [-+] [0-9]", programs)))
})


test_that("brms_log_M threshold guards stay finite below and at Tc; zero forms vanish at Tr", {
  a <- 2.82e-7; b <- 6.13; cc <- 1e-4; Tc <- 36; Tr <- 30; k_sw <- 50
  ev <- function(s, Ta) {
    env <- list(loga = log(a), logb = log(b), logc = log(cc),
                Tc = Tc, Tr = Tr, k_sw = k_sw, Ta = Ta)
    exp(eval(parse(text = s), env))
  }
  # abs(...) + the small constant keep the below-Tc and at-Tc logarithms finite
  # (no NaN, no -Inf), so the gate can suppress that region cleanly.
  for (k in c("const_threshold", "exp_threshold", "zeroexp_threshold", "zeropower_threshold")) {
    expect_true(is.finite(ev(brms_log_M[[k]], Tc - 2)), info = paste(k, "below Tc"))
    expect_true(is.finite(ev(brms_log_M[[k]], Tc)),     info = paste(k, "at Tc"))
  }
  # the "zero" forms are exactly zero at the starting temperature
  expect_equal(ev(brms_log_M[["zeroexp"]],   Tr), 0)
  expect_equal(ev(brms_log_M[["zeropower"]], Tr), 0)
})


# brms_log_h0 against the baseline hazards ------------------------------------

test_that("brms_log_h0 reproduces the baseline hazards, and llogis is the hazard not the density", {
  shape <- 1.8; scale <- 20; gam <- 0.05; tau <- 12; tau_ref <- 45
  env <- list(logshape = log(shape), logscale = log(scale), gam = gam,
              tau = tau, log_tau = log(tau), tau_ref = tau_ref)
  ev <- function(s) exp(eval(parse(text = s), env))

  z <- (tau / scale)^shape
  # the Weibull scale is the reference time, not a parameter
  h_weibull  <- (shape / tau_ref) * (tau / tau_ref)^(shape - 1)
  h_llogis   <- (shape / scale) * (tau / scale)^(shape - 1) / (1 + z)   # hazard
  h_gompertz <- exp(gam * tau)                                          # amplitude in loga
  f_llogis   <- h_llogis / (1 + z)                                      # density (factor-2 error)

  expect_equal(ev(brms_log_h0[["weibull"]]),  h_weibull)
  expect_false(grepl("logscale", brms_log_h0[["weibull"]]))
  expect_equal(h0_catalogue$weibull$dpars, "logshape")
  expect_equal(h0_catalogue$llogis$dpars, c("logshape", "logscale"))
  expect_equal(ev(brms_log_h0[["llogis"]]),   h_llogis)
  expect_equal(ev(brms_log_h0[["gompertz"]]), h_gompertz)
  # mem is the empty string: a memoryless baseline adds nothing to the predictor
  expect_identical(brms_log_h0[["mem"]], "")
  # the llogis string is the hazard (matches h), not the density (differs from f)
  expect_false(isTRUE(all.equal(ev(brms_log_h0[["llogis"]]), f_llogis)))
})


# ===========================================================================
# 2. Build the row-per-interval table (build_hazard_data)
# ===========================================================================

# build_hazard_data() reshapes one row per individual into one row per at-risk
# interval. The checks below use a single linear temperature ramp, so every derived
# quantity has a closed form written independently of the function: for the ramp
# temp(t) = 20 + 0.25 t (20 -> 35 over 60 minutes) the temperature at an interval
# midpoint is 20 + 0.25 * mid, and with checks every 10 minutes each interval is 10
# wide with midpoints at the odd multiples of 5.

ramp_traj   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
ramp_temp   <- function(t) 20 + 0.25 * t                # the ramp
ramp_checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))


test_that("interval mode gives one row per at-risk interval up to the failure check", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)

  # Failure found at the 30-min check -> intervals [0,10], [10,20], [20,30].
  expect_equal(nrow(out), 3L)
  expect_equal(out$tau, c(5, 15, 25))                   # interval midpoints
  expect_equal(out$exposure, rep(10, 3))
  expect_equal(out$Ta, ramp_temp(c(5, 15, 25)))         # from the full ramp, at midpoints
  expect_equal(out$failed, c(0L, 0L, 1L))               # only the interval it falls in
})


test_that("a failure found at a check lands in the interval ending there, not the next", {
  # t_upper is exactly the 30-min check. The -1e-9 nudge keeps the failure in
  # [20,30], the interval ENDING at 30, rather than [30,40], the one starting there.
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)
  fail_row <- out[out$failed == 1, ]
  expect_equal(nrow(fail_row), 1L)
  expect_equal(fail_row$tau, 25)                        # midpoint of [20,30]
})


test_that("a right-censored survivor never fails and runs to its last check", {
  ind <- data.frame(id = 2, path_id = "A", event = 0, t_upper = 60)
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)
  expect_equal(nrow(out), 6L)                           # [0,10] .. [50,60]
  expect_true(all(out$failed == 0L))                    # event = 0 -> no failure row
  expect_equal(max(out$tau), 55)                        # last interval midpoint
})


test_that("the reference time is the median time of failure unless given", {
  # failures found at the 20, 30 and 50 minute checks, and a survivor seen to 60
  four <- data.frame(id = 1:4, path_id = "A", event = c(1, 1, 1, 0), t_upper = c(20, 30, 50, 60))
  expect_equal(unique(build_hazard_data(four, ramp_traj, "interval",
                                        checks = ramp_checks, grid_dT = Inf)$tau_ref), 30)
  expect_equal(.typical_failure_time(build_hazard_data(four, ramp_traj, "interval",
                                                       checks = ramp_checks, grid_dT = Inf)), 30)
  # with no failures, the longest exposure stands in
  none <- build_hazard_data(transform(four, event = 0), ramp_traj, "interval",
                            checks = ramp_checks, grid_dT = Inf)
  expect_equal(unique(none$tau_ref), 60)
  ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 50))
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)
  expect_equal(unique(out$tau_ref), 30)                 # the one failure, not the survivor
  # a known-time failure: its exact time, not the midpoint of its last interval
  ind_c <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  expect_equal(unique(build_hazard_data(ind_c, ramp_traj, "continuous",
                                        checks = ramp_checks, grid_dT = Inf)$tau_ref), 42)
  given <- build_hazard_data(ind, ramp_traj, "interval", tau_ref = 60, checks = ramp_checks, grid_dT = Inf)
  expect_equal(unique(given$tau_ref), 60)
  for (bad in list(0, -5, c(30, 60), "60", NA_real_, Inf))
    expect_error(build_hazard_data(ind, ramp_traj, "interval", tau_ref = bad,
                                   checks = ramp_checks, grid_dT = Inf), class = "rTFT_bad_tau_ref")
})


test_that("a table with no reference time is refused", {
  ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 50))
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)
  no_ref <- out[setdiff(names(out), "tau_ref")]
  expect_error(.hazard_brm_args(no_ref, "exp", "weibull"), class = "rTFT_no_tau_ref")
  # the other baselines need it too, to scale their default priors
  expect_error(.hazard_brm_args(no_ref, "exp", "mem"), class = "rTFT_no_tau_ref")
  expect_no_error(.hazard_brm_args(out, "exp", "weibull"))
})


test_that("the output columns and their derived quantities are correct", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)

  expect_setequal(names(out),
                  c("id", "path_id", "tau", "exposure", "Ta", "failed",
                    "Tr", "log_tau", "tau_ref", "log_exposure"))
  expect_equal(out$log_tau, log(out$tau))
  expect_equal(out$log_exposure, log(out$exposure))
  expect_false("exact" %in% names(out))                 # exact is mixed-mode only
})


# The time grid --------------------------------------------------------------
# `ramp_traj` warms from 20 to 35 degrees over 60 minutes (0.25 degrees per minute)
# and `ramp_checks` are every 10 minutes, so each interval between checks spans 2.5
# degrees: 25 parts at the default of 0.1 degrees per interval.

test_that(".subdivide cuts each interval into equal parts, in order", {
  cut <- .subdivide(ts = c(0, 10), te = c(10, 30), parts = c(2, 4))
  expect_equal(cut$ts, c(0, 5, 10, 15, 20, 25))
  expect_equal(cut$te, c(5, 10, 15, 20, 25, 30))
  whole <- .subdivide(c(0, 10), c(10, 30), c(1, 1))
  expect_equal(whole, list(ts = c(0, 10), te = c(10, 30)))
})


test_that(".temp_span is the range of temperature within an interval, readings inside it included", {
  ramp <- stats::approxfun(c(0, 60), c(20, 35), rule = 2)
  expect_equal(.temp_span(ramp, c(0, 60), ts = c(0, 10), te = c(10, 30)), c(2.5, 5))
  # a trace that peaks between the two ends of an interval: the ends alone would
  # say it did not change
  peak <- stats::approxfun(c(0, 5, 10), c(20, 24, 20), rule = 2)
  expect_equal(.temp_span(peak, c(0, 5, 10), ts = 0, te = 10), 4)
  expect_equal(.temp_span(peak, c(0, 5, 10), ts = c(0, 5), te = c(5, 10)), c(4, 4))
  flat <- stats::approxfun(c(0, 60), c(25, 25), rule = 2)
  expect_equal(.temp_span(flat, c(0, 60), 0, 60), 0)
})


test_that("intervals a fish survived are divided by temperature, and the one it was found failed in is not", {
  ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
  expect_message(
    expect_warning(out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks),
                   class = "rTFT_wide_failure_interval"),
    class = "rTFT_grid_refined")
  one <- out[out$id == 1, ]; two <- out[out$id == 2, ]
  # fish 1: two survived intervals of 25 parts each, then the interval of failure whole
  expect_equal(nrow(one), 2L * 25L + 1L)
  expect_equal(one$exposure, c(rep(0.4, 50), 10))
  expect_equal(one$failed, c(rep(0L, 50), 1L))
  expect_equal(one$tau[51], 25)                           # the midpoint of [20, 30]
  # fish 2 survived all six intervals: every one is divided
  expect_equal(nrow(two), 6L * 25L)
  expect_true(all(two$failed == 0L))
  # the time at risk is what it was, and no row but the failure's spans more than 0.1
  expect_equal(sum(one$exposure), 30)
  expect_equal(sum(two$exposure), 60)
  expect_true(all(0.25 * out$exposure[out$failed == 0] <= 0.1 + 1e-9))
  expect_equal(two$Ta, 20 + 0.25 * two$tau)
})


test_that("a failure at a known time is closed with a short interval ending at it", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  out <- suppressMessages(build_hazard_data(ind, ramp_traj, "continuous", checks = ramp_checks))
  n <- nrow(out)
  expect_equal(out$failed, c(rep(0L, n - 1), 1L))
  expect_equal(out$exposure[n], 1e-4 * 42)                # the closing interval
  expect_equal(out$tau[n] + out$exposure[n] / 2, 42)      # which ends at the failure
  expect_equal(out$Ta[n], 20 + 0.25 * 42, tolerance = 1e-4)   # the temperature of the failure
  expect_equal(sum(out$exposure), 42)
  expect_true(all(0.25 * out$exposure <= 0.1 + 1e-9))
  # four whole intervals of 25 parts, the last (40 to 42, half a degree) in 5, and
  # the closing one
  expect_equal(n, 4L * 25L + 5L + 1L)
  # a survivor has no failure to close at
  alive <- suppressMessages(build_hazard_data(transform(ind, event = 0), ramp_traj,
                                              "continuous", checks = ramp_checks))
  expect_equal(nrow(alive), n - 1L)
  expect_true(all(alive$failed == 0L))
})


test_that("with no checks and no grid, a known-time trajectory is still followed closely", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  # a ramp given by its two ends: 10.5 degrees to the failure, so 105 parts
  out <- suppressMessages(build_hazard_data(ind, ramp_traj, "continuous"))
  expect_equal(nrow(out), 105L + 1L)
  expect_true(all(0.25 * out$exposure <= 0.1 + 1e-9))
  # a constant temperature gives no step in temperature: at least 20 intervals
  flat <- data.frame(path_id = "A", time = c(0, 60), temp = c(25, 25))
  still <- suppressMessages(build_hazard_data(ind, flat, "continuous"))
  expect_equal(nrow(still), 20L + 1L)
  expect_equal(still$exposure[1:19], rep(42 / 20, 19))
  expect_equal(sum(still$exposure), 42)
})


test_that("grid_dT sets the step, and Inf keeps the intervals as given", {
  ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
  # as given: one row per interval between checks, and nothing said
  expect_silent(as_given <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                                              grid_dT = Inf))
  expect_equal(nrow(as_given), 3L + 6L)
  # half a degree per interval: 5 parts where 0.1 gave 25
  half <- suppressWarnings(suppressMessages(
    build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = 0.5)))
  expect_equal(nrow(half), (2L * 5L + 1L) + 6L * 5L)
  # a known-time failure as given: no closing interval either
  ind_c <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  expect_silent(plain <- build_hazard_data(ind_c, ramp_traj, "continuous", checks = ramp_checks,
                                           grid_dT = Inf))
  expect_equal(plain$exposure, c(10, 10, 10, 10, 2))
  for (bad in list(0, -0.1, c(0.1, 0.2), "0.1", NA_real_))
    expect_error(build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                                   grid_dT = bad), class = "rTFT_bad_grid_dT")
})


test_that("nothing is divided, said or warned when temperature does not change", {
  flat <- data.frame(path_id = "A", time = c(0, 60), temp = c(25, 25))
  ind  <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
  expect_silent(out <- build_hazard_data(ind, flat, "interval", checks = ramp_checks))
  expect_equal(nrow(out), 3L + 6L)
})


test_that("the warning counts only failures found in an interval that spans much temperature", {
  # checks every 2 minutes span half a degree, the limit: no warning; every 10, 2.5
  fine <- data.frame(path_id = "A", time = seq(0, 60, by = 2))
  ind  <- data.frame(id = 1:3, path_id = "A", event = c(1, 1, 0), t_upper = c(30, 40, 60))
  expect_no_warning(suppressMessages(build_hazard_data(ind, ramp_traj, "interval", checks = fine)))
  expect_warning(suppressMessages(build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks)),
                 "2 individuals", class = "rTFT_wide_failure_interval")
})


test_that("Tr defaults to the median temperature at failure and can be overridden", {
  # the ramp is 20 + 0.25 t; failures found at the 30, 40 and 60 minute checks fall
  # in intervals whose midpoints are at 26.25, 28.75 and 33.75 degrees
  ind <- data.frame(id = 1:4, path_id = "A", event = c(1, 1, 1, 0), t_upper = c(30, 40, 60, 60))
  out_def <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)
  expect_true(all(out_def$Tr == 28.75))                 # the survivor does not count
  expect_equal(.typical_failure_temp(out_def), 28.75)
  # with no failures, the median over every interval at risk stands in
  none <- build_hazard_data(transform(ind, event = 0), ramp_traj, "interval",
                            checks = ramp_checks, grid_dT = Inf)
  expect_equal(unique(none$Tr), stats::median(none$Ta))
  # a zero-type rate cannot use that default: the reference must be given
  expect_error(.hazard_brm_args(out_def, "zeroexp"), class = "rTFT_Tr_below_data")
  cold <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, Tr = 19, grid_dT = Inf)
  expect_no_error(.hazard_brm_args(cold, "zeroexp"))

  out_15 <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, Tr = 15, grid_dT = Inf)
  expect_true(all(out_15$Tr == 15))                     # a constant column at the set value
})


test_that("covariates are replicated onto every row and moved to the front", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30, mass = 12.1)
  out <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                           covariates = "mass", grid_dT = Inf)
  expect_true(all(out$mass == 12.1))                    # the per-individual value, on each row
  expect_equal(names(out)[1:2], c("id", "mass"))        # id then covariates lead the table
})


test_that("continuous mode truncates the final interval at the exact failure time", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  out <- build_hazard_data(ind, ramp_traj, "continuous", checks = ramp_checks, grid_dT = Inf)

  # nodes 0,10,20,30,40,42 -> five intervals, the last [40,42].
  expect_equal(nrow(out), 5L)
  expect_equal(out$exposure, c(10, 10, 10, 10, 2))      # final interval is 2 wide
  expect_equal(out$tau, c(5, 15, 25, 35, 41))           # last midpoint is (40 + 42)/2
  expect_equal(out$failed, c(0L, 0L, 0L, 0L, 1L))
})


test_that("grid_dt integrates a known-time failure on a regular grid ending at the failure time", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  out <- build_hazard_data(ind, ramp_traj, "continuous", grid_dt = 5, grid_dT = Inf)

  # seq(0, 42, 5) = 0,5,..,40, plus the exact 42 -> nine intervals, the last [40,42].
  expect_equal(nrow(out), 9L)
  expect_equal(out$exposure, c(rep(5, 8), 2))
  expect_equal(out$Ta, ramp_temp(out$tau))              # still the true midpoint temperatures
})


test_that("mixed mode sets exact per obs_type and keeps a single `failed` response", {
  ind <- data.frame(id = 1:2, path_id = "A", event = 1,
                    obs_type = c("interval", "continuous"),
                    t_upper = c(30, NA), t_event = c(NA, 42))
  out <- build_hazard_data(ind, ramp_traj, "mixed", checks = ramp_checks, grid_dT = Inf)

  expect_true("exact" %in% names(out))
  expect_true("failed" %in% names(out))                 # one response name in every mode
  expect_true(all(out$exact[out$id == 1] == 0L))        # the interval-censored individual
  expect_true(all(out$exact[out$id == 2] == 1L))        # the known-time individual
})


test_that("temperature comes from the full series and is held flat past its ends", {
  # A logger finer than the checks: a flat leg to t = 25, then a ramp, so the
  # midpoint temperature is the interpolant of the readings, not a check-to-check line.
  fine <- data.frame(path_id = "A", time = c(0, 25, 60), temp = c(20, 20, 35))
  ind  <- data.frame(id = 1, path_id = "A", event = 0, t_upper = 60)
  out  <- build_hazard_data(ind, fine, "interval", checks = ramp_checks, grid_dT = Inf)
  expect_equal(out$Ta[out$tau == 5], 20)                                  # on the flat leg
  expect_equal(out$Ta[out$tau == 35], 20 + (35 - 20) / (60 - 25) * (35 - 25))  # on the ramp

  # Flat extrapolation: a check beyond the last reading still gets a finite Ta.
  short <- data.frame(path_id = "A", time = c(0, 40), temp = c(20, 30))  # series ends at 40
  out2  <- build_hazard_data(ind, short, "interval", checks = ramp_checks, grid_dT = Inf)  # checks reach 60
  expect_true(all(is.finite(out2$Ta)))
  expect_equal(max(out2$Ta), 30)                        # held flat at the last reading
})


test_that("build_hazard_data errors loudly, with cli condition classes, on malformed input", {
  ok_ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)

  # an unknown observation mode
  expect_error(build_hazard_data(ok_ind, ramp_traj, observation = "nope", grid_dT = Inf),
               class = "rTFT_bad_observation")
  # a missing required input column (path lacks temp)
  bad_path <- data.frame(path_id = "A", time = c(0, 60))
  expect_error(build_hazard_data(ok_ind, bad_path, "interval", grid_dT = Inf),
               class = "rTFT_missing_column")
  # the non-mixed timing column absent (interval needs t_upper)
  no_time <- data.frame(id = 1, path_id = "A", event = 1)
  expect_error(build_hazard_data(no_time, ramp_traj, "interval", grid_dT = Inf),
               class = "rTFT_missing_column")
  # a non-finite timing value
  na_time <- data.frame(id = 1, path_id = "A", event = 1, t_upper = NA_real_)
  expect_error(build_hazard_data(na_time, ramp_traj, "interval", grid_dT = Inf),
               class = "rTFT_nonfinite_timing")
})


test_that("mixed mode guards obs_type and its per-level timing columns", {
  # obs_type missing entirely
  m1 <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  expect_error(build_hazard_data(m1, ramp_traj, "mixed", grid_dT = Inf),
               class = "rTFT_missing_obs_type")
  # an unrecognised obs_type level
  m2 <- data.frame(id = 1, path_id = "A", event = 1, obs_type = "typo", t_upper = 30)
  expect_error(build_hazard_data(m2, ramp_traj, "mixed", grid_dT = Inf),
               class = "rTFT_bad_obs_type")
  # a valid level whose timing column is absent
  m3 <- data.frame(id = 1, path_id = "A", event = 1, obs_type = "continuous")
  expect_error(build_hazard_data(m3, ramp_traj, "mixed", grid_dT = Inf),
               class = "rTFT_missing_timing_col")
  # a valid level whose timing value is non-finite
  m4 <- data.frame(id = 1, path_id = "A", event = 1, obs_type = "interval",
                   t_upper = NA_real_)
  expect_error(build_hazard_data(m4, ramp_traj, "mixed", grid_dT = Inf),
               class = "rTFT_nonfinite_timing")
})


# ===========================================================================
# 3. Formula and family builder (build_hazard_spec)
# ===========================================================================

# build_hazard_spec() assembles a nonlinear brms formula and a family without
# fitting, so the checks read the pieces off the returned object: the main formula's
# right-hand side (the linear predictor), the per-parameter sub-formulas, and the
# family. The expected linear predictors are written out by hand from the KU and
# baseline strings, so the test checks the assembly -- paste order, log_exposure, the
# k_sw substitution, the response and the vint() flag -- rather than re-running the
# function's own logic.

# The RHS of the assembled main formula, as a language object (spacing-independent),
# and a family's name whether it is a custom family (`name`) or a built-in (`family`).
main_rhs <- function(ff) ff$formula$formula[[3]]
main_lhs <- function(ff) deparse1(ff$formula$formula[[2]])
fam_name <- function(ff) if (!is.null(ff$family$name)) ff$family$name else ff$family$family


test_that("the linear predictor is log h0 + log M + log_exposure, in that order", {
  # mem contributes no baseline term, so the predictor is log M + log_exposure
  eta_exp_mem <- "loga + exp(logb) * (Ta - Tr) + log_exposure"
  expect_equal(main_rhs(build_hazard_spec("exp", "mem")), str2lang(eta_exp_mem))

  # a weibull baseline prepends its log h0; exp_threshold contributes its gate with
  # the catalogue k_sw = 50 substituted in for the placeholder
  eta_thr_wb <- paste(
    "logshape - log(tau_ref) + (exp(logshape) - 1) * (log_tau - log(tau_ref))",
    "+ loga + exp(logb) * (Ta - Tc) - log1p(exp(-50 * (Ta - Tc)))",
    "+ log_exposure")
  expect_equal(main_rhs(build_hazard_spec("exp_threshold", "weibull")),
               str2lang(eta_thr_wb))
})


test_that("a memoryless baseline contributes no term to the predictor", {
  d_mem <- deparse1(main_rhs(build_hazard_spec("exp", "mem")))
  expect_false(grepl("log_tau", d_mem, fixed = TRUE))     # no baseline-hazard term
  expect_true(grepl("log_exposure", d_mem, fixed = TRUE)) # but still the exposure term
  # a weibull baseline, by contrast, brings in a log_tau term
  expect_true(grepl("log_tau", deparse1(main_rhs(build_hazard_spec("exp", "weibull"))),
                    fixed = TRUE))
})


test_that("k_sw defaults to the catalogue value, and an override substitutes for it", {
  d <- deparse1(main_rhs(build_hazard_spec("exp_threshold")))
  expect_true(grepl("-50 * (Ta - Tc)", d, fixed = TRUE))  # catalogue default
  expect_false(grepl("k_sw", d, fixed = TRUE))            # the placeholder is gone

  o <- deparse1(main_rhs(build_hazard_spec("exp_threshold", k_sw = 20)))
  expect_true(grepl("-20 * (Ta - Tc)", o, fixed = TRUE))  # the override
  expect_false(grepl("-50", o, fixed = TRUE))
})


test_that("k_sw on a threshold-less function is ignored, with a warning", {
  expect_warning(build_hazard_spec("exp", k_sw = 20), class = "rTFT_k_sw_ignored")
  f <- suppressWarnings(build_hazard_spec("exp", k_sw = 20))
  expect_false(grepl("20", deparse1(main_rhs(f)), fixed = TRUE))  # predictor unchanged
})


test_that("a non-positive or non-scalar k_sw errors", {
  expect_error(build_hazard_spec("exp_threshold", k_sw = -5), class = "rTFT_bad_k_sw")
  expect_error(build_hazard_spec("exp_threshold", k_sw = 0),  class = "rTFT_bad_k_sw")
  expect_error(build_hazard_spec("exp_threshold", k_sw = c(10, 20)), class = "rTFT_bad_k_sw")
})


test_that("observation sets the response LHS and the family", {
  fi <- build_hazard_spec("exp", observation = "interval")
  expect_equal(fi$observation, "interval")
  expect_equal(fam_name(fi), "cloglog_bern")
  expect_equal(main_lhs(fi), "failed")

  fc <- build_hazard_spec("exp", observation = "continuous")
  expect_equal(fam_name(fc), "poisson")           # Poisson piecewise-exponential
  expect_equal(fc$family$link, "log")
  expect_equal(main_lhs(fc), "failed")

  fm <- build_hazard_spec("exp", observation = "mixed")
  expect_equal(fam_name(fm), "hazard_mixed")
  expect_equal(main_lhs(fm), "failed | vint(exact)")  # the per-row known-time flag
})


test_that("each needed parameter gets a sub-formula, defaulting to intercept-only", {
  f <- build_hazard_spec("exp_threshold", "weibull")
  # exp_threshold needs loga, logb, Tc; weibull adds logshape (its scale is fixed)
  expect_setequal(names(f$formula$pforms), c("loga", "logb", "Tc", "logshape"))
  # the log-logistic keeps its scale
  expect_setequal(names(build_hazard_spec("exp", "llogis")$formula$pforms),
                  c("loga", "logb", "logshape", "logscale"))
  # every parameter defaults to ~ 1 (RHS is the bare intercept)
  expect_true(all(vapply(f$formula$pforms, \(p) deparse1(p[[3]]) == "1", logical(1))))
})


test_that("dpar_formulas put covariates on a parameter, and extras warn", {
  f <- build_hazard_spec("exp",
         dpar_formulas = list(loga = ~ 1 + acclim + (1 | trial)))
  expect_equal(deparse1(f$formula$pforms$loga[[3]]), "1 + acclim + (1 | trial)")
  expect_equal(deparse1(f$formula$pforms$logb[[3]]), "1")   # untouched -> intercept only

  # a formula for a parameter the chosen functions do not use is dropped with a warning
  expect_warning(build_hazard_spec("exp", dpar_formulas = list(Tc = ~ 1 + acclim)),
                 class = "rTFT_unused_dpar_formulas")
})


test_that("unknown m_form, baseline, or observation error with their condition classes", {
  expect_error(build_hazard_spec("nope"), class = "rTFT_unknown_m_form")
  expect_error(build_hazard_spec("exp", baseline = "nope"), class = "rTFT_unknown_baseline")
  expect_error(build_hazard_spec("exp", observation = "nope"), class = "rTFT_bad_observation")
})


test_that("build_hazard_spec returns the formula, family, and observation", {
  f <- build_hazard_spec("exp")
  expect_named(f, c("formula", "family", "observation"))
  expect_s3_class(f$formula, "brmsformula")
  expect_s3_class(f$family, "brmsfamily")
  expect_true(isTRUE(attr(f$formula$formula, "nl")))   # a nonlinear formula
})


# ===========================================================================
# 4. Prior builder (build_hazard_prior) and its formula helpers
# ===========================================================================

# build_hazard_prior() returns a brmsprior without fitting, so the checks read the
# rows off it. The expected prior strings are written out by hand, so the test
# checks the defaults and the override wiring rather than re-running the builder.
# The tree helpers (.has_re / .has_slopes) that decide which slope and
# standard-deviation priors are added are tested directly, including a grouping-
# helper term with inner parentheses.

# The prior string for a parameter's intercept, covariate slope, or random-effect sd
# (an intercept has coef "Intercept"; a slope and an sd have the empty coef).
icpt_of  <- function(p, dp) p$prior[p$class == "b"  & p$coef == "Intercept" & p$nlpar == dp]
# The least a table needs for the default priors: its reference time and
# temperature, and a failure at the temperature failures are typical of.
prior_long <- function(tau_ref = 60, Tr = 28, T_fail = 28)
  data.frame(Ta = T_fail, failed = 1L, Tr = Tr, tau_ref = tau_ref)
slope_of <- function(p, dp) p$prior[p$class == "b"  & p$coef == ""          & p$nlpar == dp]
sd_of    <- function(p, dp) p$prior[p$class == "sd" & p$nlpar == dp]


test_that("build_hazard_prior returns a brmsprior with an intercept per needed parameter", {
  p <- build_hazard_prior("exp", "mem", long = prior_long())
  expect_s3_class(p, "brmsprior")
  # exp + mem needs loga and logb, and nothing else
  ib <- p[p$class == "b" & p$coef == "Intercept", ]
  expect_setequal(ib$nlpar, c("loga", "logb"))
})


test_that("the default intercept priors match the documented values", {
  p1 <- build_hazard_prior("power", "mem", long = prior_long())       # loga, logb, logc
  expect_equal(icpt_of(p1, "loga"), "normal(-5.094, 2)")      # normal(-1 - log(60), 2)
  expect_equal(icpt_of(p1, "logb"), "normal(-0.7, 1)")
  expect_equal(icpt_of(p1, "logc"), "normal(-8.7, 2)")        # normal(log(0.01) - log(60), 2)

  p2 <- build_hazard_prior("exp_threshold", "llogis", long = prior_long())   # + Tc, shape, scale
  expect_equal(icpt_of(p2, "Tc"), "normal(28, 3)")            # at the failures
  expect_equal(icpt_of(p2, "logshape"), "normal(0, 1)")
  expect_equal(icpt_of(p2, "logscale"), "normal(4.094, 1)")   # normal(log(60), 1)
  # the Weibull scale is fixed, so it has no prior
  expect_false("logscale" %in% build_hazard_prior("exp", "weibull", long = prior_long())$nlpar)

  p3 <- build_hazard_prior("const", "gompertz", long = prior_long())  # loga, gam
  expect_equal(icpt_of(p3, "gam"), "normal(0, 0.05)")         # normal(0, 3 / 60)
})


centre <- function(p, par) as.numeric(sub("normal\\(([^,]+),.*", "\\1", icpt_of(p, par)))
spread <- function(p, par) as.numeric(sub(".*, ([^)]+)\\)", "\\1", icpt_of(p, par)))

test_that("the priors that depend on time are scaled by the reference time", {
  # an hour-long exposure timed in minutes, in hours, and a ten-hour one in minutes
  for (tau_ref in c(60, 1, 600)) {
    mem <- build_hazard_prior("power", "mem", long = prior_long(tau_ref))
    expect_equal(centre(mem, "loga"), -1 - log(tau_ref), tolerance = 1e-3)
    expect_equal(centre(mem, "logc"), log(0.01) - log(tau_ref), tolerance = 1e-3)
    gom <- build_hazard_prior("exp", "gompertz", long = prior_long(tau_ref))
    expect_equal(centre(gom, "loga"), -1 - log(tau_ref), tolerance = 1e-3)
    expect_equal(spread(gom, "gam"), 3 / tau_ref, tolerance = 1e-3)
    llo <- build_hazard_prior("exp", "llogis", long = prior_long(tau_ref))
    expect_equal(centre(llo, "logscale"), log(tau_ref), tolerance = 1e-3)
    # with the Weibull and log-logistic, a has no time unit and its prior does not move
    expect_equal(icpt_of(llo, "loga"), "normal(-1, 2)")
    expect_equal(icpt_of(build_hazard_prior("exp", "weibull", long = prior_long(tau_ref)),
                         "loga"), "normal(-1, 2)")
    expect_equal(icpt_of(build_hazard_prior("power", "weibull", long = prior_long(tau_ref)),
                         "logc"), "normal(-4.605, 2)")           # normal(log(0.01), 2)
  }
  # what does not depend on time is left alone
  expect_equal(icpt_of(build_hazard_prior("exp", "mem", long = prior_long(600)), "logb"),
               "normal(-0.7, 1)")
})


test_that(".loga_prior carries the prior from the failures back to the reference", {
  b_median <- exp(-0.7)
  b_sd     <- sqrt((exp(1) - 1) * exp(2 * -0.7 + 1))
  # at the failures it is the prior stated there
  expect_equal(.loga_prior("exp", -1, 0), list(centre = -1, sd = 2))
  # the exponential rate: log M at the failures is loga + b d
  far <- .loga_prior("exp", -1, 14.3)
  expect_equal(far$centre, -1 - b_median * 14.3)
  expect_equal(far$sd, sqrt(4 + (14.3 * b_sd)^2))
  # a reference above the failures moves the centre the other way, as wide
  above <- .loga_prior("exp", -1, -5)
  expect_equal(above$centre, -1 + b_median * 5)
  expect_equal(above$sd, sqrt(4 + (5 * b_sd)^2))
  # the zero-type exponential: loga + log(exp(b d) - 1), close to exp far away
  zero <- .loga_prior("zeroexp", -1, 14.3)
  expect_equal(zero$centre, -1 - log(expm1(b_median * 14.3)))
  expect_equal(zero$centre, far$centre, tolerance = 1e-3)
  expect_gt(zero$sd, far$sd)
  # the power-type rates: loga + b log(d), a much gentler dependence
  pw <- .loga_prior("zeropower", -1, 14.3)
  expect_equal(pw$centre, -1 - b_median * log(14.3))
  expect_equal(pw$sd, sqrt(4 + (log(14.3) * b_sd)^2))
  expect_equal(.loga_prior("power", -1, 14.3), pw)
  # the rates that do not use the reference are left as stated
  for (form in c("const", "exp_threshold", "zeroexp_threshold", "zeropower_threshold"))
    expect_equal(.loga_prior(form, -1, 14.3), list(centre = -1, sd = 2))
})


test_that("the priors on loga and Tc follow where the failures are", {
  # failures at 28 degrees, the reference 14 degrees below them
  far <- build_hazard_prior("exp", "weibull", long = prior_long(Tr = 14, T_fail = 28))
  expect_equal(centre(far, "loga"), -1 - exp(-0.7) * 14, tolerance = 1e-3)
  expect_gt(spread(far, "loga"), 15)
  # at the failures, the usual prior
  near <- build_hazard_prior("exp", "weibull", long = prior_long(Tr = 28, T_fail = 28))
  expect_equal(icpt_of(near, "loga"), "normal(-1, 2)")
  # the critical temperature is centred at the failures, whatever the organism
  for (T_fail in c(22, 28, 40))
    expect_equal(icpt_of(build_hazard_prior("exp_threshold",
                                            long = prior_long(Tr = 14, T_fail = T_fail)), "Tc"),
                 sprintf("normal(%d, 3)", T_fail))
  # and loga of a threshold function does not depend on the reference
  expect_equal(icpt_of(build_hazard_prior("exp_threshold", "weibull",
                                          long = prior_long(Tr = 14, T_fail = 28)), "loga"),
               "normal(-1, 2)")
})


test_that("build_hazard_prior needs the table the prior is for", {
  expect_error(build_hazard_prior("exp"), class = "rTFT_no_prior_data")
  expect_error(build_hazard_prior("exp", long = 60), class = "rTFT_no_prior_data")
  expect_error(build_hazard_prior("exp", long = prior_long()[c("Ta", "failed", "Tr")]),
               class = "rTFT_no_prior_data")
  # the table from build_hazard_data() is all it needs
  ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
  expect_no_error(build_hazard_prior(
    "exp", long = build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)))
})


test_that("fit_hazard scales the default priors by the reference time of the data", {
  ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 50))
  long <- build_hazard_data(ind, ramp_traj, "interval", tau_ref = 600, checks = ramp_checks, grid_dT = Inf)
  args <- .hazard_brm_args(long, "exp", "gompertz")
  expect_equal(icpt_of(args$prior, "gam"), "normal(0, 0.005)")
  expect_equal(icpt_of(args$prior, "loga"), "normal(-7.397, 2)")    # -1 - log(600)
  # a prior given by the user is still used as given
  args <- .hazard_brm_args(long, "exp", "gompertz",
                           intercept_priors = list(gam = "normal(0, 0.2)"))
  expect_equal(icpt_of(args$prior, "gam"), "normal(0, 0.2)")
})


test_that("slope and sd priors are added only when the parameter has those terms", {
  p <- build_hazard_prior("exp", "mem",
         dpar_formulas = list(loga = ~ 1 + acclim, logb = ~ 1 + (1 | trial)), long = prior_long())
  # loga has a covariate: a slope prior, no sd prior
  expect_equal(slope_of(p, "loga"), "normal(0, 1)")
  expect_length(sd_of(p, "loga"), 0)
  # logb has a random intercept: an sd prior, no slope prior
  expect_equal(sd_of(p, "logb"), "exponential(2)")
  expect_length(slope_of(p, "logb"), 0)
  # with neither, there are no slope or sd rows at all (only intercepts)
  p0 <- build_hazard_prior("exp", "mem", long = prior_long())
  expect_equal(nrow(p0[p0$coef == "" & p0$class %in% c("b", "sd"), ]), 0L)
})


test_that("intercept_priors override the named centres and leave the rest default", {
  p <- build_hazard_prior("exp_threshold", "llogis",
         intercept_priors = list(Tc = "normal(26, 4)", loga = "normal(-6, 1.5)"), long = prior_long())
  expect_equal(icpt_of(p, "Tc"),   "normal(26, 4)")       # overridden
  expect_equal(icpt_of(p, "loga"), "normal(-6, 1.5)")     # overridden
  expect_equal(icpt_of(p, "logb"), "normal(-0.7, 1)")     # untouched default
  expect_equal(icpt_of(p, "logscale"), "normal(4.094, 1)")# untouched default
})


test_that("slope_prior and sd_prior override the generic slope and sd priors", {
  p <- build_hazard_prior("exp", "mem",
         dpar_formulas = list(loga = ~ 1 + acclim, logb = ~ 1 + (1 | trial)),
         slope_prior = "normal(0, 0.5)", sd_prior = "exponential(3)", long = prior_long())
  expect_equal(slope_of(p, "loga"), "normal(0, 0.5)")
  expect_equal(sd_of(p, "logb"), "exponential(3)")
})


test_that("intercept_priors is validated: unnamed errors, an unused name warns", {
  expect_error(build_hazard_prior("exp", intercept_priors = list("normal(0, 1)"), long = prior_long()),
               class = "rTFT_bad_intercept_priors")
  # Tc is not a parameter of the plain exponential form
  expect_warning(build_hazard_prior("exp", intercept_priors = list(Tc = "normal(1, 1)"), long = prior_long()),
                 class = "rTFT_unused_intercept_prior")
})


test_that("build_hazard_prior rejects an unknown m_form or baseline", {
  expect_error(build_hazard_prior("nope", long = prior_long()), class = "rTFT_unknown_m_form")
  expect_error(build_hazard_prior("exp", baseline = "nope", long = prior_long()), class = "rTFT_unknown_baseline")
})


test_that(".has_re and .has_slopes read the formula parse tree, not its text", {
  expect_true(.has_slopes(~ 1 + a + b))            # several predictors
  expect_false(.has_slopes(~ 1))                   # intercept only
  expect_false(.has_slopes(~ 1 + (1 | g)))         # random effect only, no fixed slope
  expect_true(.has_slopes(~ 1 + a * b))            # an interaction is a slope
  expect_true(.has_re(~ 1 + (1 | g) + (1 | h)))    # several random effects
  expect_true(.has_re(~ 1 + (1 + a || g)))           # uncorrelated effects, a double bar
  expect_false(.has_re(~ 1 + a))
  expect_false(.has_re(NULL))
  expect_false(.has_slopes(NULL))
})


test_that("a grouping helper term with inner parentheses is still one random effect", {
  # (1 | gr(site, cor = TRUE)) has an inner ")" that stopping at the first ")" would trip on.
  f <- ~ 1 + acclim + (1 | gr(site, cor = TRUE))
  expect_true(.has_re(f))                          # the group is recognised
  expect_true(.has_slopes(f))                      # acclim is still a fixed slope
  # and it flows through to the prior: both an sd prior and the slope prior appear
  p <- build_hazard_prior("exp", dpar_formulas = list(loga = f), long = prior_long())
  expect_equal(sd_of(p, "loga"), "exponential(2)")
  expect_equal(slope_of(p, "loga"), "normal(0, 1)")
})


# ===========================================================================
# 5. Model fitting (fit_hazard: the assembled brm() arguments and the Tr guard)
# ===========================================================================

# fit_hazard() itself calls brms::brm(), which compiles a model and cannot run under
# R CMD check, so the checks target .hazard_brm_args(), which includes everything
# fit_hazard does except the fit, and the .check_Tr() guard. A minimal interval table
# stands in for real data; nothing here fits a model.

fit_ind  <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
long_int <- build_hazard_data(fit_ind, ramp_traj, "interval", checks = ramp_checks, grid_dT = Inf)


test_that(".hazard_brm_args assembles the formula, data, family, prior and stanvars", {
  a <- .hazard_brm_args(long_int, "exp", observation = "interval")
  expect_setequal(names(a), c("formula", "data", "family", "prior", "stanvars"))
  expect_s3_class(a$formula, "brmsformula")
  expect_identical(a$data, long_int)                     # the data is passed through unchanged
  expect_equal(a$family$name, "cloglog_bern")
  expect_s3_class(a$prior, "brmsprior")
  expect_identical(a$stanvars, .cloglog_bern_stanvars)   # the interval family's Stan code
})


test_that("the family and stanvars follow the observation mode", {
  # continuous -> built-in Poisson, no stanvars
  long_c <- build_hazard_data(data.frame(id = 1, path_id = "A", event = 1, t_event = 42),
                              ramp_traj, "continuous", checks = ramp_checks, grid_dT = Inf)
  ac <- .hazard_brm_args(long_c, "exp", observation = "continuous")
  expect_equal(ac$family$family, "poisson")
  expect_null(ac$stanvars)

  # mixed -> the hazard_mixed family and its Stan code
  long_m <- build_hazard_data(
    data.frame(id = 1:2, path_id = "A", event = 1, obs_type = c("interval", "continuous"),
               t_upper = c(30, NA), t_event = c(NA, 42)),
    ramp_traj, "mixed", checks = ramp_checks, grid_dT = Inf)
  am <- .hazard_brm_args(long_m, "exp", observation = "mixed")
  expect_equal(am$family$name, "hazard_mixed")
  expect_identical(am$stanvars, .hazard_mixed_stanvars)
})


test_that("missing values in a variable the model uses are refused", {
  with_na <- function(column, rows = 2) {
    out <- long_int
    if (!column %in% names(out)) out[[column]] <- seq_len(nrow(out))
    out[[column]][rows] <- NA
    out
  }
  # the temperature and the exposure are used by every model
  expect_error(.hazard_brm_args(with_na("Ta"), "exp"), class = "rTFT_missing_values")
  expect_error(.hazard_brm_args(with_na("log_exposure"), "exp"), class = "rTFT_missing_values")
  # elapsed time only by a baseline that depends on it
  expect_no_error(.hazard_brm_args(with_na("log_tau"), "exp", "mem"))
  expect_error(.hazard_brm_args(with_na("log_tau"), "exp", "weibull"),
               class = "rTFT_missing_values")
  # a covariate or a grouping only where a parameter's formula has it
  expect_no_error(.hazard_brm_args(with_na("mass"), "exp"))
  expect_error(.hazard_brm_args(with_na("mass"), "exp", dpar_formulas = list(loga = ~ 1 + mass)),
               class = "rTFT_missing_values")
  expect_error(.hazard_brm_args(with_na("trial"), "exp",
                                dpar_formulas = list(loga = ~ 1 + (1 | trial))),
               class = "rTFT_missing_values")
  # the message names the columns and counts the rows, each row once
  both <- with_na("Ta", rows = 1:2); both$log_exposure[2:3] <- NA
  expect_error(.hazard_brm_args(both, "exp"), "3 of [0-9]+ rows", class = "rTFT_missing_values")
  expect_error(.hazard_brm_args(both, "exp"), "Ta.*log_exposure", class = "rTFT_missing_values")
})


test_that("the default prior is built when none is supplied, and a supplied prior is used as-is", {
  a <- .hazard_brm_args(long_int, "exp")                 # default prior built
  expect_equal(a$prior$prior[a$prior$nlpar == "loga" & a$prior$coef == "Intercept"],
               paste0("normal(", signif(-1 - log(long_int$tau_ref[1]), 4), ", 2)"))

  mine <- brms::set_prior("normal(-3, 1)", class = "b", coef = "Intercept", nlpar = "loga")
  a2 <- .hazard_brm_args(long_int, "exp", prior = mine)  # supplied prior passed straight through
  expect_identical(a2$prior, mine)
})


test_that("prior-override arguments are forwarded to build_hazard_prior", {
  a <- .hazard_brm_args(long_int, "exp",
         dpar_formulas = list(loga = ~ 1 + acclim),
         intercept_priors = list(loga = "normal(-6, 1)"),
         slope_prior = "normal(0, 0.5)")
  expect_equal(a$prior$prior[a$prior$nlpar == "loga" & a$prior$coef == "Intercept"],
               "normal(-6, 1)")
  expect_equal(a$prior$prior[a$prior$nlpar == "loga" & a$prior$coef == ""],
               "normal(0, 0.5)")
})


test_that(".check_Tr stops the forms that need Ta above Tr, and passes the others", {
  # Tr set above the data, so the log / fractional-power forms are undefined
  hot <- build_hazard_data(fit_ind, ramp_traj, "interval", checks = ramp_checks, Tr = 25, grid_dT = Inf)
  for (m in c("zeroexp", "zeropower", "power"))
    expect_error(.check_Tr(hot, m), class = "rTFT_Tr_below_data")
  # forms without a logarithm or fractional power of (Ta - Tr) are unaffected
  for (m in c("exp", "const", "exp_threshold"))
    expect_silent(.check_Tr(hot, m))
})


test_that(".check_Tr is strict for zeroexp/zeropower but allows power at Ta == Tr", {
  # a flat cold start so the coldest interval midpoint equals the default Tr
  flat <- data.frame(path_id = "A", time = c(0, 25, 60), temp = c(20, 20, 35))
  eq   <- build_hazard_data(fit_ind, flat, "interval", checks = ramp_checks, grid_dT = Inf)  # min(Ta) == Tr == 20
  expect_error(.check_Tr(eq, "zeroexp"),   class = "rTFT_Tr_below_data")       # -Inf at equality
  expect_error(.check_Tr(eq, "zeropower"), class = "rTFT_Tr_below_data")
  expect_silent(.check_Tr(eq, "power"))                                        # 0^b = 0 is fine
})


# ===========================================================================
# 6. Fixed-ramp helper (ramp_path)
# ===========================================================================

# ramp_path() builds a linear ramp trajectory to stand in for the `path` argument.
# Expected temperatures are written out from temp(t) = T_start + rate * t.

test_that("ramp_path builds a two-point linear path by default", {
  p <- ramp_path("A", T_start = 20, rate = 0.25, t_max = 60)
  expect_equal(names(p), c("path_id", "time", "temp"))
  expect_equal(p$time, c(0, 60))                     # just the endpoints
  expect_equal(p$temp, c(20, 35))                    # 20 + 0.25 * c(0, 60)
  expect_true(all(p$path_id == "A"))
})


test_that("a `by` step samples the ramp and always closes on t_max", {
  p <- ramp_path("A", T_start = 20, rate = 0.2, t_max = 75, by = 10)
  # seq(0, 75, 10) = 0..70, plus the exact 75 (not a multiple of the step)
  expect_equal(p$time, c(seq(0, 70, 10), 75))
  expect_equal(p$temp, 20 + 0.2 * p$time)            # linear throughout
})


test_that("ramp_path is vectorised over several ramps and recycles a length-one arg", {
  p <- ramp_path(c("slow", "fast"), T_start = 20, rate = c(0.2, 0.4),
                 t_max = c(75, 40), by = 10)
  expect_setequal(unique(p$path_id), c("slow", "fast"))
  expect_equal(nrow(p), 14L)                         # slow: 0..70,75 (9); fast: 0..40 (5)
  # T_start = 20 was recycled to both ramps
  expect_equal(p$temp[p$path_id == "slow"][1], 20)
  expect_equal(p$temp[p$path_id == "fast"][1], 20)
  # each ramp carries its own realised slope
  expect_equal(coef(lm(temp ~ time, p[p$path_id == "slow", ]))[["time"]], 0.2)
})


test_that("a ramp_path output feeds build_hazard_data and reproduces the ramp temperature", {
  p    <- ramp_path("A", T_start = 20, rate = 0.25, t_max = 60, by = 10)
  ind  <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  long <- build_hazard_data(ind, p, "interval", grid_dT = Inf)      # p's `by` grid is the check schedule
  # the failure interval [20, 30] has midpoint 25; the ramp gives 20 + 0.25 * 25 = 26.25
  expect_equal(long$Ta[long$tau == 25], 26.25)
})


# ===========================================================================
# 7. assess_ramp_linearity
# ===========================================================================

# assess_ramp_linearity() regresses temp ~ time per path and reports the realised
# ramp and its deviations. Expected values come from the line temp = T_start +
# rate * time written out by hand; the deviation summaries are checked by their
# qualitative signature (a cold-end vs a hot-end departure), which needs no fitted
# residuals.

test_that("a perfectly linear trace recovers the ramp and shows no deviation", {
  p <- data.frame(path_id = "r1", time = 0:20, temp = 20 + 0.3 * (0:20))
  a <- assess_ramp_linearity(p)
  expect_equal(nrow(a), 1L)
  expect_equal(a$T_start, 20)
  expect_equal(a$rate, 0.3)
  expect_equal(a$r_squared, 1)
  expect_true(a$rmse < 1e-8)
  expect_true(a$max_abs_dev < 1e-8)
})


test_that("nominal adds the realised-minus-nominal rate bias", {
  p <- data.frame(path_id = "r1", time = 0:20, temp = 20 + 0.3 * (0:20))
  a <- assess_ramp_linearity(p, nominal = c(r1 = 0.25))
  expect_equal(a$rate_nominal, 0.25)
  expect_equal(a$rate_bias, 0.05)              # 0.30 realised - 0.25 nominal
})


test_that("max_dev_hot isolates the hot quartile from a cold-end departure", {
  base <- 20 + 0.3 * (0:20)
  cold <- base; cold[3] <- cold[3] + 2         # a bump at t = 2, in the cold quarter
  a <- assess_ramp_linearity(data.frame(path_id = "r1", time = 0:20, temp = cold))
  expect_lt(a$r_squared, 1)                    # no longer a perfect line
  expect_lt(a$max_dev_hot, a$max_abs_dev)      # the worst deviation is cold, not hot

  hot <- base; hot[21] <- hot[21] + 2          # a bump at t = 20, the hottest point
  b <- assess_ramp_linearity(data.frame(path_id = "r1", time = 0:20, temp = hot))
  expect_equal(b$max_dev_hot, b$max_abs_dev)   # the worst deviation is in the hot quartile
})


test_that("one row per path, and nominal fills only the paths it names", {
  p <- rbind(
    data.frame(path_id = "a", time = 0:10, temp = 20 + 0.3 * (0:10)),
    data.frame(path_id = "b", time = 0:10, temp = 22 + 0.5 * (0:10))
  )
  a <- assess_ramp_linearity(p, nominal = c(a = 0.3))   # only "a" has a nominal rate
  expect_setequal(a$path_id, c("a", "b"))
  expect_equal(a$rate_bias[a$path_id == "a"], 0)        # 0.3 realised - 0.3 nominal
  expect_true(is.na(a$rate_bias[a$path_id == "b"]))     # no nominal for "b"
})


test_that("assess_ramp_linearity errors on a missing column", {
  expect_error(assess_ramp_linearity(data.frame(path_id = "a", time = 0:5)),
               class = "rTFT_missing_column")
})


# ===========================================================================
# 8. Co-stressor layer (costressor and its wiring)
# ===========================================================================

# costressor() declares a non-temperature stressor once; the same object flows
# into build_hazard_data (which adds the standardised column z_<name>),
# build_hazard_spec (which appends c<name> and dT<name> to log M), and
# build_hazard_prior (which sets their priors). Dissolved oxygen is the worked
# case: z_DO = -(do - 8.5), so a fully oxygenated trial (do = 8.5) sits at z = 0
# and reduces to the pure M(T). Expected z values and predictor strings are written
# out by hand, not read back from the function.

do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                    interaction = TRUE)


test_that("costressor stores the declaration and carries its class", {
  expect_s3_class(do_cs, "rTFT_costressor")
  expect_equal(do_cs$name, "DO")
  expect_equal(do_cs$column, "do")
  expect_equal(do_cs$ref, 8.5)
  expect_equal(do_cs$scale, 1)
  expect_equal(do_cs$direction, -1)
  expect_true(do_cs$interaction)            # the temperature interaction, asked for
  # without asking, the co-stressor enters through its main effect alone
  expect_false(costressor("DO", "do", 8.5, 1, -1)$interaction)
  expect_null(do_cs$trajectory)
  expect_null(do_cs$transform)
})


test_that("costressor validates each field with its own condition class", {
  ok <- list(column = "do", ref = 8.5, scale = 1, direction = -1)
  expect_error(costressor("D_O", "do", 8.5, 1, -1), class = "rTFT_bad_costressor_name")
  expect_error(costressor("2do", "do", 8.5, 1, -1), class = "rTFT_bad_costressor_name")
  expect_error(costressor(c("a", "b"), "do", 8.5, 1, -1), class = "rTFT_bad_costressor_name")
  expect_error(costressor("DO", "", 8.5, 1, -1), class = "rTFT_bad_costressor_column")
  expect_error(costressor("DO", "do", NA_real_, 1, -1), class = "rTFT_bad_costressor_ref")
  expect_error(costressor("DO", "do", 8.5, 0, -1), class = "rTFT_bad_costressor_scale")
  expect_error(costressor("DO", "do", 8.5, -1, -1), class = "rTFT_bad_costressor_scale")
  expect_error(costressor("DO", "do", 8.5, 1, 2), class = "rTFT_bad_costressor_direction")
  expect_error(costressor("DO", "do", 8.5, 1, -1, interaction = NA), class = "rTFT_bad_costressor_interaction")
  expect_error(costressor("DO", "do", 8.5, 1, -1, trajectory = 1:3), class = "rTFT_bad_costressor_trajectory")
  expect_error(costressor("DO", "do", 8.5, 1, -1, transform = "sq"), class = "rTFT_bad_costressor_transform")
})


test_that(".check_costressors rejects a non-costressor list and duplicate names", {
  expect_error(.check_costressors(list(list(name = "DO"))), class = "rTFT_bad_costressors")
  do2 <- costressor("DO", "do2", 8.5, 1, -1)
  expect_error(.check_costressors(list(do_cs, do2)), class = "rTFT_duplicate_costressor")
  expect_silent(.check_costressors(list()))          # none is fine
})


test_that(".costressor_dpars and .costressor_logM emit the parameters and terms", {
  # main effect plus interaction by default
  expect_equal(.costressor_dpars(list(do_cs)), c("cDO", "dTDO"))
  expect_equal(.costressor_logM(list(do_cs)),
               "cDO * z_DO + dTDO * (Ta - Tr) * z_DO")
  # interaction off drops the dT term and its parameter
  do_noint <- costressor("DO", "do", 8.5, 1, -1, interaction = FALSE)
  expect_equal(.costressor_dpars(list(do_noint)), "cDO")
  expect_equal(.costressor_logM(list(do_noint)), "cDO * z_DO")
  # two stressors accumulate
  sal <- costressor("Sal", "sal", 0, 5, 1, interaction = TRUE)
  expect_equal(.costressor_dpars(list(do_cs, sal)), c("cDO", "dTDO", "cSal", "dTSal"))
  # empty
  expect_equal(.costressor_dpars(list()), character(0))
  expect_equal(.costressor_logM(list()), "")
})


test_that("build_hazard_data adds z_<name> for a fixed-per-trial co-stressor", {
  ind  <- data.frame(id = 1:2, path_id = "A", event = c(1, 0),
                     t_upper = c(30, 60), do = c(6, 8.5))
  out  <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                            costressors = list(do_cs), grid_dT = Inf)
  expect_true("z_DO" %in% names(out))
  # z = -(do - 8.5): individual 1 (do = 6) -> 2.5 on every row; individual 2 -> 0
  expect_true(all(out$z_DO[out$id == 1] == 2.5))
  expect_true(all(out$z_DO[out$id == 2] == 0))       # a fully oxygenated trial sits at z = 0
})


test_that("build_hazard_data interpolates a time-varying co-stressor at the midpoints", {
  # dissolved oxygen falls linearly from 9 to 5 over the hour: do(t) = 9 - (4/60) t
  do_path <- data.frame(path_id = "A", time = c(0, 60), do = c(9, 5))
  ind     <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  cs_tv   <- costressor("DO", "do", ref = 8.5, scale = 1, direction = -1,
                        trajectory = do_path)
  out     <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                               costressors = list(cs_tv), grid_dT = Inf)
  # z at each midpoint, computed a second way with approx()
  do_mid  <- stats::approx(do_path$time, do_path$do, xout = out$tau)$y
  expect_equal(out$z_DO, -(do_mid - 8.5))
})


test_that("build_hazard_data applies a supplied transform in place of the signed-linear map", {
  ind  <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30, sal = 10)
  # a U-shaped salinity response: the squared standardised deviation
  cs_u <- costressor("Sal", "sal", ref = 5, scale = 5, direction = 1,
                     transform = function(u) u^2)
  out  <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                            costressors = list(cs_u), grid_dT = Inf)
  expect_true(all(out$z_Sal == ((10 - 5) / 5)^2))    # ((X - ref)/scale)^2 = 1
})


test_that("build_hazard_data errors on a missing co-stressor column or trajectory path", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)   # no `do` column
  expect_error(build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                                 costressors = list(do_cs), grid_dT = Inf),
               class = "rTFT_missing_column")
  # a trajectory that names a path the data never visits
  do_path <- data.frame(path_id = "B", time = c(0, 60), do = c(9, 5))
  cs_tv   <- costressor("DO", "do", 8.5, 1, -1, trajectory = do_path)
  expect_error(build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                                 costressors = list(cs_tv), grid_dT = Inf),
               class = "rTFT_missing_costressor_path")
})


test_that("build_hazard_spec appends the co-stressor terms to the predictor", {
  eta <- "loga + exp(logb) * (Ta - Tr) + cDO * z_DO + dTDO * (Ta - Tr) * z_DO + log_exposure"
  expect_equal(main_rhs(build_hazard_spec("exp", costressors = list(do_cs))),
               str2lang(eta))
  # interaction off leaves only the main effect
  do_noint <- costressor("DO", "do", 8.5, 1, -1, interaction = FALSE)
  eta0 <- "loga + exp(logb) * (Ta - Tr) + cDO * z_DO + log_exposure"
  expect_equal(main_rhs(build_hazard_spec("exp", costressors = list(do_noint))),
               str2lang(eta0))
})


test_that("co-stressor coefficients get a sub-formula, and dpar_formulas can target one", {
  f <- build_hazard_spec("exp", costressors = list(do_cs))
  expect_setequal(names(f$formula$pforms), c("loga", "logb", "cDO", "dTDO"))
  # a covariate on the oxygen main effect, carried through without a warning
  f2 <- build_hazard_spec("exp", costressors = list(do_cs),
                          dpar_formulas = list(cDO = ~ 1 + acclim))
  expect_equal(deparse1(f2$formula$pforms$cDO[[3]]), "1 + acclim")
})


test_that("build_hazard_prior sets the co-stressor coefficient priors", {
  p <- build_hazard_prior("exp", costressors = list(do_cs), long = prior_long())
  expect_equal(icpt_of(p, "cDO"), "normal(0, 1)")      # main effect, wider
  expect_equal(icpt_of(p, "dTDO"), "normal(0, 0.5)")   # interaction, tighter
  # interaction off: no dTDO prior at all
  do_noint <- costressor("DO", "do", 8.5, 1, -1, interaction = FALSE)
  p0 <- build_hazard_prior("exp", costressors = list(do_noint), long = prior_long())
  expect_length(icpt_of(p0, "dTDO"), 0)
  # and intercept_priors overrides a co-stressor coefficient
  p2 <- build_hazard_prior("exp", costressors = list(do_cs),
                           intercept_priors = list(cDO = "normal(0, 2)"), long = prior_long())
  expect_equal(icpt_of(p2, "cDO"), "normal(0, 2)")
})


test_that(".hazard_brm_args threads a co-stressor into the formula and the default prior", {
  ind  <- data.frame(id = 1:2, path_id = "A", event = c(1, 0),
                     t_upper = c(30, 60), do = c(6, 8.5))
  long <- build_hazard_data(ind, ramp_traj, "interval", checks = ramp_checks,
                            costressors = list(do_cs), grid_dT = Inf)
  a    <- .hazard_brm_args(long, "exp", costressors = list(do_cs))
  expect_true(grepl("cDO * z_DO", deparse1(main_rhs(a)), fixed = TRUE))
  expect_equal(a$prior$prior[a$prior$nlpar == "cDO" & a$prior$coef == "Intercept"],
               "normal(0, 1)")
})


test_that("fit_hazard passes a custom k_sw into the formula", {
  # A steepness set on fit_hazard() must reach the model, not be forwarded to
  # brms::brm() through `...`.
  expect_true("k_sw" %in% names(formals(fit_hazard)))
  a <- .hazard_brm_args(long_int, "exp_threshold", k_sw = 20)
  expect_true(grepl("-20 * (Ta - Tc)", deparse1(main_rhs(a)), fixed = TRUE))
})


test_that(".hazard_settings records what the reporting functions need", {
  st <- .hazard_settings(long_int, "exp_threshold", "mem", "interval", k_sw = NULL,
                         costressors = list())
  expect_equal(st[c("m_form", "baseline", "observation")],
               list(m_form = "exp_threshold", baseline = "mem", observation = "interval"))
  expect_equal(st$Tr, long_int$Tr[1])                    # the reference of the data
  expect_equal(st$tau_ref, long_int$tau_ref[1])          # and the reference time
  expect_equal(st$Ta_range, range(long_int$Ta))          # the ranges the data cover
  expect_equal(st$tau_range, range(long_int$tau))
  expect_equal(st$id, long_int$id)                       # the individual of each row
  expect_equal(st$k_sw, 50)                              # the catalogue steepness
  expect_equal(.hazard_settings(long_int, "exp_threshold", "mem", "interval", 20,
                                list())$k_sw, 20)        # a custom one
  # a k_sw ignored for a threshold-less function is recorded as the catalogue 0,
  # without a second warning
  expect_silent(st2 <- .hazard_settings(long_int, "exp", "mem", "interval", 20, list()))
  expect_equal(st2$k_sw, 0)
})


test_that(".fit_settings refuses a fit that carries no settings", {
  expect_error(.fit_settings(list()), class = "rTFT_not_brmsfit")
  bare <- structure(list(), class = "brmsfit")
  expect_error(.fit_settings(bare), class = "rTFT_no_settings")
})
