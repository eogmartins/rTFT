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
# 0b. Unified family for mixed interval-censored + known-time endpoints
# ===========================================================================

# Mixed interval-censored and known-time response ------------------------------
#
# The family declaration (hazard_mixed_family), its Stan log-mass and density
# function (.hazard_mixed_stanvars), and the three R methods brms calls after
# fitting.

# The Stan log-mass and density definition, written out in R as this family's
# reference. Three branches, keyed on the survived indicator y and the known-time
# flag kt: a survived interval contributes -Lambda, an interval-censored endpoint
# log(1 - exp(-Lambda)), and a known-time endpoint the log density mu - Lambda,
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
  # fit contain both kinds of endpoint.
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
  y  <- c(0, 1, 1)                           # survived, endpoint, endpoint
  kt <- c(0, 0, 1)                           # kt is ignored for the survived row
  prep <- stub_prep(mu, y, vint1 = kt)

  for (i in seq_len(3)) {
    ll <- log_lik_hazard_mixed(i, prep)
    expect_equal(ll, ref_loglik_mixed(mu[, i], y[i], kt[i]))
    expect_length(ll, 200)
  }
})


test_that("vint1 selects the endpoint term, and is ignored when the animal survived", {
  set.seed(5)
  mu  <- matrix(rnorm(50), nrow = 50, ncol = 1)   # 50 draws, 1 observation
  lam <- exp(mu[, 1])

  # On an endpoint row (y = 1) the flag picks interval vs known-time, and they differ.
  ll_interval <- log_lik_hazard_mixed(1L, stub_prep(mu, y = 1, vint1 = 0))
  ll_known    <- log_lik_hazard_mixed(1L, stub_prep(mu, y = 1, vint1 = 1))
  expect_equal(ll_interval, log(1 - exp(-lam)))    # interval-censored term
  expect_equal(ll_known,    mu[, 1] - lam)         # known-time term
  expect_false(isTRUE(all.equal(ll_interval, ll_known)))
  # The interval-censored term is exactly the cloglog family's endpoint term.
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


test_that("brms_log_M threshold guards stay finite below and at Tc; zero forms vanish at Tr", {
  a <- 2.82e-7; b <- 6.13; cc <- 1e-4; Tc <- 36; Tr <- 30; k_sw <- 50
  ev <- function(s, Ta) {
    env <- list(loga = log(a), logb = log(b), logc = log(cc),
                Tc = Tc, Tr = Tr, k_sw = k_sw, Ta = Ta)
    exp(eval(parse(text = s), env))
  }
  # abs(...) + the 1e-12 literal keep the below-Tc and at-Tc logarithms finite
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
  shape <- 1.8; scale <- 20; gam <- 0.05; tau <- 12
  env <- list(logshape = log(shape), logscale = log(scale), gam = gam,
              tau = tau, log_tau = log(tau))
  ev <- function(s) exp(eval(parse(text = s), env))

  z <- (tau / scale)^shape
  h_weibull  <- (shape / scale) * (tau / scale)^(shape - 1)
  h_llogis   <- (shape / scale) * (tau / scale)^(shape - 1) / (1 + z)   # hazard
  h_gompertz <- exp(gam * tau)                                          # amplitude in loga
  f_llogis   <- h_llogis / (1 + z)                                      # density (factor-2 error)

  expect_equal(ev(brms_log_h0[["weibull"]]),  h_weibull)
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

ramp_path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
ramp_temp   <- function(t) 20 + 0.25 * t                # the ramp
ramp_checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))


test_that("interval mode gives one row per at-risk interval up to the failure check", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  out <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks)

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
  out <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks)
  fail_row <- out[out$failed == 1, ]
  expect_equal(nrow(fail_row), 1L)
  expect_equal(fail_row$tau, 25)                        # midpoint of [20,30]
})


test_that("a right-censored survivor never fails and runs to its last check", {
  ind <- data.frame(id = 2, path_id = "A", event = 0, t_upper = 60)
  out <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks)
  expect_equal(nrow(out), 6L)                           # [0,10] .. [50,60]
  expect_true(all(out$failed == 0L))                    # event = 0 -> no failure row
  expect_equal(max(out$tau), 55)                        # last interval midpoint
})


test_that("the output columns and their derived quantities are correct", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  out <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks)

  expect_setequal(names(out),
                  c("id", "path_id", "tau", "exposure", "Ta", "failed",
                    "Tr", "log_tau", "log_exposure"))
  expect_equal(out$log_tau, log(out$tau))
  expect_equal(out$log_exposure, log(out$exposure))
  expect_false("exact" %in% names(out))                 # exact is mixed-mode only
})


test_that("Tr defaults to the coldest path temperature and can be overridden", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  out_def <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks)
  expect_true(all(out_def$Tr == 20))                    # min(ramp_path$temp)

  out_15 <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks, Tr = 15)
  expect_true(all(out_15$Tr == 15))                     # a constant column at the set value
})


test_that("covariates are replicated onto every row and moved to the front", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30, mass = 12.1)
  out <- build_hazard_data(ind, ramp_path, "interval", checks = ramp_checks,
                           covariates = "mass")
  expect_true(all(out$mass == 12.1))                    # the per-individual value, on each row
  expect_equal(names(out)[1:2], c("id", "mass"))        # id then covariates lead the table
})


test_that("continuous mode truncates the final interval at the exact failure time", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  out <- build_hazard_data(ind, ramp_path, "continuous", checks = ramp_checks)

  # nodes 0,10,20,30,40,42 -> five intervals, the last [40,42].
  expect_equal(nrow(out), 5L)
  expect_equal(out$exposure, c(10, 10, 10, 10, 2))      # final interval is 2 wide
  expect_equal(out$tau, c(5, 15, 25, 35, 41))           # last midpoint is (40 + 42)/2
  expect_equal(out$failed, c(0L, 0L, 0L, 0L, 1L))
})


test_that("grid_dt integrates a known-time failure on a regular grid ending at the failure time", {
  ind <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
  out <- build_hazard_data(ind, ramp_path, "continuous", grid_dt = 5)

  # seq(0, 42, 5) = 0,5,..,40, plus the exact 42 -> nine intervals, the last [40,42].
  expect_equal(nrow(out), 9L)
  expect_equal(out$exposure, c(rep(5, 8), 2))
  expect_equal(out$Ta, ramp_temp(out$tau))              # still the true midpoint temperatures
})


test_that("mixed mode sets exact per obs_type and keeps a single `failed` response", {
  ind <- data.frame(id = 1:2, path_id = "A", event = 1,
                    obs_type = c("interval", "continuous"),
                    t_upper = c(30, NA), t_event = c(NA, 42))
  out <- build_hazard_data(ind, ramp_path, "mixed", checks = ramp_checks)

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
  out  <- build_hazard_data(ind, fine, "interval", checks = ramp_checks)
  expect_equal(out$Ta[out$tau == 5], 20)                                  # on the flat leg
  expect_equal(out$Ta[out$tau == 35], 20 + (35 - 20) / (60 - 25) * (35 - 25))  # on the ramp

  # Flat extrapolation: a check beyond the last reading still gets a finite Ta.
  short <- data.frame(path_id = "A", time = c(0, 40), temp = c(20, 30))  # series ends at 40
  out2  <- build_hazard_data(ind, short, "interval", checks = ramp_checks)  # checks reach 60
  expect_true(all(is.finite(out2$Ta)))
  expect_equal(max(out2$Ta), 30)                        # held flat at the last reading
})


test_that("build_hazard_data errors loudly, with cli condition classes, on malformed input", {
  ok_ind <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)

  # an unknown observation mode
  expect_error(build_hazard_data(ok_ind, ramp_path, observation = "nope"),
               class = "rTFT_bad_observation")
  # a missing required input column (path lacks temp)
  bad_path <- data.frame(path_id = "A", time = c(0, 60))
  expect_error(build_hazard_data(ok_ind, bad_path, "interval"),
               class = "rTFT_missing_column")
  # the non-mixed timing column absent (interval needs t_upper)
  no_time <- data.frame(id = 1, path_id = "A", event = 1)
  expect_error(build_hazard_data(no_time, ramp_path, "interval"),
               class = "rTFT_missing_column")
  # a non-finite timing value
  na_time <- data.frame(id = 1, path_id = "A", event = 1, t_upper = NA_real_)
  expect_error(build_hazard_data(na_time, ramp_path, "interval"),
               class = "rTFT_nonfinite_timing")
})


test_that("mixed mode guards obs_type and its per-level timing columns", {
  # obs_type missing entirely
  m1 <- data.frame(id = 1, path_id = "A", event = 1, t_upper = 30)
  expect_error(build_hazard_data(m1, ramp_path, "mixed"),
               class = "rTFT_missing_obs_type")
  # an unrecognised obs_type level
  m2 <- data.frame(id = 1, path_id = "A", event = 1, obs_type = "typo", t_upper = 30)
  expect_error(build_hazard_data(m2, ramp_path, "mixed"),
               class = "rTFT_bad_obs_type")
  # a valid level whose timing column is absent
  m3 <- data.frame(id = 1, path_id = "A", event = 1, obs_type = "continuous")
  expect_error(build_hazard_data(m3, ramp_path, "mixed"),
               class = "rTFT_missing_timing_col")
  # a valid level whose timing value is non-finite
  m4 <- data.frame(id = 1, path_id = "A", event = 1, obs_type = "interval",
                   t_upper = NA_real_)
  expect_error(build_hazard_data(m4, ramp_path, "mixed"),
               class = "rTFT_nonfinite_timing")
})
