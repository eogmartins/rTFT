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
