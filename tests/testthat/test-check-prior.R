# ===========================================================================
# 1. Helpers
# ===========================================================================

# Two fitted models: `fit_exp` has 40 individuals on one ramp, with body mass as a
# covariate of log a, and `fit_trials` has eight trials of six individuals on one
# ramp, with a random intercept per trial. Both are interval-censored.
fit_exp    <- readRDS(test_path("fixtures", "hazard-fit-exp.rds"))
fit_trials <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))

# A table of priors as brms writes one, a row per prior, with the columns the check
# reads.
prior_table <- function(prior, coef, nlpar, class = "b", group = "", lb = "", ub = "")
  data.frame(prior = prior, class = class, coef = coef, group = group, nlpar = nlpar,
             lb = lb, ub = ub)


# ===========================================================================
# 2. The log prior and the log-likelihood of each draw
# ===========================================================================

test_that("the prior in force on a coefficient is its own, else its parameter's", {
  pt <- prior_table(c("normal(0, 1)", "normal(-7, 6)", "", "normal(0, 5)", "", ""),
                    coef = c("", "Intercept", "mass", "speciesB", "", "Intercept"),
                    nlpar = c("loga", "loga", "loga", "loga", "logb", "logb"))
  expect_equal(.coef_prior(pt, "loga", "Intercept"), "normal(-7, 6)")
  expect_equal(.coef_prior(pt, "loga", "speciesB"), "normal(0, 5)")
  expect_equal(.coef_prior(pt, "loga", "mass"), "normal(0, 1)")   # inherited
  expect_equal(.coef_prior(pt, "logb", "Intercept"), "")          # flat
  # a bound on the coefficients adds a correction the check does not evaluate
  pt$lb[1] <- "0"
  expect_error(.coef_prior(pt, "loga", "mass"), class = "rTFT_prior_unsupported")
  # as in the fitted models
  expect_equal(.coef_prior(as.data.frame(fit_exp$prior), "loga", "mass"), "normal(0, 1)")
})


test_that("the log density of a prior is that brms adds to the log prior", {
  x <- c(-3, 0, 1.5)
  expect_equal(.prior_lpdf(x, "normal(-0.7, 1)", "b"), dnorm(x, -0.7, 1, log = TRUE))
  expect_equal(.prior_lpdf(x, "normal(log(60), 2)", "b"), dnorm(x, log(60), 2, log = TRUE))
  # brms's student_t(nu, mu, sigma)
  expect_equal(.prior_lpdf(x, "student_t(3, 1, 2.5)", "b"),
               dt((x - 1) / 2.5, df = 3, log = TRUE) - log(2.5))
  expect_equal(.prior_lpdf(x, "", "b"), c(0, 0, 0))               # flat: left out
  expect_error(.prior_lpdf(x, "cauchy(0, 1)", "b"), class = "rTFT_prior_unsupported")
  expect_error(.prior_lpdf(x, "normal(0, 0)", "b"), class = "rTFT_prior_unsupported")
})


test_that("the log prior by block adds up to the log prior brms kept", {
  # no group-level terms: a block per parameter, nothing left over
  d <- as.data.frame(brms::as_draws_df(fit_exp))
  pb <- .prior_blocks(fit_exp, d)
  expect_named(pb$blocks, c("loga", "logb"))
  expect_equal(pb$coefs$variable, c("b_loga_Intercept", "b_loga_mass", "b_logb_Intercept"))
  expect_equal(pb$coefs$block, c("loga", "loga", "logb"))
  expect_equal(pb$coefs$prior[2], "normal(0, 1)")
  expect_equal(Reduce(`+`, pb$blocks), d$lprior, tolerance = 1e-6)
  expect_equal(pb$blocks$logb, dnorm(d$b_logb_Intercept, -0.7, 1, log = TRUE))
  # a trial SD: its prior is the block sd, the rest of the log prior
  d <- as.data.frame(brms::as_draws_df(fit_trials))
  pb <- .prior_blocks(fit_trials, d)
  expect_named(pb$blocks, c("loga", "logb", "sd"))
  expect_equal(pb$blocks$sd, dexp(d$sd_trial__loga_Intercept, 2, log = TRUE),
               tolerance = 1e-6)
})


test_that("a co-stressor's main effect and interaction make one block", {
  fit <- structure(list(
    prior = prior_table(c("normal(-7, 6)", "normal(0, 1)", "normal(0, 0.5)"),
                        coef = "Intercept", nlpar = c("loga", "cDO", "dTDO")),
    rtft = list(costressors = list(costressor("DO", "do", 8.5, 1, -1, interaction = TRUE)))),
    class = "brmsfit")
  d <- data.frame(b_loga_Intercept = c(-6, -8), b_cDO_Intercept = c(0.2, -1),
                  b_dTDO_Intercept = c(0.1, 0.3))
  d$lprior <- dnorm(d$b_loga_Intercept, -7, 6, log = TRUE) +
    dnorm(d$b_cDO_Intercept, 0, 1, log = TRUE) + dnorm(d$b_dTDO_Intercept, 0, 0.5, log = TRUE)
  pb <- .prior_blocks(fit, d)
  expect_named(pb$blocks, c("loga", "cDO"))
  expect_equal(pb$coefs$block, c("loga", "cDO", "cDO"))
  expect_equal(pb$blocks$cDO, dnorm(d$b_cDO_Intercept, 0, 1, log = TRUE) +
                 dnorm(d$b_dTDO_Intercept, 0, 0.5, log = TRUE))
})


test_that("the log prior by block refuses what it cannot read", {
  d <- as.data.frame(brms::as_draws_df(fit_exp))
  # no log prior among the draws
  expect_error(.prior_blocks(fit_exp, d[setdiff(names(d), "lprior")]), class = "rTFT_no_lprior")
  # a log prior the priors read from the fit do not add up to
  off <- d; off$lprior <- off$lprior + 0.01
  expect_error(.prior_blocks(fit_exp, off), class = "rTFT_prior_mismatch")
  # a coefficient of a parameter missing from the table of priors
  stray <- d; stray$b_logc_Intercept <- 0
  expect_error(.prior_blocks(fit_exp, stray), class = "rTFT_prior_mismatch")
  # a monotonic effect, whose priors the check does not read
  mono <- d; mono$bsp_loga_moage <- 0
  expect_error(.prior_blocks(fit_exp, mono), class = "rTFT_prior_unsupported")
  # a prior family the check cannot evaluate
  other <- fit_exp
  other$prior$prior[other$prior$nlpar == "loga" & other$prior$coef == ""] <- "cauchy(0, 1)"
  expect_error(.prior_blocks(other, d), class = "rTFT_prior_unsupported")
})


test_that("the log-likelihood summed over rows is the same in any batch of draws", {
  ll <- .summed_log_lik(fit_trials)
  expect_length(ll, brms::ndraws(fit_trials))
  expect_equal(ll, rowSums(brms::log_lik(fit_trials)))
  expect_equal(.summed_log_lik(fit_trials, batch = 7L), ll)
  expect_true(all(is.finite(ll) & ll < 0))
})
