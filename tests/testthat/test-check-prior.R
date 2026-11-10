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
  expect_equal(.coef_prior(as.data.frame(fit_exp$prior), "loga", "mass"),
               "normal(0, 1)")
})


test_that("the log density of a prior is that brms adds to the log prior", {
  x <- c(-3, 0, 1.5)
  expect_equal(.prior_lpdf(x, "normal(-0.7, 1)", "b"), dnorm(x, -0.7, 1, log = TRUE))
  expect_equal(.prior_lpdf(x, "normal(log(60), 2)", "b"),
               dnorm(x, log(60), 2, log = TRUE))
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
  expect_equal(pb$coefs$variable,
               c("b_loga_Intercept", "b_loga_mass", "b_logb_Intercept"))
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
    rtft = list(costressors = list(
      costressor("DO", "do", 8.5, 1, -1, interaction = TRUE)))),
    class = "brmsfit")
  d <- data.frame(b_loga_Intercept = c(-6, -8), b_cDO_Intercept = c(0.2, -1),
                  b_dTDO_Intercept = c(0.1, 0.3))
  d$lprior <- dnorm(d$b_loga_Intercept, -7, 6, log = TRUE) +
    dnorm(d$b_cDO_Intercept, 0, 1, log = TRUE) +
      dnorm(d$b_dTDO_Intercept, 0, 0.5, log = TRUE)
  pb <- .prior_blocks(fit, d)
  expect_named(pb$blocks, c("loga", "cDO"))
  expect_equal(pb$coefs$block, c("loga", "cDO", "cDO"))
  expect_equal(pb$blocks$cDO, dnorm(d$b_cDO_Intercept, 0, 1, log = TRUE) +
                 dnorm(d$b_dTDO_Intercept, 0, 0.5, log = TRUE))
})


test_that("the log prior by block refuses what it cannot read", {
  d <- as.data.frame(brms::as_draws_df(fit_exp))
  # no log prior among the draws
  expect_error(.prior_blocks(fit_exp, d[setdiff(names(d), "lprior")]),
               class = "rTFT_no_lprior")
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
  wide <- other$prior$nlpar == "loga" & other$prior$coef == ""
  other$prior$prior[wide] <- "cauchy(0, 1)"
  expect_error(.prior_blocks(other, d), class = "rTFT_prior_unsupported")
})


test_that("the log-likelihood summed over rows is the same in any batch of draws", {
  ll <- .summed_log_lik(fit_trials)
  expect_length(ll, brms::ndraws(fit_trials))
  expect_equal(ll, rowSums(brms::log_lik(fit_trials)))
  expect_equal(.summed_log_lik(fit_trials, batch = 7L), ll)
  expect_true(all(is.finite(ll) & ll < 0))
})


# ===========================================================================
# 3. The check
# ===========================================================================

test_that("the quantities checked are the coefficients, then each row's parameters", {
  d  <- as.data.frame(brms::as_draws_df(fit_exp))
  pb <- .prior_blocks(fit_exp, d)
  nd <- data.frame(mass = c(8, 12))
  qs <- .check_quantities(fit_exp, d, pb$coefs, nd, p = 0.5, t_ref = 60, levels = NULL)
  expect_equal(qs$index$quantity,
               c(pb$coefs$variable, rep(c("a", "b", "z", "sCTmax"), 2)))
  expect_equal(qs$index$row, c(NA, NA, NA, rep(1:2, each = 4)))
  expect_equal(qs$values[[1]], d$b_loga_Intercept)
  # a at mass 12, from the same draws in the same order
  expect_equal(qs$values[[8]], exp(d$b_loga_Intercept + 12 * d$b_loga_mass))
  # sCTmax as derive_tdt_ttl() reports it, from every draw
  sct <- derive_tdt_ttl(fit_exp, newdata = nd, p = 0.5, t_ref = 60,
                        n_draws = NULL)$sctmax
  expect_equal(c(median(qs$values[[7]]), median(qs$values[[11]])), sct$median)
  # without newdata, the coefficients alone
  expect_length(.check_quantities(fit_exp, d, pb$coefs, NULL, 0.5, NULL, NULL)$values,
                3)
})


test_that("each prior is set beside the posterior of its coefficient", {
  d  <- as.data.frame(brms::as_draws_df(fit_exp))
  pb <- .prior_blocks(fit_exp, d)
  pr <- .prior_table(pb$coefs, d)
  expect_named(pr, c("parameter", "coefficient", "prior", "prior_centre", "prior_scale",
                     "median", "sd", "distance", "sd_ratio"))
  expect_equal(pr$coefficient, c("Intercept", "mass", "Intercept"))
  expect_equal(pr$prior_centre, c(-7.012, 0, -0.7))
  expect_equal(pr$prior_scale, c(6.36, 1, 1))
  expect_equal(pr$median[2], median(d$b_loga_mass))
  expect_equal(pr$distance, (pr$median - pr$prior_centre) / pr$prior_scale)
  expect_equal(pr$sd_ratio, pr$sd / pr$prior_scale)
  # a Student-t gives its location and scale, a flat prior neither
  coefs <- data.frame(variable = c("b_loga_mass", "b_loga_mass"), nlpar = "loga",
                      coef = "mass", prior = c("student_t(3, 1, 2)", ""))
  pr <- .prior_table(coefs, d)
  expect_equal(pr$prior_centre, c(1, NA))
  expect_equal(pr$prior_scale, c(2, NA))
  expect_true(is.na(pr$distance[2]))
})


test_that("the check gives each block's sensitivity for each quantity", {
  skip_if_not_installed("priorsense")
  chk <- check_prior(fit_exp, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
  expect_named(chk, c("sensitivity", "priors"))
  s <- chk$sensitivity
  expect_named(s, c("block", "quantity", "mass", "prior", "likelihood", "diagnosis",
                    "k_prior", "k_likelihood"))
  expect_equal(unique(s$block), c("all", "loga", "logb"))
  expect_equal(nrow(s), 3 * 11)                 # 3 coefficients, 2 rows of 4
  expect_true(all(is.na(s$mass[startsWith(s$quantity, "b_")])))
  expect_equal(s$mass[s$quantity == "sCTmax"], rep(c(8, 12), 3))
  expect_true(all(is.finite(s$prior) & s$prior >= 0))
  expect_true(all(is.finite(s$likelihood) & s$likelihood >= 0))
  # the likelihood is scaled the same way whichever block of priors is
  expect_equal(s$likelihood[s$block == "loga"], s$likelihood[s$block == "all"])
  expect_equal(nrow(chk$priors), 3)
  # the coefficients alone, without newdata
  only <- check_prior(fit_exp)$sensitivity
  expect_named(only, c("block", "quantity", "prior", "likelihood", "diagnosis",
                       "k_prior", "k_likelihood"))
  expect_equal(nrow(only), 9)
  # a model with a trial SD has the block sd; a model without covariates takes one row
  tr <- check_prior(fit_trials, newdata = data.frame(row.names = 1),
                    t_ref = 30)$sensitivity
  expect_equal(unique(tr$block), c("all", "loga", "logb", "sd"))
  expect_true("sCTmax" %in% tr$quantity)
})


test_that("a quantity with infinite draws has no sensitivity", {
  skip_if_not_installed("priorsense")
  d  <- as.data.frame(brms::as_draws_df(fit_exp))
  pb <- .prior_blocks(fit_exp, d)
  never <- d$b_loga_Intercept
  never[1] <- Inf
  expect_message(
    s <- .sensitivity_table(list(d$b_logb_Intercept, never), pb$blocks, d$lprior,
                            .summed_log_lik(fit_exp), d),
    class = "rTFT_unreached")
  expect_true(all(is.na(s$prior[s$q == 2])))
  expect_true(all(is.finite(s$prior[s$q == 1])))
})


test_that("the check refuses what it cannot do", {
  old <- fit_exp
  old$rtft <- NULL
  expect_error(check_prior(old), class = "rTFT_no_settings")
  local_mocked_bindings(.has_priorsense = function() FALSE)
  expect_error(check_prior(fit_exp), class = "rTFT_needs_priorsense")
  local_mocked_bindings(.has_priorsense = function() TRUE)
  nd <- data.frame(mass = 10)
  expect_error(check_prior(fit_exp, t_ref = 60), class = "rTFT_needs_newdata")
  expect_error(check_prior(fit_exp, nd, t_ref = -1), class = "rTFT_bad_t_ref")
  expect_error(check_prior(fit_exp, nd, t_ref = c(30, 60)), class = "rTFT_bad_t_ref")
  expect_error(check_prior(fit_exp, nd, t_ref = 60, p = 1), class = "rTFT_bad_p")
  expect_error(check_prior(fit_exp, nd, t_ref = 60, p = c(0.1, 0.5)),
               class = "rTFT_bad_p")
  expect_error(check_prior(fit_exp, data.frame(size = 10)),
               class = "rTFT_missing_covariates")
  # p is checked whether or not t_ref is given
  expect_error(check_prior(fit_exp, nd, p = 2), class = "rTFT_bad_p")
  # a column of newdata may not take a name of the table returned, nor hold NA
  expect_error(check_prior(fit_exp, data.frame(mass = 10, prior = "x")),
               class = "rTFT_newdata_names")
  expect_error(check_prior(fit_exp, data.frame(mass = NA_real_)),
               class = "rTFT_missing_values")
  # levels without newdata reach nothing
  expect_warning(check_prior(fit_exp, levels = c(do = 50)),
                 class = "rTFT_unused_levels")
})


test_that("the tolerance on the rebuilt log prior follows the rounding of the draws", {
  d <- as.data.frame(brms::as_draws_df(fit_exp))
  # draws stored to six significant figures, as CmdStan writes them at its coarsest
  rounded <- d
  for (v in c("b_loga_Intercept", "b_loga_mass", "b_logb_Intercept", "lprior"))
    rounded[[v]] <- signif(d[[v]], 6)
  expect_no_error(.prior_blocks(fit_exp, rounded))
  # the tolerance is 1e-5 of the size of the log prior (here about 4.7)
  near <- d; near$lprior <- d$lprior + 3e-5
  expect_no_error(.prior_blocks(fit_exp, near))
  far <- d; far$lprior <- d$lprior + 1e-4
  expect_error(.prior_blocks(fit_exp, far), class = "rTFT_prior_mismatch")
})


test_that("a block of flat priors is left out of the scaling, with a message", {
  skip_if_not_installed("priorsense")
  d  <- as.data.frame(brms::as_draws_df(fit_exp))
  pb <- .prior_blocks(fit_exp, d)
  ll <- .summed_log_lik(fit_exp)
  flat <- pb$blocks
  flat$logb <- rep(0, nrow(d))
  expect_message(
    s <- .sensitivity_table(list(d$b_loga_Intercept), flat, d$lprior - pb$blocks$logb,
                            ll, d),
    class = "rTFT_flat_block")
  expect_equal(unique(s$block), c("all", "loga"))
  # with every prior flat there is nothing to scale at all
  expect_message(
    none <- .sensitivity_table(list(d$b_loga_Intercept), list(loga = rep(0, nrow(d))),
                               rep(0, nrow(d)), ll, d),
    class = "rTFT_flat_block")
  expect_equal(nrow(none), 0)
  expect_named(none, c("block", "q", "prior", "likelihood", "diagnosis", "k_prior",
                       "k_likelihood"))
})


test_that("each scaling reports its Pareto k and warns when it is unreliable", {
  skip_if_not_installed("priorsense")
  d  <- as.data.frame(brms::as_draws_df(fit_exp))
  pb <- .prior_blocks(fit_exp, d)
  ll <- .summed_log_lik(fit_exp)
  s  <- .sensitivity_table(list(d$b_loga_Intercept), pb$blocks, d$lprior, ll, d)
  threshold <- min(1 - 1 / log10(nrow(d)), 0.7)
  expect_true(all(is.finite(s$k_prior)))
  expect_equal(length(unique(s$k_likelihood)), 1L)
  expect_true(all(c(s$k_prior, s$k_likelihood) < threshold))
  # a block so steep across the draws that a few of them carry all the weight
  steep <- pb$blocks
  steep$loga <- 1e3 * d$b_loga_Intercept
  expect_warning(
    w <- .sensitivity_table(list(d$b_loga_Intercept), steep, d$lprior, ll, d),
    class = "rTFT_unreliable_scaling")
  expect_true(all(w$k_prior[w$block == "loga"] > threshold))
  expect_true(all(w$k_prior[w$block == "logb"] < threshold))
})
