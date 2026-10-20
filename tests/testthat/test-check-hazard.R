# ===========================================================================
# 1. Helpers
# ===========================================================================

# Two fitted models: `fit_exp` has 40 individuals on one ramp, with body mass as a
# covariate, and `fit_trials` has eight trials of six individuals on one ramp, with a
# random intercept per trial. Both are interval-censored.
fit_exp    <- readRDS(test_path("fixtures", "hazard-fit-exp.rds"))
fit_trials <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))


# ===========================================================================
# 2. Following individuals through the fitted data
# ===========================================================================

test_that(".fit_record gives each individual its rows, time and outcome", {
  rec <- .fit_record(fit_exp)
  st  <- fit_exp$rtft
  expect_equal(nrow(rec$individuals), 40L)
  expect_equal(rec$individuals$id, unique(st$id))
  expect_true(all(rec$individuals$path_id == "ramp"))
  # the time is the end of the individual's last interval, and the outcome is in it
  expect_equal(rec$individuals$time, as.vector(tapply(st$tau_end, st$id, max)))
  expect_equal(rec$individuals$event, as.vector(tapply(fit_exp$data$failed, st$id, max)))
  expect_equal(st$tau_end[rec$individuals$last], rec$individuals$time)
  # each individual's rows are in order of time, and together they are every row
  expect_true(all(vapply(rec$rows, \(r) !is.unsorted(st$tau_end[r]), logical(1))))
  expect_setequal(unlist(rec$rows), seq_along(st$id))
  # a fit with no record of the intervals' times cannot be checked
  old <- fit_exp
  old$rtft$tau_end <- NULL
  expect_error(.fit_record(old), class = "rTFT_no_check_record")
  expect_error(.fit_record(list()), class = "rTFT_not_brmsfit")
})


test_that(".individual_columns tells what an individual carries from what changes along its exposure", {
  rec  <- .fit_record(fit_trials)
  cols <- .individual_columns(rec$data, rec$rows)
  expect_true(all(c("trial", "Tr") %in% cols))
  # the temperature changes along a ramp, and the response is never carried
  expect_false(any(c("Ta", "failed") %in% cols))
  # on a hand-made table: `mass` is one value per individual, `Ta` is not
  d <- data.frame(mass = c(8, 8, 12, 12), Ta = c(20, 21, 20, 21), failed = c(0, 1, 0, 0))
  expect_equal(.individual_columns(d, list(1:2, 3:4)), "mass")
})


# ===========================================================================
# 3. The posterior predictive check of survival
# ===========================================================================

test_that(".km_survival is the Kaplan-Meier estimate, allowing for those no longer observed", {
  time  <- c(5, 10, 10, 20)
  event <- c(1, 1, 0, 1)
  # at 5: 3 of 4 remain; at 10: 1 of the 3 at risk fails; at 20: the last one fails
  expect_equal(.km_survival(time, event, at = c(1, 5, 10, 15, 20)),
               c(1, 0.75, 0.5, 0.5, 0))
  expect_equal(.km_survival(time, rep(0, 4), at = c(5, 20)), c(1, 1))   # no failures
  expect_equal(.km_survival(c(3, 3), c(1, 1), at = c(2, 3)), c(1, 0))
})


test_that(".replicate_survival follows each individual until it fails", {
  # one draw, two individuals over three intervals: the first is certain to fail in
  # its second interval, the second cannot fail
  p <- matrix(c(0, 1, 0, 0, 0, 0), nrow = 1)
  alive <- .replicate_survival(p, n_individuals = 2, n_intervals = 3)
  expect_equal(dim(alive), c(1L, 2L, 3L))
  expect_equal(alive[1, 1, ], c(TRUE, FALSE, FALSE))     # failed stays failed
  expect_equal(alive[1, 2, ], c(TRUE, TRUE, TRUE))
  expect_true(all(.replicate_survival(matrix(0, 5, 12), 3, 4)))
  expect_false(any(.replicate_survival(matrix(1, 5, 12), 3, 4)))
  # with chance in it, survival never comes back
  set.seed(1)
  some <- .replicate_survival(matrix(0.3, 50, 40), 4, 10)
  expect_true(all(apply(some, c(1, 2), \(x) !is.unsorted(rev(x)))))
  expect_equal(mean(some[, , 1]), 0.7, tolerance = 0.1)
})


test_that("check_survival returns observed and simulated survival through time", {
  set.seed(1)
  chk <- check_survival(fit_exp, n_draws = 100)
  expect_named(chk, c("path_id", "n", "time", "Ta", "observed", "median", "lower",
                      "upper", "prob"))
  expect_true(all(chk$path_id == "ramp"))
  expect_true(all(chk$n == 40))
  # it starts with everyone unaffected, and follows the ramp's checks to the last
  # time any individual was observed
  expect_equal(unlist(chk[1, c("time", "observed", "median", "lower", "upper")],
                      use.names = FALSE), c(0, 1, 1, 1, 1))
  rec <- .fit_record(fit_exp)
  expect_equal(max(chk$time), max(rec$individuals$time))
  expect_equal(diff(chk$time), rep(5, nrow(chk) - 1))
  # the observed proportion is the Kaplan-Meier estimate of the individuals' times
  expect_equal(chk$observed[-1],
               .km_survival(rec$individuals$time, rec$individuals$event, chk$time[-1]))
  expect_true(all(diff(chk$observed) <= 0))
  expect_true(all(chk$lower <= chk$median & chk$median <= chk$upper))
  expect_true(all(diff(chk$median) <= 1e-12))
  # the model the data were simulated from covers what was observed
  expect_gt(mean(chk$observed >= chk$lower & chk$observed <= chk$upper), 0.9)
  # the temperature is that of the ramp in each interval
  expect_equal(chk$Ta[-1], 20 + 0.25 * (chk$time[-1] - 2.5))
})


test_that("check_survival is the same again under the same seed, and narrower with prob", {
  set.seed(7); a <- check_survival(fit_exp, n_draws = 50)
  set.seed(7); b <- check_survival(fit_exp, n_draws = 50)
  expect_equal(a, b)
  set.seed(7); half <- check_survival(fit_exp, prob = 0.5, n_draws = 50)
  expect_true(all(half$upper - half$lower <= a$upper - a$lower + 1e-12))
  expect_equal(unique(half$prob), 0.5)
})


test_that("by splits the individuals of a trajectory into groups", {
  set.seed(1)
  chk <- check_survival(fit_trials, by = "trial", n_draws = 50)
  expect_true("trial" %in% names(chk))
  groups <- unique(chk[c("path_id", "trial", "n")])
  expect_equal(nrow(groups), 8L)
  expect_true(all(groups$n == 6))
  # every group is followed over the same times, those of its trajectory
  expect_equal(length(unique(table(chk$trial))), 1L)
  # a group of six can only take sevenths... of six: proportions are multiples of 1/6
  expect_true(all(abs(chk$median * 6 - round(chk$median * 6)) < 1e-8 |
                    abs(chk$median * 12 - round(chk$median * 12)) < 1e-8))
  rec <- .fit_record(fit_trials)
  first <- rec$individuals[rec$data$trial[rec$individuals$last] == groups$trial[1], ]
  one <- chk[chk$trial == groups$trial[1], ]
  expect_equal(one$observed[-1], .km_survival(first$time, first$event, one$time[-1]))
})


test_that("check_survival refuses bad arguments", {
  expect_error(check_survival(fit_exp, by = "sex"), class = "rTFT_bad_by")
  expect_error(check_survival(fit_exp, by = "Ta"), class = "rTFT_bad_by")     # changes in time
  expect_error(check_survival(fit_exp, by = c("mass", "mass")), class = "rTFT_bad_by")
  expect_error(check_survival(fit_exp, prob = 1.5), class = "rTFT_bad_prob")
  expect_error(check_survival(list()), class = "rTFT_not_brmsfit")
})


# ===========================================================================
# 4. Residuals by individual
# ===========================================================================

test_that("hazard_residuals gives one residual per individual, with what it carries", {
  set.seed(1)
  res <- hazard_residuals(fit_exp, n_draws = 200)
  expect_named(res, c("id", "path_id", "mass", "time", "Ta", "event", "randomisation",
                      "residual"))
  expect_true(all(res$randomisation == 1))
  expect_equal(nrow(res), 40L)
  rec <- .fit_record(fit_exp)
  expect_equal(res$id, rec$individuals$id)
  expect_equal(res$time, rec$individuals$time)
  expect_equal(res$event, rec$individuals$event)
  expect_true(all(res$residual > 0 & res$residual < 1))
  with_trial <- hazard_residuals(fit_trials, n_draws = 100)
  expect_true("trial" %in% names(with_trial))
  expect_equal(nrow(with_trial), 48L)
})


test_that("each residual lies within what was observed of the individual", {
  # the fraction the model predicts to have failed by the end of each individual's
  # last interval and by its start, by hand from the linear predictor
  ids <- .thin_ids(brms::ndraws(fit_exp), 200)
  accrued <- exp(brms::posterior_linpred(fit_exp, draw_ids = ids))
  rec <- .fit_record(fit_exp)
  F_end <- vapply(rec$rows, \(r)
    1 - mean(exp(-rowSums(accrued[, r, drop = FALSE]))), numeric(1))
  F_start <- vapply(rec$rows, \(r)
    1 - mean(exp(-rowSums(accrued[, r[-length(r)], drop = FALSE]))), numeric(1))
  # the fraction predicted to have failed rises along an individual's exposure
  expect_true(all(F_start < F_end))
  # every individual of the fixture failed; the first five are made survivors, last
  # seen at the end of the interval in which they in fact failed
  some_survive <- fit_exp
  some_survive$data$failed[rec$individuals$last[1:5]] <- 0L
  set.seed(3)
  res <- hazard_residuals(some_survive, n_draws = 200)
  failed <- res$event == 1
  expect_equal(sum(!failed), 5L)
  # a failure known to an interval: between the fractions at its two ends
  expect_true(all(res$residual[failed] >= F_start[failed] &
                    res$residual[failed] <= F_end[failed]))
  # a survivor: above the fraction at its last observation
  expect_true(all(res$residual[!failed] >= F_end[!failed]))
  # a failure at a known time: the fraction at that time itself
  exact <- some_survive
  exact$rtft$observation <- "continuous"
  set.seed(3)
  res_exact <- hazard_residuals(exact, n_draws = 200)
  expect_equal(res_exact$residual[failed], unname(F_end[failed]))
  expect_true(all(res_exact$residual[!failed] >= F_end[!failed]))
})


test_that("the residuals are the same again under the same seed, and uniform for the right model", {
  set.seed(11); a <- hazard_residuals(fit_exp)
  set.seed(11); b <- hazard_residuals(fit_exp)
  expect_equal(a, b)
  # the fixture was simulated from the model fitted to it
  expect_gt(stats::ks.test(a$residual, "punif")$p.value, 0.05)
  expect_error(hazard_residuals(list()), class = "rTFT_not_brmsfit")
})


test_that("the residuals can be drawn several times, the first as a single draw gives it", {
  set.seed(5); one   <- hazard_residuals(fit_exp, n_draws = 200)
  set.seed(5); three <- hazard_residuals(fit_exp, n_draws = 200, n_randomisations = 3)
  # every individual once in each randomisation, with what it carries
  expect_equal(nrow(three), 3L * nrow(one))
  expect_equal(three$randomisation, rep(1:3, each = nrow(one)))
  for (k in 1:3) {
    of_k <- three[three$randomisation == k, ]
    expect_equal(of_k[setdiff(names(of_k), c("randomisation", "residual"))],
                 one[setdiff(names(one), c("randomisation", "residual"))])
  }
  # under the same seed, the first randomisation is the single one
  expect_equal(three$residual[three$randomisation == 1], one$residual)
  # failures known only to an interval are drawn again in each randomisation
  r <- matrix(three$residual, ncol = 3)
  expect_true(all(r[, 1] != r[, 2] & r[, 2] != r[, 3]))
  # a failure at a known time is the same in every randomisation, and a survivor is
  # drawn again
  rec <- .fit_record(fit_exp)
  exact <- fit_exp
  exact$rtft$observation <- "continuous"
  exact$data$failed[rec$individuals$last[1:5]] <- 0L
  set.seed(5)
  r <- matrix(hazard_residuals(exact, n_draws = 200, n_randomisations = 3)$residual,
              ncol = 3)
  expect_equal(r[-(1:5), 1], r[-(1:5), 2])
  expect_equal(r[-(1:5), 1], r[-(1:5), 3])
  expect_true(all(r[1:5, 1] != r[1:5, 2]))
  # the number of randomisations is a whole number, 1 or more
  for (bad in list(0, 1.5, NA_real_, Inf, c(1, 2), "2"))
    expect_error(hazard_residuals(fit_exp, n_randomisations = bad),
                 class = "rTFT_bad_n_randomisations")
})


# ===========================================================================
# 5. A fitted model of static and ramping trials together
# ===========================================================================

# A small fitted model saved as a fixture (made by fixtures/make-hazard-fit-mixed.R):
# three static trials checked every 5 minutes, whose failures are known to the
# interval between two checks, and two ramps watched throughout, whose failures are
# at known times, with ten individuals in each.
fit_mixed <- readRDS(test_path("fixtures", "hazard-fit-mixed.rds"))

test_that("check_survival follows the static trials at their checks and the ramps throughout", {
  set.seed(1)
  chk <- check_survival(fit_mixed, n_draws = 100)
  rec <- .fit_record(fit_mixed)
  expect_setequal(chk$path_id, c("T23", "T25", "T27", "r0.25", "r0.1"))
  expect_true(all(chk$n == 10))
  # every trial starts with everyone unaffected
  expect_true(all(chk$observed[chk$time == 0] == 1))
  expect_true(all(chk$lower <= chk$median & chk$median <= chk$upper))
  # a static trial is followed at its checks and stays at its temperature
  of_static <- grepl("^T", chk$path_id)
  expect_true(all(chk$time[of_static] %% 5 == 0))
  expect_equal(chk$Ta[of_static & chk$time > 0],
               as.numeric(sub("T", "", chk$path_id[of_static & chk$time > 0])))
  # a ramp warms through the times at which it is followed
  expect_true(all(diff(chk$Ta[chk$path_id == "r0.25" & chk$time > 0]) > 0))
  # in each trial the observed proportion is the Kaplan-Meier estimate of the times
  # of its own individuals, whichever way their failures were observed
  for (p in unique(chk$path_id)) {
    rows <- chk$path_id == p & chk$time > 0
    its  <- rec$individuals[rec$individuals$path_id == p, ]
    expect_equal(chk$observed[rows], .km_survival(its$time, its$event, chk$time[rows]))
  }
})


test_that("hazard_residuals gives one residual per individual of either kind of assay", {
  set.seed(1)
  res <- hazard_residuals(fit_mixed, n_draws = 200)
  rec <- .fit_record(fit_mixed)
  expect_equal(nrow(res), 50L)
  expect_equal(res$id, rec$individuals$id)
  expect_equal(res$event, rec$individuals$event)
  expect_true(all(res$residual > 0 & res$residual < 1))
  # the fixture was simulated from the model fitted to it
  expect_gt(stats::ks.test(res$residual, "punif")$p.value, 0.05)
})
