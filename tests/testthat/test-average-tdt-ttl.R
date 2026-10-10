# ===========================================================================
# 1. Helpers
# ===========================================================================

# Two models fitted to the same data (40 individuals, log a depending on mass): one
# with the exponential failure rate and one with the zero-threshold exponential.
fits <- list(exp     = readRDS(test_path("fixtures", "hazard-fit-exp.rds")),
             zeroexp = readRDS(test_path("fixtures", "hazard-fit-zeroexp.rds")))

# A hand-made cross-validation result: the log predictive density of each unit, in
# the layout loo_hazard() (column elpd_loo) or kfold_hazard() (elpd_kfold) returns.
make_cv <- function(elpd, unit = "id", column = "elpd_loo") {
  pointwise <- matrix(elpd, ncol = 1, dimnames = list(names(elpd), column))
  structure(list(pointwise = pointwise, unit = unit), class = "loo")
}

# Three hand-made models over 30 units. Each is the best predictor for a third of
# the units, most often `good` and least often `poor`, so all three earn weight.
set.seed(20261005)
units <- paste0("u", 1:30)
base  <- stats::setNames(rnorm(30, -2, 0.3), units)
cv3 <- list(poor   = make_cv(base - c(rep(1.5, 24), rep(0, 6))),
            good   = make_cv(base - c(rep(0, 16), rep(1.5, 14))),
            middle = make_cv(base - c(rep(1.5, 16), rep(0, 8), rep(1.5, 6))))

quietly <- function(expr) suppressMessages(expr)


# ===========================================================================
# 2. Model weights
# ===========================================================================

test_that(".pointwise_elpd reads either kind of cross-validation result", {
  elpd <- c(a = -1, b = -2)
  expect_equal(.pointwise_elpd(make_cv(elpd)), elpd)
  expect_equal(.pointwise_elpd(make_cv(elpd, column = "elpd_kfold")), elpd)
  expect_null(.pointwise_elpd(list(pointwise = matrix(1))))        # not a result
  expect_null(.pointwise_elpd(make_cv(unname(elpd))))              # units not named
  expect_equal(names(.pointwise_elpd(loo_hazard(fits$exp))),
               as.character(sort(unique(fits$exp$rtft$id))))
})


test_that("hazard_weights ranks the models and reports what it weighed", {
  expect_message(w <- hazard_weights(cv3), class = "rTFT_weights")
  expect_named(w, c("model", "elpd", "weight_all", "cumulative", "kept", "weight"))
  expect_equal(w$model, c("good", "middle", "poor"))
  expect_equal(sum(w$weight_all), 1)
  expect_true(all(diff(w$weight_all) <= 0))
  expect_equal(w$cumulative, cumsum(w$weight_all))
  expect_true(all(w$kept))
  expect_equal(w$weight, w$weight_all)
  expect_equal(w$elpd[w$model == "good"], sum(cv3$good$pointwise))
  expect_equal(attr(w, "unit"), "id")
  expect_equal(attr(w, "n_units"), 30)
  expect_equal(attr(w, "method"), "stacking")
})


test_that("the weights are those of loo for the same pointwise values", {
  lpd <- vapply(cv3, \(x) x$pointwise[, 1], numeric(30))
  w <- quietly(hazard_weights(cv3))
  expect_equal(w$weight_all, as.numeric(loo::stacking_weights(lpd))[c(2, 3, 1)])
  set.seed(1); w_bma <- quietly(hazard_weights(cv3, method = "pseudobma"))
  set.seed(1); direct <- as.numeric(loo::pseudobma_weights(lpd))
  expect_equal(w_bma$weight_all, direct[match(w_bma$model, names(cv3))])
  expect_equal(attr(w_bma, "method"), "pseudobma")
})


test_that("units are matched by name, whatever their order", {
  shuffled <- cv3
  shuffled$good <- make_cv(cv3$good$pointwise[30:1, 1])
  expect_equal(quietly(hazard_weights(shuffled)), quietly(hazard_weights(cv3)))
})


test_that("min_weight keeps the leading models and weights them among themselves", {
  all3 <- quietly(hazard_weights(cv3))
  # just above the leading weight: the second model is needed, the third is not
  w <- quietly(hazard_weights(cv3, min_weight = all3$weight_all[1] + 0.01))
  expect_equal(w$kept, c(TRUE, TRUE, FALSE))
  expect_equal(w$weight_all, all3$weight_all)          # the ranking weights are unchanged
  expect_equal(w$weight[3], 0)
  expect_equal(sum(w$weight), 1)
  lpd <- vapply(cv3[c("good", "middle")], \(x) x$pointwise[, 1], numeric(30))
  expect_equal(w$weight[1:2], as.numeric(loo::stacking_weights(lpd)))
  # recomputed among the two, not rescaled from the weights among the three
  expect_false(isTRUE(all.equal(w$weight[1:2], all3$weight_all[1:2] / sum(all3$weight_all[1:2]))))
  # a threshold the leading model meets alone leaves it with all the weight
  one <- quietly(hazard_weights(cv3, min_weight = all3$weight_all[1]))
  expect_equal(one$kept, c(TRUE, FALSE, FALSE))
  expect_equal(one$weight, c(1, 0, 0))
})


test_that("the message says how many units the weights rest on, and of what kind", {
  few <- purrr::map(cv3, \(x) make_cv(x$pointwise[c(1:3, 17:19, 25:27), 1], unit = "trial"))
  expect_message(hazard_weights(few), "9 units left out \\(one trial each\\)")
  expect_message(hazard_weights(cv3), "30 units left out \\(one individual each\\)")
})


test_that("hazard_weights refuses what cannot be weighted", {
  expect_error(hazard_weights(cv3[1]), class = "rTFT_bad_cv")                 # one model
  expect_error(hazard_weights(unname(cv3)), class = "rTFT_bad_cv")            # no names
  expect_error(hazard_weights(cv3[c(1, 1)]), class = "rTFT_bad_cv")           # a name twice
  expect_error(hazard_weights(list(a = cv3$good, b = fits$exp)), class = "rTFT_bad_cv")
  expect_error(hazard_weights(cv3, min_weight = 0), class = "rTFT_bad_min_weight")
  expect_error(hazard_weights(cv3, min_weight = 1.2), class = "rTFT_bad_min_weight")
  expect_error(hazard_weights(cv3, method = "aic"))
  # a different unit, and the same unit over different units left out
  other_unit <- cv3; other_unit$good$unit <- "trial"
  expect_error(hazard_weights(other_unit), class = "rTFT_cv_mismatch")
  fewer <- cv3; fewer$good <- make_cv(cv3$good$pointwise[-1, 1])
  expect_error(hazard_weights(fewer), class = "rTFT_cv_mismatch")
  # a single unit left out cannot weight models
  one_unit <- lapply(cv3, \(x) make_cv(stats::setNames(x$pointwise[1, 1], rownames(x$pointwise)[1])))
  expect_error(hazard_weights(one_unit), class = "rTFT_too_few_units")
})


test_that("hazard_weights weighs fitted models, by PSIS or by refitting", {
  cv <- purrr::map(fits, loo_hazard)
  w <- quietly(hazard_weights(cv))
  expect_setequal(w$model, c("exp", "zeroexp"))
  expect_equal(sum(w$weight), 1)
  expect_equal(attr(w, "n_units"), 40)
  mixed <- list(a = make_cv(base), b = make_cv(base - 0.1, column = "elpd_kfold"))
  expect_equal(quietly(hazard_weights(mixed))$model[1], "a")
})


test_that("models cross-validated with fewer folds than the others are named", {
  with_K <- function(x, K) { x$K <- K; x }
  k_fold <- purrr::map(cv3, \(x) make_cv(x$pointwise[, 1], column = "elpd_kfold"))
  # leave-one-out for two models, ten folds over the 30 units for the third
  mixed <- list(poor = cv3$poor, good = cv3$good, middle = with_K(k_fold$middle, 10))
  expect_message(suppressMessages(hazard_weights(mixed), classes = "rTFT_weights"),
                 "middle", class = "rTFT_mixed_cv")
  # one fold per unit estimates what leave-one-out does: nothing to say
  as_loo <- list(poor = cv3$poor, good = cv3$good, middle = with_K(k_fold$middle, 30))
  expect_false(any(grepl("number of folds", testthat::capture_messages(hazard_weights(as_loo)))))
  # all refitted with the same number of folds: nothing to say either
  same_K <- purrr::map(k_fold, with_K, K = 10)
  expect_false(any(grepl("number of folds", testthat::capture_messages(hazard_weights(same_K)))))
  expect_equal(suppressMessages(hazard_weights(mixed))$weight,
               suppressMessages(hazard_weights(cv3))$weight)       # the weights are not altered
})


# ===========================================================================
# 3. Averaging what is derived
# ===========================================================================

test_that(".allocate_draws shares the draws in proportion, as whole numbers", {
  expect_equal(.allocate_draws(c(a = 0.6, b = 0.3, c = 0.1), 1000), c(a = 600L, b = 300L, c = 100L))
  thirds <- .allocate_draws(c(a = 1, b = 1, c = 1) / 3, 10)
  expect_equal(sum(thirds), 10L)
  expect_true(all(thirds %in% 3:4))
  expect_equal(.allocate_draws(c(a = 1, b = 0), 7), c(a = 7L, b = 0L))
  expect_equal(.allocate_draws(c(a = 2, b = 6), 8), c(a = 2L, b = 6L))    # need not sum to 1
})


test_that(".pool_draws puts the models' draws one after another", {
  one_model <- function(value, n) list(
    newdata = data.frame(mass = 10), p = 0.5, t_ref = 60, Ta = c(20, 22), tau = c(1, 10, 100),
    rows = list(list(params  = list(a = rep(value, n)),
                     sctmax  = list(rep(value, n)),
                     tdt     = list(matrix(value, n, 2)),
                     surface = list(matrix(value, n, 3), matrix(value + 0.5, n, 3)))))
  pooled <- .pool_draws(list(one_model(1, 3), one_model(2, 2)))
  expect_equal(pooled$rows[[1]]$sctmax[[1]], c(1, 1, 1, 2, 2))
  expect_equal(dim(pooled$rows[[1]]$tdt[[1]]), c(5, 2))
  expect_equal(pooled$rows[[1]]$surface[[2]][, 1], c(1.5, 1.5, 1.5, 2.5, 2.5))
  expect_null(pooled$rows[[1]]$params)
  expect_equal(pooled$Ta, c(20, 22))
})


nd <- data.frame(mass = c(8, 12))
w_fits <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0.75, 0.25))

test_that("average_tdt_ttl returns the tables of derive_tdt_ttl, without parameters", {
  avg <- average_tdt_ttl(fits, w_fits, newdata = nd, t_ref = 60, n_draws = 400)
  one <- derive_tdt_ttl(fits$exp, newdata = nd, t_ref = 60, n_draws = 400)
  expect_named(avg, c("sctmax", "tdt_curves", "ttl_surface", "weights"))
  expect_named(avg$sctmax, names(one$sctmax))
  expect_named(avg$tdt_curves, names(one$tdt_curves))
  expect_named(avg$ttl_surface, names(one$ttl_surface))
  expect_equal(nrow(avg$tdt_curves), nrow(one$tdt_curves))
  expect_equal(avg$ttl_surface[c("mass", "Ta", "tau")], one$ttl_surface[c("mass", "Ta", "tau")])
  expect_equal(avg$weights$n_draws, c(300L, 100L))
  expect_true(all(avg$sctmax$n_draws == 400))
})


test_that("a model with all the weight gives back its own results", {
  alone <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0, 1))
  avg <- average_tdt_ttl(fits, alone, newdata = nd, t_ref = 60, n_draws = 300)
  one <- derive_tdt_ttl(fits$zeroexp, newdata = nd, t_ref = 60, n_draws = 300,
                        Ta = unique(avg$tdt_curves$Ta), tau = unique(avg$ttl_surface$tau))
  expect_equal(avg$sctmax, one$sctmax)
  expect_equal(avg$tdt_curves, one$tdt_curves)
  expect_equal(avg$ttl_surface, one$ttl_surface)
  expect_equal(avg$weights$n_draws, c(0L, 300L))
})


test_that("the average is the summary of the models' draws pooled by weight", {
  avg <- quietly(average_tdt_ttl(fits, w_fits, newdata = nd[1, , drop = FALSE], p = 0.5,
                                 t_ref = 60, n_draws = 400))
  grid <- .fit_settings(fits$exp)
  Ta  <- seq(grid$Ta_range[1], grid$Ta_range[2], length.out = 50)
  tau <- exp(seq(log(grid$tau_range[1]), log(grid$tau_range[2]), length.out = 50))
  d_exp  <- .derive_draws(fits$exp,     nd[1, , drop = FALSE], 0.5, 60, Ta, tau, NULL, 300)
  d_zero <- .derive_draws(fits$zeroexp, nd[1, , drop = FALSE], 0.5, 60, Ta, tau, NULL, 100)
  pooled <- c(d_exp$rows[[1]]$sctmax[[1]], d_zero$rows[[1]]$sctmax[[1]])
  expect_equal(avg$sctmax$median, stats::median(pooled))
  expect_equal(avg$sctmax$lower, unname(stats::quantile(pooled, 0.025)))
})


test_that("the pool shrinks when a model has too few draws for its share", {
  # each fit holds 600 draws; three quarters of 1,000 would need 750 of one of them
  expect_message(
    avg <- average_tdt_ttl(fits, w_fits, newdata = nd[1, , drop = FALSE], p = 0.5, t_ref = 60,
                           n_draws = 1000),
    class = "rTFT_pool_reduced")
  expect_equal(avg$weights$n_draws, c(600L, 200L))
  expect_true(all(avg$sctmax$n_draws == 800))
})


test_that("weights from hazard_weights, and extra fits, are accepted", {
  w <- quietly(hazard_weights(purrr::map(fits, loo_hazard), min_weight = 0.5))
  avg <- average_tdt_ttl(c(fits, list(spare = fits$exp)), w, newdata = nd, t_ref = 60,
                         n_draws = 200)
  expect_equal(avg$weights$model, w$model)
  expect_equal(sum(avg$weights$n_draws), 200L)
  expect_equal(avg$weights$n_draws[!avg$weights$kept], rep(0L, sum(!w$kept)))
})


test_that("the plots draw an average as they draw a single model", {
  avg <- average_tdt_ttl(fits, w_fits, newdata = nd, t_ref = 60, n_draws = 200)
  expect_s3_class(plot_tdt(avg), "ggplot")
  expect_s3_class(plot_ttl(avg), "ggplot")
  expect_s3_class(plot_ttl(avg, type = "curves"), "ggplot")
})


test_that("average_tdt_ttl refuses what cannot be averaged", {
  expect_error(average_tdt_ttl(fits, w_fits, newdata = nd), class = "rTFT_bad_t_ref")
  expect_error(average_tdt_ttl(unname(fits), w_fits, newdata = nd, t_ref = 60),
               class = "rTFT_bad_fits")
  expect_error(average_tdt_ttl(fits["exp"], w_fits, newdata = nd, t_ref = 60),
               class = "rTFT_bad_fits")
  expect_error(average_tdt_ttl(fits, c(exp = 0.5, zeroexp = 0.5), newdata = nd, t_ref = 60),
               class = "rTFT_bad_weights")
  expect_error(average_tdt_ttl(fits, w_fits, newdata = nd, t_ref = 60, n_draws = 0),
               class = "rTFT_bad_n_draws")
  expect_error(average_tdt_ttl(fits, w_fits, newdata = nd, t_ref = 60, prob = 1.5))
  # the same models on other data: one individual's rows relabelled
  other <- fits
  other$zeroexp$rtft$id[other$zeroexp$rtft$id == other$zeroexp$rtft$id[1]] <- -99
  expect_error(average_tdt_ttl(other, w_fits, newdata = nd, t_ref = 60),
               class = "rTFT_fits_mismatch")
})


test_that("a warning every model would give is given once", {
  given <- 0
  avg <- withCallingHandlers(
    average_tdt_ttl(fits, w_fits, newdata = nd, t_ref = 60, n_draws = 200,
                    levels = list(DO = 5)),
    rTFT_unused_levels = function(w) {
      given <<- given + 1
      invokeRestart("muffleWarning")
    })
  expect_equal(given, 1)
  expect_named(avg, c("sctmax", "tdt_curves", "ttl_surface", "weights"))
})


# ===========================================================================
# 4. dCTmax averaged over models
# ===========================================================================

test_that("average_dctmax returns the table of derive_dctmax, with the weights", {
  avg <- suppressMessages(average_dctmax(fits, w_fits, rate = c(0.1, 0.25), T_start = 20,
                                         newdata = nd, p = c(0.5, 0.9), n_draws = 400))
  expect_named(avg, c("dctmax", "weights"))
  one <- suppressMessages(derive_dctmax(fits$exp, rate = c(0.1, 0.25), T_start = 20,
                                        newdata = nd, p = c(0.5, 0.9), n_draws = 400))
  expect_named(avg$dctmax, names(one))
  expect_equal(avg$dctmax[c("mass", "rate", "T_start", "p")], one[c("mass", "rate", "T_start", "p")])
  expect_true(all(avg$dctmax$n_draws == 400))
  expect_equal(avg$weights$n_draws, c(300L, 100L))
  expect_true(all(avg$dctmax$lower <= avg$dctmax$median & avg$dctmax$median <= avg$dctmax$upper))
})


test_that("a model with all the weight gives back its own dCTmax", {
  alone <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0, 1))
  avg <- suppressMessages(average_dctmax(fits, alone, rate = c(0.1, 0.25), T_start = 20,
                                         newdata = nd, n_draws = 300))
  one <- suppressMessages(derive_dctmax(fits$zeroexp, rate = c(0.1, 0.25), T_start = 20,
                                        newdata = nd, n_draws = 300))
  expect_equal(avg$dctmax, one)
})


test_that("the averaged dCTmax is the summary of the models' draws pooled by weight", {
  one_row <- nd[1, , drop = FALSE]
  avg <- suppressMessages(average_dctmax(fits, w_fits, rate = 0.25, T_start = 20,
                                         newdata = one_row, n_draws = 400))
  d_exp  <- .dctmax_fit_draws(fits$exp,     0.25, 20, one_row, 0.5, NULL, 300, NULL)
  d_zero <- .dctmax_fit_draws(fits$zeroexp, 0.25, 20, one_row, 0.5, NULL, 100, NULL)
  pooled <- c(d_exp$rows[[1]][[1]][, 1], d_zero$rows[[1]][[1]][, 1])
  expect_length(pooled, 400)
  expect_equal(avg$dctmax$median, stats::median(pooled))
  expect_equal(avg$dctmax$lower, unname(stats::quantile(pooled, 0.025)))
  expect_equal(avg$dctmax$upper, unname(stats::quantile(pooled, 0.975)))
})


test_that("average_dctmax refuses what derive_dctmax and the averaging refuse", {
  expect_error(quietly(average_dctmax(fits, w_fits, rate = -0.1, T_start = 20, newdata = nd)),
               class = "rTFT_bad_rate")
  expect_error(quietly(average_dctmax(fits, w_fits, rate = 0.1, T_start = c(18, 20, 22),
                                      newdata = nd)),
               class = "rTFT_bad_T_start")
  expect_error(average_dctmax(unname(fits), w_fits, rate = 0.1, T_start = 20, newdata = nd),
               class = "rTFT_bad_fits")
  expect_error(average_dctmax(fits, w_fits, rate = 0.1, T_start = 20, newdata = nd, n_draws = 0),
               class = "rTFT_bad_n_draws")
})
