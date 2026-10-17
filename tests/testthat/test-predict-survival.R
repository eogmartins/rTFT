# ===========================================================================
# 1. Helpers
# ===========================================================================

# Two models fitted to the same data (40 individuals, log a depending on mass,
# memoryless baseline, reference temperature 20): the exponential failure rate and
# the zero-threshold exponential.
fits <- list(exp     = readRDS(test_path("fixtures", "hazard-fit-exp.rds")),
             zeroexp = readRDS(test_path("fixtures", "hazard-fit-zeroexp.rds")))

# Hand-made settings and draws, so the hazard along a record can be checked against
# the formula for any baseline without a fitted model.
settings <- function(m_form = "exp", baseline = "mem", costressors = list())
  list(m_form = m_form, baseline = baseline, Tr = 20, tau_ref = 60,
       k_sw = M_catalogue[[m_form]]$k_sw, costressors = costressors,
       Ta_range = c(22, 30))
d2 <- list(loga = c(-5, -4), logb = log(c(0.4, 0.5)), logc = c(-6, -6),
           logshape = log(c(2, 0.7)), logscale = log(c(40, 80)), gam = c(0.01, -0.02),
           cDO = c(0.3, 0.5), dTDO = c(0.02, 0.01))

# A half-hour record read every five minutes, warming then cooling, with the day
# changing after the fourth reading.
rec <- data.frame(time = seq(100, 130, by = 5),
                  temp = c(22, 23, 25, 26, 25, 24, 23),
                  day  = c(1, 1, 1, 1, 2, 2, 2))
flat <- data.frame(time = c(0, 10, 25, 60), temp = 25)

quietly <- function(expr) suppressMessages(expr)


# ===========================================================================
# 2. The record and the hazard accrued along it
# ===========================================================================

test_that(".check_record refuses a record the prediction cannot follow", {
  expect_no_error(.check_record(rec, c("time", "temp")))
  expect_error(.check_record(as.matrix(rec), c("time", "temp")), class = "rTFT_bad_record")
  expect_error(.check_record(rec["time"], c("time", "temp")), class = "rTFT_bad_record")
  expect_error(.check_record(rec[1, ], c("time", "temp")), class = "rTFT_bad_record")
  gap <- rec; gap$temp[3] <- NA
  expect_error(.check_record(gap, c("time", "temp")), class = "rTFT_bad_record")
  back <- rec; back$time[3] <- back$time[2]
  expect_error(.check_record(back, c("time", "temp")), class = "rTFT_bad_record")
  text <- rec; text$temp <- as.character(text$temp)
  expect_error(.check_record(text, c("time", "temp")), class = "rTFT_bad_record")
  # a missing value in a column the prediction does not use is no obstacle
  spare <- rec; spare$note <- NA
  expect_no_error(.check_record(spare, c("time", "temp")))
})


test_that(".record_intervals counts exposure from the first reading and takes midpoints", {
  iv <- .record_intervals(rec, "day")
  expect_equal(iv$tau_start, seq(0, 25, by = 5))
  expect_equal(iv$tau_end, seq(5, 30, by = 5))
  expect_equal(iv$Ta, c(22.5, 24, 25.5, 25.5, 24.5, 23.5))
  expect_equal(iv$day, c(1, 1, 1, 1.5, 2, 2))
  expect_named(.record_intervals(rec), c("tau_start", "tau_end", "Ta"))
})


test_that("with the memoryless baseline the hazard accrued is the rate times the time", {
  iv <- .record_intervals(rec)
  accrued <- .record_hazard(settings(), d2, iv)
  expect_equal(dim(accrued), c(2L, 6L))
  rate <- function(i) exp(d2$loga[i]) * exp(exp(d2$logb[i]) * (iv$Ta - 20))
  expect_equal(accrued[1, ], rate(1) * 5)
  expect_equal(accrued[2, ], rate(2) * 5)
})


test_that("at a constant temperature the hazard accrued adds up to M(T) H0(tau), whatever the baseline", {
  iv <- .record_intervals(flat)
  M  <- exp(d2$loga) * exp(exp(d2$logb) * 5)
  shape <- exp(d2$logshape)
  total <- function(baseline) rowSums(.record_hazard(settings(baseline = baseline), d2, iv))
  expect_equal(total("mem"), M * 60)
  # the Weibull scale is the reference time, 60, so H0(60) = 1
  expect_equal(total("weibull"), M)
  expect_equal(total("llogis"), M * log1p((60 / exp(d2$logscale))^shape))
  expect_equal(total("gompertz"), M * expm1(d2$gam * 60) / d2$gam)
  # and each interval takes the growth of H0 over it, not its length
  weibull <- .record_hazard(settings(baseline = "weibull"), d2, iv)
  expect_equal(weibull[1, ], M[1] * diff((c(0, 10, 25, 60) / 60)^shape[1]))
})


test_that("the zero-type rates are zero at and below the reference, and power is refused", {
  cold <- data.frame(time = c(0, 10, 20, 30), temp = c(16, 18, 22, 26))
  iv <- .record_intervals(cold)                        # midpoints 17, 20, 24
  zero <- .record_hazard(settings("zeroexp"), d2, iv)
  expect_equal(zero[, 1:2], matrix(0, 2, 2))
  expect_equal(zero[, 3], exp(d2$loga) * (exp(exp(d2$logb) * 4) - 1) * 10)
  zp <- .record_hazard(settings("zeropower"), d2, iv)
  expect_equal(zp[, 1:2], matrix(0, 2, 2))
  expect_true(all(zp[, 3] > 0))
  expect_error(.record_hazard(settings("power"), d2, iv), class = "rTFT_below_Tr")
  expect_no_error(.record_hazard(settings("power"), d2, .record_intervals(flat)))
  # the exponential rate has no such floor
  expect_true(all(.record_hazard(settings("exp"), d2, iv) > 0))
})


test_that("a co-stressor in the record is followed through time", {
  do_cs <- costressor("DO", column = "do", ref = 8, scale = 1, direction = -1,
                      interaction = TRUE)
  st <- settings(costressors = list(do_cs))
  steady <- transform(flat, do = 6)
  varying <- transform(flat, do = c(8, 8, 4, 4))       # midpoints 8, 6, 4
  at_level <- .record_hazard(st, d2, .record_intervals(flat), levels = c(DO = 6))
  from_record <- .record_hazard(st, d2, .record_intervals(steady, "do"), levels = c(DO = 0))
  expect_equal(from_record, at_level)                  # the record's value, not the level
  moving <- .record_hazard(st, d2, .record_intervals(varying, "do"), levels = c(DO = 0))
  expect_equal(moving[, 2], at_level[, 2])             # the interval at 6
  expect_true(all(moving[, 1] < at_level[, 1]))        # more oxygen, lower rate
  expect_true(all(moving[, 3] > at_level[, 3]))        # less oxygen, higher rate
  # each interval at the rate for its own oxygen
  z <- -(c(8, 6, 4) - 8)
  rate <- exp(d2$loga[1] + exp(d2$logb[1]) * 5 + d2$cDO[1] * z + d2$dTDO[1] * 5 * z)
  expect_equal(moving[1, ], rate * c(10, 15, 35))
})


test_that(".failure_temperatures takes the temperatures at which individuals failed", {
  # a ramp from 14 degrees in which the two failures came at 27.5 and 28.5: the
  # temperatures the individuals were exposed to start at 14, those at which the
  # data inform the failure rate at 27.5
  ramp <- list(data = data.frame(Ta = c(14.5, 20, 27.5, 14.5, 20, 28.5),
                                 failed = c(0, 0, 1, 0, 0, 1)))
  expect_equal(.failure_temperatures(ramp, list(Ta_range = c(14.5, 28.5))), c(27.5, 28.5))
  # with no failure at all there is nothing to take, and the whole range stands in
  none <- list(data = data.frame(Ta = c(14.5, 20), failed = c(0, 0)))
  expect_equal(.failure_temperatures(none, list(Ta_range = c(14.5, 20))), c(14.5, 20))
  # in the fixture individuals failed at both ends of the range
  expect_equal(.failure_temperatures(fits$exp, fits$exp$rtft), fits$exp$rtft$Ta_range)
})


test_that(".record_zone and .record_coverage say how much of a record, and of its risk, is outside what the data inform", {
  spans <- data.frame(time = c(0, 10, 30, 60, 100), temp = c(18, 20, 26, 32, 36))
  iv <- .record_intervals(spans)                         # midpoints 19, 23, 29, 34
  zone <- .record_zone(iv, c(22, 30))
  expect_equal(as.character(zone), c("below", "within", "within", "above"))
  expect_equal(levels(zone), c("below", "within", "above"))
  duration <- iv$tau_end - iv$tau_start                  # 10, 20, 30, 40
  # two draws: the hazard accrued in each interval
  accrued <- rbind(c(1, 1, 2, 4), c(0, 3, 3, 2))
  after <- c(FALSE, FALSE, TRUE, TRUE)                   # the last two come after tau 30
  cover <- .record_coverage(accrued, zone, duration, c(22, 30), after, 30)
  expect_equal(c(cover$below, cover$within, cover$above), c(0.1, 0.5, 0.4))
  expect_equal(c(cover$Ta_failed_min, cover$Ta_failed_max), c(22, 30))
  # shares of the risk: draw 1 gives 1/8, 3/8, 4/8; draw 2 gives 0, 6/8, 2/8
  expect_equal(cover$risk_below,  stats::median(c(1 / 8, 0)))
  expect_equal(cover$risk_within, stats::median(c(3 / 8, 6 / 8)))
  expect_equal(cover$risk_above,  stats::median(c(4 / 8, 2 / 8)))
  # in time: 70 of 100 by duration; of the risk, draw 1 gives 6/8 and draw 2 gives 5/8
  expect_equal(cover$tau_max, 30)
  expect_equal(cover$beyond, 0.7)
  expect_equal(cover$risk_beyond, stats::median(c(6 / 8, 5 / 8)))
  # a record no longer than the data has nothing beyond them
  short <- .record_coverage(accrued, zone, duration, c(22, 30), rep(FALSE, 4), 200)
  expect_equal(c(short$beyond, short$risk_beyond), c(0, 0))
  # most of the time below, most of the risk within
  warm_hours <- rbind(c(0.01, 5, 0, 0))
  spent <- .record_coverage(warm_hours, factor(c("below", "within", "above", "above"),
                                               levels = levels(zone)), c(90, 10, 0, 0), c(22, 30),
                            rep(FALSE, 4), 200)
  expect_equal(spent$below, 0.9)
  expect_gt(spent$risk_within, 0.99)
  # a draw with no hazard has no shares; with none at all, they are missing
  none <- .record_coverage(matrix(0, 2, 4), zone, duration, c(22, 30), after, 30)
  expect_true(all(is.na(c(none$risk_below, none$risk_within, none$risk_above,
                          none$risk_beyond))))
  expect_equal(none$within, 0.5)
  one <- .record_coverage(rbind(accrued[1, ], 0), zone, duration, c(22, 30), after, 30)
  expect_equal(one$risk_above, 0.5)
})


test_that(".risk_by_period gives each period's share, leaving out draws with no hazard", {
  accrued <- rbind(c(1, 1, 2, 4), c(0, 0, 0, 0), c(3, 0, 0, 1))
  risk <- .risk_by_period(accrued, period = c("a", "a", "b", "b"), by = "day", prob = 0.9)
  expect_named(risk, c("day", "median", "lower", "upper", "prob", "n_draws"))
  expect_equal(risk$day, c("a", "b"))
  expect_equal(risk$n_draws, c(2L, 2L))                # the draw with no hazard is left out
  expect_equal(risk$median, c(mean(c(0.25, 0.75)), mean(c(0.75, 0.25))))
  expect_equal(sum(risk$median), 1)
  none <- .risk_by_period(accrued[2, , drop = FALSE], c("a", "a", "b", "b"), "day", 0.9)
  expect_true(all(is.na(none$median)))
  expect_equal(none$n_draws, c(0L, 0L))
})


# ===========================================================================
# 3. Survival along a record, from a fitted model
# ===========================================================================

nd <- data.frame(mass = c(8, 12))

test_that("the mean of survival is the mean over the draws, and can differ from the median", {
  # three draws over two intervals, two of which keep nearly everyone alive and one
  # of which does not: survival at the end is 0.98, 0.97 and 0.2
  accrued <- -log(rbind(c(0.99, 0.98 / 0.99), c(0.98, 0.97 / 0.98), c(0.5, 0.2 / 0.5)))
  draws <- list(newdata = data.frame(row.names = 1L), time = c(0, 10, 20), rows = list(accrued),
                by = NULL, zone = factor(c("within", "within"), levels = c("below", "within", "above")),
                duration = c(10, 10), Ta_failed = c(22, 30), beyond_tau = c(FALSE, FALSE),
                tau_max = 60)
  out <- .summarise_record(draws, 0.95)
  expect_equal(out$survival$mean, c(1, mean(c(0.99, 0.98, 0.5)), mean(c(0.98, 0.97, 0.2))))
  expect_equal(out$survival$median, c(1, 0.98, 0.97))
  expect_equal(out$mortality$mean, 1 - mean(c(0.98, 0.97, 0.2)))
  expect_equal(out$mortality$median, 0.03)
})


test_that("predict_survival returns survival through time and mortality at the end", {
  pred <- quietly(predict_survival(fits$exp, rec, newdata = nd, n_draws = 300))
  expect_named(pred, c("survival", "mortality", "coverage", "record"))
  expect_named(pred$coverage, c("mass", "Ta_failed_min", "Ta_failed_max", "below", "within",
                                "above", "risk_below", "risk_within", "risk_above",
                                "tau_max", "beyond", "risk_beyond"))
  expect_equal(pred$record, tibble::tibble(time = rec$time, temp = rec$temp))   # not `day`
  expect_named(pred$survival, c("mass", "time", "mean", "median", "lower", "upper", "prob",
                                "n_draws", "n_unreached"))
  expect_named(pred$mortality, names(pred$survival))
  expect_equal(nrow(pred$survival), 2L * nrow(rec))
  expect_equal(pred$survival$time, rep(rec$time, 2))
  first <- pred$survival[pred$survival$time == 100, ]
  expect_equal(first[c("median", "lower", "upper")],
               tibble::tibble(median = c(1, 1), lower = c(1, 1), upper = c(1, 1)))
  for (m in nd$mass)
    expect_true(all(diff(pred$survival$median[pred$survival$mass == m]) < 0))
  expect_true(all(pred$survival$median > 0 & pred$survival$median <= 1))
  expect_equal(pred$mortality$time, c(130, 130))
  expect_equal(pred$mortality$median,
               1 - pred$survival$median[pred$survival$time == 130])
  # the mean of mortality is one minus the mean of survival at the last time
  expect_equal(pred$mortality$mean, 1 - pred$survival$mean[pred$survival$time == 130])
  expect_equal(pred$survival$mean[pred$survival$time == 100], c(1, 1))
  expect_true(all(pred$survival$mean >= pred$survival$lower &
                    pred$survival$mean <= pred$survival$upper))
  # heavier individuals have the higher rate in these data
  expect_gt(pred$mortality$median[2], pred$mortality$median[1])
})


test_that("the prediction is the model's rate integrated along the record, draw by draw", {
  one <- nd[1, , drop = FALSE]
  dr <- .fit_draws(fits$exp, c("loga", "logb"), one, draw_ids = .thin_ids(600, 300))
  a <- exp(dr$loga[, 1]); b <- exp(dr$logb[, 1])
  T0 <- rec$temp[-7]; T1 <- rec$temp[-1]                 # the two ends of each interval
  # between two readings the temperature is linear, so the integral of
  # a exp(b (T - 20)) over an interval of 5 minutes has a closed form
  exact <- vapply(seq_along(a), \(i)
    sum(a[i] * 5 * (exp(b[i] * (T1 - 20)) - exp(b[i] * (T0 - 20))) / (b[i] * (T1 - T0))),
    numeric(1))
  pred <- quietly(predict_survival(fits$exp, rec, newdata = one, n_draws = 300))
  expect_equal(pred$mortality$median, stats::median(1 - exp(-exact)), tolerance = 1e-3)
  expect_equal(pred$mortality$lower, unname(stats::quantile(1 - exp(-exact), 0.025)),
               tolerance = 1e-3)
  # one step between readings takes the rate at the midpoint of each interval,
  # which falls short of the integral where the rate curves upward
  mid <- (T0 + T1) / 2
  by_midpoint <- vapply(seq_along(a), \(i) sum(a[i] * exp(b[i] * (mid - 20)) * 5), numeric(1))
  coarse <- quietly(predict_survival(fits$exp, rec, newdata = one, n_draws = 300, grid_dT = Inf))
  expect_equal(coarse$mortality$median, stats::median(1 - exp(-by_midpoint)))
  expect_lt(coarse$mortality$median, pred$mortality$median)
  expect_error(predict_survival(fits$exp, rec, newdata = one, grid_dT = 0),
               class = "rTFT_bad_grid_dT")
})


test_that(".refine_intervals divides the record where temperature changes, linearly between readings", {
  two <- data.frame(time = c(100, 110, 130), temp = c(20, 20.3, 20.3), do = c(8, 5, 5))
  fine <- .refine_intervals(two, "do", grid_dT = 0.1)
  # the first interval changes by 0.3 degrees: three parts; the second not at all
  expect_equal(fine$parent, c(1, 1, 1, 2))
  expect_equal(fine$tau_start, c(0, 10 / 3, 20 / 3, 10))
  expect_equal(fine$tau_end, c(10 / 3, 20 / 3, 10, 30))
  expect_equal(fine$Ta, c(20.05, 20.15, 20.25, 20.3))
  expect_equal(fine$do, c(7.5, 6.5, 5.5, 5))             # the co-stressor follows along
  # as given: the layout of .record_intervals()
  whole <- .refine_intervals(two, "do", grid_dT = Inf)
  plain <- .record_intervals(two, "do")
  expect_equal(whole[names(plain)], plain)
  expect_equal(whole$parent, 1:2)
  # a cooling interval is divided as a warming one is
  expect_equal(.refine_intervals(data.frame(time = c(0, 10), temp = c(21, 20.6)))$parent, rep(1, 4))
})


test_that("with the memoryless baseline it does not matter when the record begins", {
  later <- transform(rec, time = time + 5000)
  a <- quietly(predict_survival(fits$exp, rec, newdata = nd, n_draws = 200))
  b <- quietly(predict_survival(fits$exp, later, newdata = nd, n_draws = 200))
  expect_equal(a$mortality$median, b$mortality$median)
  expect_equal(a$survival$median, b$survival$median)
})


test_that("by gives each period's share of the risk", {
  pred <- quietly(predict_survival(fits$exp, rec, newdata = nd, by = "day", n_draws = 200))
  expect_named(pred, c("survival", "mortality", "risk", "coverage", "record"))
  expect_equal(pred$risk$day, c(1, 2, 1, 2))
  expect_equal(pred$risk$n_draws, rep(200L, 4))
  # intervals belong to the period of their first reading: four in day 1, two in day 2
  expect_true(all(pred$risk$median[pred$risk$day == 1] > 0.5))
  shares <- tapply(pred$risk$median, pred$risk$mass, sum)
  expect_equal(as.vector(shares), c(1, 1), tolerance = 0.02)
  expect_error(predict_survival(fits$exp, rec, newdata = nd, by = "week"),
               class = "rTFT_bad_by")
  expect_error(predict_survival(fits$exp, rec, newdata = nd, by = c("day", "time")),
               class = "rTFT_bad_by")
})


test_that("a record longer than the exposures in the data is said to be", {
  # the longest exposure in the fixture is 45; `rec` spans 30, a record of 180 does not
  expect_equal(quietly(predict_survival(fits$exp, rec, newdata = nd, n_draws = 100))$
                 coverage$beyond, c(0, 0))
  long_rec <- data.frame(time = seq(0, 180, by = 15), temp = 25)
  expect_message(pred <- predict_survival(fits$exp, long_rec, newdata = nd, n_draws = 100),
                 "after the longest exposure in the data \\(45\\)", class = "rTFT_record_notes")
  expect_equal(pred$coverage$tau_max, c(45, 45))
  expect_equal(pred$coverage$beyond, c(0.75, 0.75))        # from 45 to 180, of 180
  # at a constant temperature and with a memoryless baseline the risk accrues evenly
  expect_equal(pred$coverage$risk_beyond, c(0.75, 0.75))
  expect_message(predict_survival(fits$exp, long_rec, newdata = nd, n_draws = 100),
                 "75% of the predicted risk is accrued after it", class = "rTFT_record_notes")
})


test_that("a record outside the fitted temperatures is said to be", {
  # the fits cover 20.6 to 30.6; `rec` lies within, a cooler record does not
  expect_no_message(predict_survival(fits$exp, rec, newdata = nd, n_draws = 100))
  cool <- transform(rec, temp = temp - 6)
  expect_message(pred <- predict_survival(fits$exp, cool, newdata = nd, n_draws = 100),
                 "outside the temperatures", class = "rTFT_record_notes")
  expect_equal(pred$coverage$below, c(1, 1))               # a row per row of newdata
  expect_equal(pred$coverage$mass, nd$mass)
  expect_equal(pred$coverage$within + pred$coverage$below + pred$coverage$above, c(1, 1))
  expect_equal(pred$coverage$risk_below, c(1, 1))          # and all of the risk with it
  expect_message(predict_survival(fits$exp, cool, newdata = nd, n_draws = 100),
                 "100% of the predicted risk", class = "rTFT_record_notes")
  # a record that only dips below the fitted temperatures takes little risk from there
  dips <- transform(rec, temp = ifelse(seq_along(temp) <= 2, 19, temp + 3))
  mostly <- quietly(predict_survival(fits$exp, dips, newdata = nd, n_draws = 100))
  expect_gt(mostly$coverage$below[1], 0.1)
  expect_lt(mostly$coverage$risk_below[1], mostly$coverage$below[1])
  expect_equal(with(mostly$coverage, risk_below + risk_within + risk_above), c(1, 1),
               tolerance = 0.02)
  # below the reference the zero-type rate gives no risk at all, and says so
  zero <- quietly(predict_survival(fits$zeroexp, cool, newdata = nd, by = "day",
                                   n_draws = 100))
  expect_true(all(zero$mortality$median < 1e-3))
  cold <- transform(rec, temp = temp - 10)
  none <- quietly(predict_survival(fits$zeroexp, cold, newdata = nd, by = "day",
                                   n_draws = 100))
  expect_equal(none$mortality$median, c(0, 0))
  expect_true(all(is.na(none$risk$median)))
  expect_equal(none$risk$n_draws, rep(0L, 4))
})


test_that("predict_survival refuses bad arguments", {
  expect_error(predict_survival(fits$exp, rec["time"], newdata = nd), class = "rTFT_bad_record")
  expect_error(predict_survival(fits$exp, rec, newdata = nd, prob = 2), class = "rTFT_bad_prob")
  expect_error(predict_survival(fits$exp, rec), class = "rTFT_missing_covariates")   # mass
  expect_error(predict_survival(list(), rec), class = "rTFT_not_brmsfit")
  expect_warning(quietly(predict_survival(fits$exp, rec, newdata = nd, levels = c(DO = 5),
                                          n_draws = 50)),
                 class = "rTFT_unused_levels")
})


# ===========================================================================
# 4. Survival along a record, averaged over models
# ===========================================================================

w_fits <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0.75, 0.25))

test_that("average_survival returns the tables of predict_survival, with the weights", {
  avg <- quietly(average_survival(fits, w_fits, rec, newdata = nd, by = "day", n_draws = 400))
  one <- quietly(predict_survival(fits$exp, rec, newdata = nd, by = "day", n_draws = 400))
  expect_named(avg, c("survival", "mortality", "risk", "coverage", "record", "weights"))
  expect_named(avg$survival, names(one$survival))
  expect_named(avg$risk, names(one$risk))
  expect_named(avg$coverage, names(one$coverage))
  expect_equal(avg$coverage[c("mass", "below", "within", "above")],
               one$coverage[c("mass", "below", "within", "above")])
  expect_equal(avg$weights$n_draws, c(300L, 100L))
  expect_true(all(avg$survival$n_draws == 400))
})


test_that("a model with all the weight gives back its own prediction", {
  alone <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0, 1))
  avg <- quietly(average_survival(fits, alone, rec, newdata = nd, by = "day", n_draws = 300))
  one <- quietly(predict_survival(fits$zeroexp, rec, newdata = nd, by = "day", n_draws = 300))
  expect_equal(avg[c("survival", "mortality", "risk", "coverage", "record")], one)
})


test_that("the average is the summary of the models' draws pooled by weight", {
  one <- nd[1, , drop = FALSE]
  avg <- quietly(average_survival(fits, w_fits, rec, newdata = one, n_draws = 400))
  d_exp  <- .record_draws(fits$exp,     rec, one, NULL, NULL, 300)$rows[[1]]
  d_zero <- .record_draws(fits$zeroexp, rec, one, NULL, NULL, 100)$rows[[1]]
  pooled <- 1 - exp(-c(rowSums(d_exp), rowSums(d_zero)))
  expect_equal(avg$mortality$median, stats::median(pooled))
  expect_equal(avg$mortality$upper, unname(stats::quantile(pooled, 0.975)))
})


test_that("average_survival refuses what cannot be averaged, and warns once", {
  expect_error(average_survival(unname(fits), w_fits, rec, newdata = nd), class = "rTFT_bad_fits")
  expect_error(average_survival(fits, c(exp = 1), rec, newdata = nd), class = "rTFT_bad_weights")
  expect_error(average_survival(fits, w_fits, rec["time"], newdata = nd, n_draws = 100),
               class = "rTFT_bad_record")
  given <- 0
  quietly(withCallingHandlers(
    average_survival(fits, w_fits, rec, newdata = nd, levels = c(DO = 5), n_draws = 100),
    rTFT_unused_levels = function(w) {
      given <<- given + 1
      invokeRestart("muffleWarning")
    }))
  expect_equal(given, 1)
})


# ===========================================================================
# 5. A fitted model of static and ramping trials together
# ===========================================================================

# A small fitted model saved as a fixture (made by fixtures/make-hazard-fit-mixed.R):
# three static trials and two ramps of ten individuals each, with the exponential
# failure rate, the memoryless baseline, no covariate and a reference of 20 degC.
fit_mixed <- readRDS(test_path("fixtures", "hazard-fit-mixed.rds"))

test_that("survival along a record is predicted from a fit of both kinds of assay", {
  # An hour at 25 degC: with the memoryless baseline the hazard accrued is the rate
  # at that temperature times the hour, draw by draw.
  cf <- as.data.frame(posterior::as_draws_df(fit_mixed))[.thin_ids(600, 300), ]
  rate    <- exp(cf$b_loga_Intercept + exp(cf$b_logb_Intercept) * (25 - 20))
  by_hand <- 1 - exp(-rate * 60)
  pred <- quietly(predict_survival(fit_mixed, flat, n_draws = 300))
  expect_named(pred, c("survival", "mortality", "coverage", "record"))
  expect_equal(pred$mortality$time, 60)
  expect_equal(pred$mortality$median, stats::median(by_hand))
  expect_equal(pred$mortality$mean, mean(by_hand))
  expect_equal(pred$mortality$lower, unname(stats::quantile(by_hand, 0.025)))
  # 25 degC lies within the temperatures at which individuals failed, in the static
  # trials and under the ramps alike
  expect_equal(pred$coverage$within, 1)
})
