# ===========================================================================
# 1. Baseline cumulative hazard and its inverse (H0_of / H0inv_of)
# ===========================================================================

# H0_of() and H0inv_of() return closed-form functions of (tau/x, shape, scale). The
# expected values are written out by hand from the four formulas, and the inverse is
# checked by round-tripping H0inv(H0(tau)) back to tau, so the test pins the algebra
# rather than re-running the functions' own expressions.

test_that("H0_of reproduces each baseline's cumulative hazard", {
  # mem: H0(tau) = tau
  expect_equal(H0_of("mem")(c(10, 30, 60), shape = 1, scale = 1), c(10, 30, 60))
  # weibull: (tau / scale)^shape = (30 / 60)^2 = 0.25
  expect_equal(H0_of("weibull")(30, shape = 2, scale = 60), 0.25)
  # llogis: log(1 + (tau / scale)^shape) = log(1 + 0.25) = log(1.25)
  expect_equal(H0_of("llogis")(30, shape = 2, scale = 60), log(1.25))
  # gompertz: (exp(g tau) - 1) / g = (exp(1.5) - 1) / 0.05
  expect_equal(H0_of("gompertz")(30, shape = 0.05, scale = 1),
               (exp(1.5) - 1) / 0.05)
})


test_that("H0inv_of inverts H0_of for every baseline", {
  for (b in c("mem", "weibull", "llogis", "gompertz")) {
    H0  <- H0_of(b)
    H0i <- H0inv_of(b)
    shape <- if (b == "gompertz") 0.05 else 2
    scale <- 60
    tau <- 30
    expect_equal(H0i(H0(tau, shape, scale), shape, scale), tau,
                 info = paste("round-trip failed for", b))
  }
})


test_that("the Gompertz forms pass smoothly through gamma = 0 (the memoryless limit)", {
  # At a rate constant far below the tolerance, H0 and its inverse must agree with
  # the memoryless case (H0(tau) = tau, H0^-1(x) = x) rather than return NaN from 0/0.
  expect_equal(H0_of("gompertz")(30, shape = 1e-12, scale = 1), 30, tolerance = 1e-6)
  expect_equal(H0inv_of("gompertz")(30, shape = 1e-12, scale = 1), 30, tolerance = 1e-6)
  # exactly zero is the singular point, and it must still be finite, not NaN
  expect_equal(H0_of("gompertz")(30, shape = 0, scale = 1), 30)
  expect_equal(H0inv_of("gompertz")(30, shape = 0, scale = 1), 30)
})


test_that("a negative Gompertz rate constant is a bounded hazard with an honest cure fraction", {
  # gamma = -0.05: H0(tau) rises to the bound -1/gamma = 20 and no higher.
  big <- H0_of("gompertz")(1e6, shape = -0.05, scale = 1)
  expect_lte(big, 20)                                 # never exceeds the bound
  expect_gt(big, 19.9)                                # but approaches it
  # a target below the bound is reached at a finite time
  expect_true(is.finite(H0inv_of("gompertz")(10, shape = -0.05, scale = 1)))
  # a target above the bound is reached at no finite time: NaN, not a large number
  expect_true(is.nan(H0inv_of("gompertz")(25, shape = -0.05, scale = 1)))
  # a target exactly at the bound (gamma = -0.5, bound 2) is reached only after
  # infinite time: Inf, not NaN
  expect_equal(H0inv_of("gompertz")(2, shape = -0.5, scale = 1), Inf)
})


test_that("H0 returns one value per draw even when the time is a single number", {
  # The memoryless baseline, H0(tau) = tau, ignores shape and scale, so a single tau
  # with a per-draw shape and scale must still come back as one value per draw.
  out <- H0_of("mem")(60, shape = rep(1, 5), scale = rep(1, 5))
  expect_length(out, 5)
  expect_true(all(out == 60))
  # the same for the inverse
  expect_length(H0inv_of("mem")(3, shape = rep(1, 4), scale = rep(1, 4)), 4)
})


test_that("a matrix time argument keeps its dimensions", {
  # A draws-by-times matrix must come back the same shape (the surface path relies on
  # this), so the short-argument growth must not reshape it.
  tau <- matrix(c(10, 20, 30, 40), nrow = 2)
  out <- H0_of("weibull")(tau, shape = 2, scale = 60)
  expect_equal(dim(out), c(2L, 2L))
  expect_equal(out, (tau / 60)^2)
})


test_that("an unknown baseline errors with the shared condition class", {
  expect_error(H0_of("nope"), class = "rTFT_unknown_baseline")
  expect_error(H0inv_of("nope"), class = "rTFT_unknown_baseline")
  expect_error(H0_of(factor("weibull")), class = "rTFT_bad_baseline")
})


# ===========================================================================
# 2. The failure rate for each posterior draw (.M_draws)
# ===========================================================================

# Expected rates are written out by hand from the M(T) formulas, so the tests pin
# the arithmetic rather than re-evaluating the package's own strings.

test_that(".M_draws gives the exponential rate for every draw and temperature", {
  draws <- list(loga = c(-5, -4), logb = log(c(0.5, 0.2)))
  Ta <- c(30, 32, 35)
  out <- .M_draws("exp", draws, Ta = Ta, Tr = 30)
  expect_equal(dim(out), c(2L, 3L))                      # draws x temperatures
  # M = exp(loga + b * (Ta - Tr)), row by row
  expect_equal(out, rbind(exp(-5 + 0.5 * (Ta - 30)), exp(-4 + 0.2 * (Ta - 30))))
})


test_that("a constant rate is the same at every temperature, one value per draw", {
  out <- .M_draws("const", list(loga = c(-5, -6, -7)), Ta = c(20, 30), Tr = 20)
  expect_equal(out, matrix(exp(c(-5, -6, -7)), 3, 2))
})


test_that("a threshold rate uses the catalogue steepness, or the one it was fitted with", {
  # exp_threshold just above Tc: log M = loga + b (Ta - Tc) - log(1 + exp(-k (Ta - Tc)))
  draws <- list(loga = -5, logb = log(0.5), Tc = 34)
  hand  <- function(k) exp(-5 + 0.5 * 0.1 - log1p(exp(-k * 0.1)))
  expect_equal(.M_draws("exp_threshold", draws, Ta = 34.1, Tr = 30),
               matrix(hand(50)))                         # catalogue k_sw = 50
  expect_equal(.M_draws("exp_threshold", draws, Ta = 34.1, Tr = 30, k_sw = 10),
               matrix(hand(10)))
})


test_that("a co-stressor enters at its chosen level, standardised as in the fit", {
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  draws <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = 0.03)
  # at 6 mg/L, z = -(6 - 8.5) / 1 = 2.5, so
  # log M = loga + b (Ta - Tr) + cDO z + dTDO (Ta - Tr) z
  out <- .M_draws("exp", draws, Ta = 32, Tr = 30,
                 costressors = list(do_cs), levels = c(DO = 6))
  expect_equal(out, matrix(exp(-5 + 0.5 * 2 + 0.2 * 2.5 + 0.03 * 2 * 2.5)))
})


test_that("a co-stressor with a transform enters through that transform", {
  # pH harmful either side of 7.5: z = ((pH - 7.5) / 0.5)^2, so pH 6.5 gives z = 4
  ph_cs <- costressor("pH", column = "ph", ref = 7.5, scale = 0.5, direction = 1,
                      interaction = FALSE, transform = function(u) u^2)
  draws <- list(loga = -5, logb = log(0.5), cpH = 0.1)
  out <- .M_draws("exp", draws, Ta = 32, Tr = 30,
                  costressors = list(ph_cs), levels = c(pH = 6.5))
  expect_equal(out, matrix(exp(-5 + 0.5 * 2 + 0.1 * 4)))
})


test_that("draws may be a data frame", {
  draws <- data.frame(loga = c(-5, -4), logb = log(c(0.5, 0.2)), other = 1:2)
  expect_equal(.M_draws("exp", draws, Ta = 30, Tr = 30), matrix(exp(c(-5, -4))))
})


test_that(".M_draws refuses missing or mismatched draws", {
  expect_error(.M_draws("exp", list(loga = -5), Ta = 30, Tr = 30),
               class = "rTFT_missing_draws")             # no draws for logb
  expect_error(.M_draws("exp", list(loga = c(-5, -4), logb = 0), Ta = 30, Tr = 30),
               class = "rTFT_bad_draws")
})


test_that(".M_draws needs a level for each co-stressor", {
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  draws <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = 0.03)
  expect_error(.M_draws("exp", draws, Ta = 32, Tr = 30, costressors = list(do_cs)),
               class = "rTFT_bad_levels")
  expect_error(.M_draws("exp", draws, Ta = 32, Tr = 30, costressors = list(do_cs),
                       levels = c(pH = 7)), class = "rTFT_bad_levels")
})


test_that(".M_draws warns about levels the model has no co-stressor for", {
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  draws <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = 0.03)
  # an extra level beside the one needed: warned about, and the rate is still the
  # one at DO = 6
  extra <- function() .M_draws("exp", draws, Ta = 32, Tr = 30,
                                costressors = list(do_cs), levels = c(DO = 6, pH = 7))
  expect_warning(extra(), class = "rTFT_unused_levels")
  expect_equal(suppressWarnings(extra()),
               .M_draws("exp", draws, Ta = 32, Tr = 30, costressors = list(do_cs),
                        levels = c(DO = 6)))
  # a level given to a model with no co-stressors
  expect_warning(.M_draws("exp", draws, Ta = 32, Tr = 30, levels = c(DO = 6)),
                 class = "rTFT_unused_levels")
})


test_that(".M_draws checks its function name, temperatures and reference", {
  draws <- list(loga = -5, logb = log(0.5))
  expect_error(.M_draws("nope", draws, Ta = 30, Tr = 30), class = "rTFT_unknown_m_form")
  expect_error(.M_draws("exp", draws, Ta = c(30, NA), Tr = 30), class = "rTFT_bad_Ta")
  expect_error(.M_draws("exp", draws, Ta = 30, Tr = c(20, 30)), class = "rTFT_bad_Tr")
})


test_that("the functions of (Ta - Tr) refuse a temperature grid below Tr", {
  draws <- list(loga = -5, logb = 0, logc = -8)
  for (m in c("zeroexp", "zeropower", "power"))
    expect_error(.M_draws(m, draws, Ta = c(29, 31), Tr = 30), class = "rTFT_Ta_below_Tr",
                 info = m)
  # at Ta = Tr a zero function has a rate of exactly zero, its definition
  expect_equal(.M_draws("zeropower", draws, Ta = 30, Tr = 30), matrix(0))
  # the other functions are defined below Tr
  expect_true(all(is.finite(.M_draws("exp", draws, Ta = c(25, 35), Tr = 30))))
})


test_that(".log_M_eval evaluates each draw at its own temperatures", {
  # Row i of the temperature matrix belongs to draw i. Two draws, each at two
  # temperatures of its own, so a mix-up of draws and temperatures would show.
  setup  <- .rate_setup("exp", list(loga = c(-5, -4), logb = log(c(0.5, 0.2))), Tr = 30)
  Ta_mat <- rbind(c(30, 32),                             # draw 1
                  c(34, 36))                             # draw 2
  # log M = loga + b (Ta - Tr), row by row
  expect_equal(.log_M_eval(setup, Ta_mat),
               rbind(-5 + 0.5 * c(0, 2),
                     -4 + 0.2 * c(4, 6)))
})


test_that(".M_draws reports the rate the model formula fits", {
  # Evaluate the right-hand side of the fitted model's formula at fixed parameter
  # values, and compare it with the rate .M_draws reports for the same settings.
  # Under the memoryless baseline the formula is log M + log_exposure, so the two
  # must agree exactly: a threshold function with a custom steepness and a
  # co-stressor, the cases where fitting and reporting could most easily part.
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  pars  <- list(loga = -5, logb = log(0.5), Tc = 34, cDO = 0.2, dTDO = 0.03)
  spec  <- build_hazard_spec("exp_threshold", k_sw = 20, costressors = list(do_cs))
  Ta    <- c(33, 34.5, 36)
  log_exposure <- log(c(1, 2, 0.5))
  eta <- eval(spec$formula$formula[[3]],
              c(pars, list(Ta = Ta, Tr = 30, z_DO = 2.5, log_exposure = log_exposure)),
              enclos = baseenv())
  M <- .M_draws("exp_threshold", pars, Ta = Ta, Tr = 30, k_sw = 20,
               costressors = list(do_cs), levels = c(DO = 6))   # z_DO = 2.5 at 6 mg/L
  expect_equal(eta, as.vector(log(M)) + log_exposure)
})


# ===========================================================================
# 3. Posterior summaries (.draw_summary)
# ===========================================================================

test_that(".draw_summary gives the median and 95% equal-tailed interval", {
  x   <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3)
  out <- .draw_summary(x)
  expect_equal(c(out$median, out$lower, out$upper),
               unname(stats::quantile(x, c(0.5, 0.025, 0.975))))
  expect_equal(c(out$n_draws, out$n_unreached), c(10L, 0L))
})


test_that(".draw_summary takes the interval's coverage from prob", {
  x <- 1:101                                             # quantiles fall on draws
  expect_equal(c(.draw_summary(x, prob = 0.5)$lower, .draw_summary(x, prob = 0.5)$upper),
               c(26, 76))                                # the middle 50%
  expect_equal(c(.draw_summary(x, prob = 0.9)$lower, .draw_summary(x, prob = 0.9)$upper),
               c(6, 96))                                 # the middle 90%
  expect_equal(.draw_summary(x, prob = 0.8)$prob, 0.8)   # the coverage is recorded
  for (bad in list(0, 1, 1.5, c(0.5, 0.9), NA_real_, "0.9"))
    expect_error(.draw_summary(x, prob = bad), class = "rTFT_bad_prob")
})


test_that("draws beyond reach are kept, not dropped", {
  # 40 of 100 draws never reach the quantity. Kept, the median is 50.5; dropped,
  # it would be 30.5, the median of only the draws where it was reached.
  x   <- c(1:60, rep(Inf, 40))
  out <- .draw_summary(x)
  expect_equal(out$median, 50.5)
  expect_equal(out$upper, Inf)                           # no upper bound
  expect_equal(out$n_unreached, 40L)
})


test_that("the bounds and median turn infinite as the share beyond grows", {
  few  <- .draw_summary(c(1:999, Inf))                   # 0.1% beyond
  most <- .draw_summary(c(1:40, rep(Inf, 60)))           # 60% beyond
  expect_true(is.finite(few$upper))                      # under 2.5%: still bounded
  expect_equal(most$median, Inf)                         # half or more: no median
})


test_that("draws beyond on both sides do not produce NaN", {
  # Half unattainable below and half above, nothing finite between: interpolating
  # would give Inf - Inf; the nearest draw is taken instead.
  out <- .draw_summary(c(rep(-Inf, 50), rep(Inf, 50)))
  expect_false(anyNA(c(out$median, out$lower, out$upper)))
  expect_equal(c(out$lower, out$upper), c(-Inf, Inf))
})


test_that("a missing draw is a fault, not something to drop", {
  expect_error(.draw_summary(c(1, 2, NA)),  class = "rTFT_missing_draw_values")
  expect_error(.draw_summary(c(1, 2, NaN)), class = "rTFT_missing_draw_values")
  expect_error(.draw_summary(numeric(0)),   class = "rTFT_bad_summary_draws")
  expect_error(.draw_summary("a"),          class = "rTFT_bad_summary_draws")
})


# ===========================================================================
# 4. The static critical thermal maximum (.sctmax_draws)
# ===========================================================================

# Where a failure-rate function can be solved for T by hand, the closed form is
# the oracle: log M(T) = log(target) = L solved for T. The solver itself never
# uses these; it bisects the shared expression.

test_that("sCTmax matches the closed form for the functions that have one", {
  L <- log(c(0.01, 0.05, 0.2))                           # three targets
  tgt <- exp(L)
  d <- list(loga = c(-5, -4.5, -6), logb = log(c(0.5, 0.3, 0.8)), logc = rep(log(0.001), 3))
  a <- d$loga; b <- exp(d$logb)

  # exp: loga + b (T - Tr) = L
  expect_equal(.sctmax_draws("exp", d, tgt, Tr = 30), 30 + (L - a) / b, tolerance = 1e-7)
  # zeroexp: loga + log(exp(b (T - Tr)) - 1) = L
  expect_equal(.sctmax_draws("zeroexp", d, tgt, Tr = 30),
               30 + log1p(exp(L - a)) / b, tolerance = 1e-7)
  # zeropower: loga + b log(T - Tr) = L
  expect_equal(.sctmax_draws("zeropower", d, tgt, Tr = 30),
               30 + exp((L - a) / b), tolerance = 1e-7)
  # power: log(c + exp(loga) (T - Tr)^b) = L
  expect_equal(.sctmax_draws("power", d, tgt, Tr = 30),
               30 + ((tgt - 0.001) / exp(a))^(1 / b), tolerance = 1e-7)
  # const_threshold: loga - log(1 + exp(-k (T - Tc))) = L
  dt <- list(loga = c(-1, -1.5, -0.5), Tc = c(34, 35, 33))
  expect_equal(.sctmax_draws("const_threshold", dt, tgt, Tr = 30),
               dt$Tc - log(expm1(dt$loga - L)) / 50, tolerance = 1e-7)
})


test_that("sCTmax solves the functions without a closed form", {
  # For these the check is that the rate at the answer equals the target.
  tgt <- c(0.01, 0.05, 0.2)
  d <- list(loga = c(-5, -4.5, -6), logb = log(c(0.5, 0.3, 0.8)), Tc = c(34, 35, 33))
  for (m in c("exp_threshold", "zeroexp_threshold", "zeropower_threshold")) {
    ct <- .sctmax_draws(m, d, tgt, Tr = 30)
    M_at <- vapply(1:3, \(i) .M_draws(m, lapply(d, `[`, i), Ta = ct[i], Tr = 30)[1, 1],
                   numeric(1))
    expect_equal(M_at, tgt, tolerance = 1e-7, info = m)
  }
  # the zero-threshold answers lie on the rising branch, above Tc
  expect_true(all(.sctmax_draws("zeroexp_threshold", d, tgt, Tr = 30) > d$Tc))
})


test_that("a custom steepness and a co-stressor move sCTmax as the formula says", {
  # exp with oxygen: loga + b (T - Tr) + cDO z + dTDO (T - Tr) z = L, z = 2.5 at 6 mg/L
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  d <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = 0.03)
  L <- log(0.05)
  expect_equal(.sctmax_draws("exp", d, 0.05, Tr = 30, costressors = list(do_cs),
                             levels = c(DO = 6)),
               30 + (L + 5 - 0.2 * 2.5) / (0.5 + 0.03 * 2.5), tolerance = 1e-7)
  # const_threshold with k_sw = 10 rather than 50
  dt <- list(loga = -1, Tc = 34)
  expect_equal(.sctmax_draws("const_threshold", dt, 0.05, Tr = 30, k_sw = 10),
               34 - log(expm1(-1 - L)) / 10, tolerance = 1e-7)
})


test_that("a rate that never crosses the target gives +Inf or -Inf, not NA", {
  # const: below the target everywhere -> +Inf; above it everywhere -> -Inf
  expect_equal(.sctmax_draws("const", list(loga = log(c(0.01, 0.5))), c(0.1, 0.1), Tr = 30),
               c(Inf, -Inf))
  # const_threshold whose ceiling a = exp(loga) lies below the target
  expect_equal(.sctmax_draws("const_threshold", list(loga = log(0.05), Tc = 34), 0.1, Tr = 30),
               Inf)
  # power whose floor c lies above the target
  expect_equal(.sctmax_draws("power", list(loga = -5, logb = 0, logc = log(0.5)), 0.1, Tr = 30),
               -Inf)
  # a crossing beyond the search span counts as never reached
  expect_equal(.sctmax_draws("exp", list(loga = -20, logb = log(1e-4)), 0.1, Tr = 30), Inf)
})


test_that("sCTmax refuses a co-stressor level at which the rate falls with temperature", {
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  # at 6 mg/L (z = 2.5) the thermal slope is b + dTDO z = 0.5 - 0.3 * 2.5 < 0
  falls <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = -0.3)
  expect_error(.sctmax_draws("exp", falls, 0.05, Tr = 30, costressors = list(do_cs),
                             levels = c(DO = 6)), class = "rTFT_falling_rate")
  # a milder negative interaction leaves the slope positive, 0.5 - 0.1 * 2.5 > 0
  mild <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = -0.1)
  expect_equal(.sctmax_draws("exp", mild, 0.05, Tr = 30, costressors = list(do_cs),
                             levels = c(DO = 6)),
               30 + (log(0.05) + 5 - 0.2 * 2.5) / (0.5 - 0.1 * 2.5), tolerance = 1e-7)
})


test_that("sCTmax needs one positive target per draw", {
  d <- list(loga = c(-5, -4.5), logb = log(c(0.5, 0.3)))
  expect_error(.sctmax_draws("exp", d, 0.05, Tr = 30), class = "rTFT_bad_target")
  expect_error(.sctmax_draws("exp", d, c(0.05, 0), Tr = 30), class = "rTFT_bad_target")
  expect_error(.sctmax_draws("exp", d, c(0.05, NA), Tr = 30), class = "rTFT_bad_target")
  expect_error(.sctmax_draws("exp", d["loga"], c(0.05, 0.1), Tr = 30),
               class = "rTFT_missing_draws")             # the shared checks apply
})


# ===========================================================================
# 5. Posterior draws from a fitted model (.fit_draws)
# ===========================================================================

# A small fitted model saved as a fixture (made by fixtures/make-hazard-fit.R): the
# exponential failure rate, memoryless baseline, body mass as a covariate on log a.
# 600 draws. The expected draws are rebuilt by hand from the fit's own coefficients:
# loga = b_loga_Intercept + mass * b_loga_mass, and logb = b_logb_Intercept.
fit_exp <- readRDS(test_path("fixtures", "hazard-fit-exp.rds"))
coefs   <- as.data.frame(posterior::as_draws_df(fit_exp))

test_that(".fit_draws gives one draws-by-rows matrix per parameter", {
  out <- .fit_draws(fit_exp, c("loga", "logb"), data.frame(mass = c(8, 12)))
  expect_named(out, c("loga", "logb"))
  expect_equal(dim(out$loga), c(600L, 2L))
  expect_equal(dim(out$logb), c(600L, 2L))
})


test_that(".fit_draws reads each parameter at the requested draws, paired", {
  # Draws 1, 5 and 9 of loga and of logb must come from the same three posterior
  # samples: checked against the fit's coefficients at those same rows.
  ids <- c(1, 5, 9)
  out <- .fit_draws(fit_exp, c("loga", "logb"), data.frame(mass = c(8, 12)),
                    draw_ids = ids)
  by_hand <- coefs[ids, ]
  expect_equal(unname(out$loga),
               cbind(by_hand$b_loga_Intercept + 8  * by_hand$b_loga_mass,
                     by_hand$b_loga_Intercept + 12 * by_hand$b_loga_mass))
  expect_equal(unname(out$logb[, 1]), by_hand$b_logb_Intercept)
})


test_that("thinning keeps the pairing", {
  ids <- .thin_ids(600, 50)
  out <- .fit_draws(fit_exp, c("loga", "logb"), data.frame(mass = 10), draw_ids = ids)
  expect_equal(unname(out$loga[, 1]),
               coefs$b_loga_Intercept[ids] + 10 * coefs$b_loga_mass[ids])
  expect_equal(unname(out$logb[, 1]), coefs$b_logb_Intercept[ids])
})


test_that(".thin_ids spaces the draws evenly and keeps the ends", {
  ids <- .thin_ids(600, 50)
  expect_length(ids, 50)
  expect_equal(range(ids), c(1, 600))
  expect_equal(.thin_ids(600), 1:600)                    # no thinning asked
  expect_equal(.thin_ids(600, 1000), 1:600)              # more than there are
  expect_error(.thin_ids(600, 0), class = "rTFT_bad_thin")
})


test_that("the columns filled in for brms do not reach the draws", {
  # Ta, Tr, log_exposure and failed are filled from the fitted data; supplying
  # quite different values must leave the draws unchanged.
  filled   <- .fit_draws(fit_exp, "loga", data.frame(mass = 10))
  supplied <- .fit_draws(fit_exp, "loga", data.frame(mass = 10, Ta = 50, Tr = 0,
                                                     log_exposure = 5, failed = 1))
  expect_equal(filled, supplied)
})


test_that("a covariate missing from newdata is an error, not a silent fill", {
  expect_error(.fit_draws(fit_exp, "loga", data.frame(size = 10)),
               class = "rTFT_missing_covariates")
  # logb has no covariates, so it needs none
  expect_silent(.fit_draws(fit_exp, "logb", data.frame(size = 10)))
})


test_that(".fit_draws checks the fit, the parameters and newdata", {
  expect_error(.fit_draws(list(), "loga", data.frame(mass = 10)), class = "rTFT_not_brmsfit")
  expect_error(.fit_draws(fit_exp, "Tc", data.frame(mass = 10)), class = "rTFT_missing_nlpar")
  expect_error(.fit_draws(fit_exp, "loga", data.frame()), class = "rTFT_bad_newdata")
})


test_that("a model fitted with fit_hazard carries its settings", {
  # The fixture was made by fit_hazard(), so this checks the whole path from the
  # fit to the record the reporting functions read.
  st <- .fit_settings(fit_exp)
  expect_equal(st[c("m_form", "baseline", "observation", "Tr", "k_sw", "costressors")],
               list(m_form = "exp", baseline = "mem", observation = "interval",
                    Tr = 20, k_sw = 0, costressors = list()))
  # the observed ranges: temperatures within the ramp (20 to 50 degC) and elapsed
  # times within the two-hour trial
  expect_true(st$Ta_range[1] >= 20 && st$Ta_range[2] <= 50)
  expect_true(st$tau_range[1] > 0 && st$tau_range[2] <= 120)
  # the individual of each row of the fitted data: 40 individuals
  expect_length(st$id, nrow(fit_exp$data))
  expect_equal(length(unique(st$id)), 40)
})


# ===========================================================================
# 6. The fitted failure rate, from a fitted model (M_draws)
# ===========================================================================

test_that("M_draws gives one row per draw, temperature and newdata row", {
  out <- M_draws(fit_exp, Ta = c(25, 30, 35), newdata = data.frame(mass = c(8, 12)),
                 n_draws = 10)
  expect_named(out, c("mass", ".draw", "Ta", "M"))
  expect_equal(nrow(out), 10 * 3 * 2)
  expect_equal(unique(out$.draw), .thin_ids(600, 10))    # the fit's own draw numbers
})


test_that("M_draws gives the fitted rate, computed from the fit's own draws", {
  # M = exp(loga + b (Ta - Tr)), with loga = b_loga_Intercept + mass b_loga_mass,
  # b = exp(b_logb_Intercept) and Tr = 20, all from the same draws.
  out <- M_draws(fit_exp, Ta = c(25, 35), newdata = data.frame(mass = c(8, 12)),
                 n_draws = 5)
  ids <- .thin_ids(600, 5)
  by_hand <- expand.grid(Ta = c(25, 35), draw = ids, mass = c(8, 12))
  cf <- coefs[by_hand$draw, ]
  by_hand$M <- exp(cf$b_loga_Intercept + by_hand$mass * cf$b_loga_mass +
                   exp(cf$b_logb_Intercept) * (by_hand$Ta - 20))
  expect_equal(out$M, by_hand$M)                         # same order: mass, draw, Ta
  expect_equal(out$.draw, by_hand$draw)
})


test_that("M_draws keeps every draw unless asked to thin", {
  out <- M_draws(fit_exp, Ta = 30, newdata = data.frame(mass = 10))
  expect_equal(nrow(out), 600)
})


test_that("M_draws needs the fit's covariates and a fit with its settings", {
  expect_error(M_draws(fit_exp, Ta = 30), class = "rTFT_missing_covariates")
  bare <- fit_exp
  bare$rtft <- NULL
  expect_error(M_draws(bare, Ta = 30, newdata = data.frame(mass = 10)),
               class = "rTFT_no_settings")
})


# ===========================================================================
# 7. Thermal tolerance from a fitted model (derive_tdt_ttl)
# ===========================================================================

test_that("derive_tdt_ttl gives four tables, each labelled by newdata", {
  out <- derive_tdt_ttl(fit_exp, newdata = data.frame(mass = c(8, 12)), p = c(0.5, 0.9),
                        t_ref = 60, n_draws = 50)
  expect_named(out, c("params", "sctmax", "tdt_curves", "ttl_surface"))
  expect_equal(names(out$tdt_curves),
               c("mass", "p", "Ta", "median", "lower", "upper", "prob", "n_draws", "n_unreached"))
  expect_equal(names(out$ttl_surface),
               c("mass", "Ta", "tau", "median", "lower", "upper", "prob", "n_draws", "n_unreached"))
  # the default grids: 50 temperatures, and 50 elapsed times for the landscape
  expect_equal(nrow(out$tdt_curves), 2 * 2 * 50)          # newdata rows x p x Ta
  expect_equal(nrow(out$ttl_surface), 2 * 50 * 50)        # newdata rows x Ta x tau
  expect_equal(names(out$params),
               c("mass", "quantity", "median", "lower", "upper", "prob", "n_draws", "n_unreached"))
  expect_equal(names(out$sctmax),
               c("mass", "p", "t_ref", "median", "lower", "upper", "prob", "n_draws", "n_unreached"))
  expect_equal(out$sctmax$mass, c(8, 8, 12, 12))
  expect_equal(out$sctmax$p, c(0.5, 0.9, 0.5, 0.9))
})


test_that("the parameters are summarised on their natural scale, with z where exact", {
  # exp with the memoryless baseline: a = exp(loga), b = exp(logb), z = log(10) / b,
  # from the same draws the function uses.
  out <- derive_tdt_ttl(fit_exp, newdata = data.frame(mass = 10), t_ref = 60,
                        n_draws = 50)
  cf  <- coefs[.thin_ids(600, 50), ]
  a <- exp(cf$b_loga_Intercept + 10 * cf$b_loga_mass)
  b <- exp(cf$b_logb_Intercept)
  expect_equal(out$params$quantity, c("a", "b", "z"))
  expect_equal(out$params[, 3:8],
               dplyr::bind_rows(.draw_summary(a), .draw_summary(b),
                                .draw_summary(log(10) / b)))
})


test_that("sCTmax matches the closed form for the exponential rate", {
  # Memoryless: H0(t_ref) = t_ref, so the target is -log(1 - p) / t_ref, and
  # loga + b (T - Tr) = log(target) gives T = Tr + (log(target) - loga) / b.
  out <- derive_tdt_ttl(fit_exp, newdata = data.frame(mass = 10), p = 0.5,
                        t_ref = 30, n_draws = 50)
  cf  <- coefs[.thin_ids(600, 50), ]
  loga <- cf$b_loga_Intercept + 10 * cf$b_loga_mass
  b    <- exp(cf$b_logb_Intercept)
  T_by_hand <- 20 + (log(-log(0.5) / 30) - loga) / b
  expect_equal(out$sctmax[, 4:9], .draw_summary(T_by_hand), tolerance = 1e-7)
})


test_that("derive_tdt_ttl refuses bad inputs, each with its own error class", {
  nd <- data.frame(mass = 10)
  expect_error(derive_tdt_ttl(fit_exp, nd, p = 1, t_ref = 60),           class = "rTFT_bad_p")
  expect_error(derive_tdt_ttl(fit_exp, nd, p = c(0.5, NA), t_ref = 60),  class = "rTFT_bad_p")
  expect_error(derive_tdt_ttl(fit_exp, nd),                             class = "rTFT_bad_t_ref")
  expect_error(derive_tdt_ttl(fit_exp, nd, t_ref = 0),                  class = "rTFT_bad_t_ref")
  expect_error(derive_tdt_ttl(fit_exp, nd, t_ref = 60, prob = 1),       class = "rTFT_bad_prob")
  expect_error(derive_tdt_ttl(fit_exp, t_ref = 60),                     class = "rTFT_missing_covariates")
  bare <- fit_exp
  bare$rtft <- NULL
  expect_error(derive_tdt_ttl(bare, nd, t_ref = 60),                    class = "rTFT_no_settings")
})


test_that("an unused co-stressor level is warned about once, not once per row and p", {
  n_warn <- 0
  withCallingHandlers(
    derive_tdt_ttl(fit_exp, newdata = data.frame(mass = c(8, 12)), p = c(0.5, 0.9),
                   t_ref = 60, levels = c(DO = 6), n_draws = 20),
    rTFT_unused_levels = function(w) {
      n_warn <<- n_warn + 1
      invokeRestart("muffleWarning")
    })
  expect_equal(n_warn, 1)
})


test_that(".param_draws reports z only where the curve is straight", {
  d <- list(loga = c(-5, -4), logb = log(c(0.5, 0.4)), Tc = c(34, 35),
            logshape = log(c(2, 3)), logscale = log(c(60, 70)), gam = c(0.01, -0.02))
  # exp + Weibull: the shape is reported, the scale is fixed and is not;
  # z = k log(10) / b
  q <- .param_draws("exp", "weibull", d, tau_ref = 60)
  expect_equal(names(q), c("a", "b", "shape", "z"))
  # a Weibull model with no recorded reference time cannot be reported
  expect_error(.param_draws("exp", "weibull", d), class = "rTFT_no_tau_ref")
  # exp + log-logistic: shape and scale both reported
  expect_equal(names(.param_draws("exp", "llogis", d)), c("a", "b", "shape", "scale"))
  # what is reported follows the function and baseline, not the draws present: d
  # also holds Tc and the baseline parameters, which exp + memoryless has none of
  expect_equal(names(.param_draws("exp", "mem", d)), c("a", "b", "z"))
  expect_equal(q$z, c(2, 3) * log(10) / c(0.5, 0.4))
  # exp + log-logistic: no single z
  expect_false("z" %in% names(.param_draws("exp", "llogis", d)))
  # exp + Gompertz: the rate constant, and no z
  q <- .param_draws("exp", "gompertz", d)
  expect_equal(q$gamma, c(0.01, -0.02))
  expect_false("z" %in% names(q))
  # a failure-rate function whose log is not linear in temperature: no z
  expect_false("z" %in% names(.param_draws("zeropower", "mem", d)))
})


test_that(".baseline_draws fixes the Weibull scale at the reference time", {
  d <- list(logshape = log(c(2, 3)), logscale = log(c(60, 70)), gam = c(0.01, -0.02))
  expect_equal(.baseline_draws("weibull", d, 2, tau_ref = 90),
               list(shape = c(2, 3), scale = c(90, 90)))
  expect_equal(.baseline_draws("llogis", d, 2, tau_ref = 90),
               list(shape = c(2, 3), scale = c(60, 70)))        # estimated, not fixed
  expect_equal(.baseline_draws("mem", d, 2), list(shape = c(1, 1), scale = c(1, 1)))
  expect_equal(.baseline_draws("gompertz", d, 2)$shape, c(0.01, -0.02))
  # so the cumulative baseline hazard at the reference time is 1 for every draw
  bl <- .baseline_draws("weibull", d, 2, tau_ref = 90)
  expect_equal(H0_of("weibull")(90, bl$shape, bl$scale), c(1, 1))
})


test_that(".param_draws uses the thermal sensitivity at the co-stressor levels", {
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                      interaction = TRUE)
  d <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = 0.03)
  q <- .param_draws("exp", "mem", d, list(do_cs), list(DO = 2.5))
  expect_equal(q$b_at_levels, 0.5 + 0.03 * 2.5)
  expect_equal(q$z, log(10) / (0.5 + 0.03 * 2.5))
  # without the interaction there is no b_at_levels, and z uses b itself
  no_int <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                       interaction = FALSE)
  q2 <- .param_draws("exp", "mem", d[c("loga", "logb", "cDO")], list(no_int),
                     list(DO = 2.5))
  expect_false("b_at_levels" %in% names(q2))
  expect_equal(q2$z, log(10) / 0.5)
})


test_that(".report_unreached gives one message listing the quantities, or none", {
  sct <- tibble::tibble(p = c(0.5, 0.9, 0.9), t_ref = 60, n_draws = 100,
                        n_unreached = c(0, 3, 12))
  par <- tibble::tibble(quantity = "a", n_draws = 100, n_unreached = 0)
  expect_message(.report_unreached(list(params = par, sctmax = sct)),
                 class = "rTFT_unreached")
  msg <- tryCatch(.report_unreached(list(params = par, sctmax = sct)),
                  message = conditionMessage)
  expect_match(msg, "sCTmax (p = 0.9, t_ref = 60): up to 12.0% of draws", fixed = TRUE)
  # quantities with nothing unreached are left out of it
  expect_false(grepl("p = 0.5", msg, fixed = TRUE))
  expect_false(grepl("a: up to", msg, fixed = TRUE))
  expect_silent(.report_unreached(list(params = par, sctmax = sct[1, ])))
})


test_that("the thermal death time curves match the closed form", {
  # Memoryless: H0(t) = t, so t_p = -log(1 - p) / M(T), with
  # M = exp(loga + b (Ta - 20)) from the same draws.
  out <- derive_tdt_ttl(fit_exp, newdata = data.frame(mass = 10), p = c(0.5, 0.9),
                        t_ref = 60, Ta = c(25, 30), n_draws = 50)
  cf   <- coefs[.thin_ids(600, 50), ]
  loga <- cf$b_loga_Intercept + 10 * cf$b_loga_mass
  b    <- exp(cf$b_logb_Intercept)
  by_hand <- dplyr::bind_rows(lapply(c(0.5, 0.9), \(pp) dplyr::bind_rows(lapply(c(25, 30), \(temp)
    .draw_summary(-log(1 - pp) / exp(loga + b * (temp - 20)))))))
  expect_equal(out$tdt_curves$p,  c(0.5, 0.5, 0.9, 0.9))
  expect_equal(out$tdt_curves$Ta, c(25, 30, 25, 30))
  expect_equal(out$tdt_curves[, 4:9], by_hand)
})


test_that("the tolerance landscape matches survival computed by hand", {
  # Memoryless: S = exp(-M(T) tau).
  out <- derive_tdt_ttl(fit_exp, newdata = data.frame(mass = 10), t_ref = 60,
                        Ta = c(25, 30), tau = c(10, 60), n_draws = 50)
  cf   <- coefs[.thin_ids(600, 50), ]
  loga <- cf$b_loga_Intercept + 10 * cf$b_loga_mass
  b    <- exp(cf$b_logb_Intercept)
  by_hand <- dplyr::bind_rows(lapply(c(25, 30), \(temp) dplyr::bind_rows(lapply(c(10, 60), \(tt)
    .draw_summary(exp(-exp(loga + b * (temp - 20)) * tt))))))
  expect_equal(out$ttl_surface$Ta,  c(25, 25, 30, 30))
  expect_equal(out$ttl_surface$tau, c(10, 60, 10, 60))
  expect_equal(out$ttl_surface[, 4:9], by_hand)
})


test_that("the grids must be finite temperatures and times that are not negative", {
  nd <- data.frame(mass = 10)
  expect_error(derive_tdt_ttl(fit_exp, nd, t_ref = 60, Ta = c(25, NA)),  class = "rTFT_bad_Ta")
  expect_error(derive_tdt_ttl(fit_exp, nd, t_ref = 60, tau = c(-1, 10)), class = "rTFT_bad_tau")
  expect_error(derive_tdt_ttl(fit_exp, nd, t_ref = 60, tau = "10"),     class = "rTFT_bad_tau")
})


test_that("the landscape starts at a survival of 1 at an elapsed time of zero", {
  out <- derive_tdt_ttl(fit_exp, data.frame(mass = 10), t_ref = 60, Ta = c(25, 30),
                        tau = c(0, 10), n_draws = 50)
  at_zero <- out$ttl_surface[out$ttl_surface$tau == 0, ]
  expect_equal(c(at_zero$median, at_zero$lower, at_zero$upper), rep(1, 6))
})


test_that("every baseline has no cumulative hazard at an elapsed time of zero", {
  # H0(0) = 0 for all four, which is what makes survival 1 there.
  expect_equal(H0_of("mem")(0, 1, 1), 0)
  expect_equal(H0_of("weibull")(0, 2, 60), 0)
  expect_equal(H0_of("llogis")(0, 2, 60), 0)
  expect_equal(H0_of("gompertz")(0, c(0.05, -0.05), 1), c(0, 0))
})


test_that(".tdt_draws turns a time that is never reached into Inf", {
  # Gompertz with gamma = -0.5: H0 levels off at 1 / 0.5 = 2. With p = 0.5 the
  # target is log(2) / M; a rate of 0.5 gives a target of 1.39 (reached), a rate
  # of 0.25 gives 2.77 (above the level: never reached).
  H0inv <- H0inv_of("gompertz")
  M <- matrix(c(0.5, 0.25), nrow = 1)                    # one draw, two temperatures
  tp <- .tdt_draws(M, 0.5, H0inv, shape = -0.5, scale = 1)
  expect_true(is.finite(tp[1, 1]))
  expect_equal(tp[1, 2], Inf)
  # a rate of zero (a zero function at Tr): the target is infinite, so never reached
  expect_equal(.tdt_draws(matrix(0), 0.5, H0inv_of("mem"), 1, 1)[1, 1], Inf)
})


test_that(".summarise_columns summarises each column on its own, in one pass", {
  # Four draws by three columns: ordinary numbers; one draw never reached; and half
  # the draws unattainable below, half above, where interpolating the median gives
  # Inf - Inf and the nearest draw must be taken instead, for that column only.
  mat <- cbind(c(3, 1, 4, 1), c(9, 2, Inf, 5), c(-Inf, -Inf, Inf, Inf))
  out <- .summarise_columns(mat)
  probs <- c(0.5, 0.025, 0.975)
  # column 1, against quantile() directly
  expect_equal(unlist(out[1, c("median", "lower", "upper")], use.names = FALSE),
               unname(stats::quantile(mat[, 1], probs)))
  # column 2: one of four draws unreached, so the upper bound is unbounded
  expect_equal(out$upper[2], Inf)
  expect_true(is.finite(out$median[2]))
  # column 3: no NaN, the median taken from the nearest draw
  expect_equal(c(out$median[3], out$lower[3], out$upper[3]), c(-Inf, -Inf, Inf))
  expect_equal(out$n_unreached, c(0L, 1L, 4L))
  expect_equal(out$n_draws, rep(4L, 3))
  expect_error(.summarise_columns(cbind(1, NA)), class = "rTFT_missing_draw_values")
})


test_that(".report_unreached lists a thermal death time curve once per p", {
  tdt <- tibble::tibble(p = c(0.9, 0.9), Ta = c(25, 30), n_draws = 100,
                        n_unreached = c(5, 40))
  msg <- tryCatch(.report_unreached(list(tdt_curves = tdt)), message = conditionMessage)
  expect_match(msg, "time until p = 0.9 fail, at some temperatures: up to 40.0% of draws",
               fixed = TRUE)
})


# ===========================================================================
# 8. The dynamic critical thermal maximum, under a linear ramp
# ===========================================================================

# Hand-made settings and two posterior draws, so dCTmax can be checked against the
# ramp's own integral without a fitted model.
ramp_settings <- function(m_form = "exp", baseline = "mem")
  list(m_form = m_form, baseline = baseline, Tr = 20, tau_ref = 60,
       k_sw = M_catalogue[[m_form]]$k_sw, costressors = list())
ramp_draws <- list(loga = c(-5, -4), logb = log(c(0.4, 0.5)), logc = c(-8, -8),
                   logshape = log(c(2, 0.7)), logscale = log(c(40, 80)), gam = c(0.01, -0.01))

test_that(".dctmax_draws solves the exponential ramp with the memoryless baseline exactly", {
  # H(t) = a exp(b (T0 - Tr)) (exp(b r t) - 1) / (b r), so the temperature at which
  # it reaches x is T0 + log(1 + x b r / (a exp(b (T0 - Tr)))) / b
  a <- exp(ramp_draws$loga); b <- exp(ramp_draws$logb)
  closed <- function(p, rate, T0)
    T0 + log1p(-log1p(-p) * b * rate / (a * exp(b * (T0 - 20)))) / b
  out <- .dctmax_draws(ramp_settings(), ramp_draws, rate = 0.25, T_start = 20,
                       p = c(0.1, 0.5, 0.9), levels = NULL, T_max = 45)
  expect_equal(dim(out), c(2L, 3L))
  expect_equal(out[, 2], closed(0.5, 0.25, 20), tolerance = 1e-5)
  expect_equal(out[, 1], closed(0.1, 0.25, 20), tolerance = 1e-5)
  expect_equal(out[, 3], closed(0.9, 0.25, 20), tolerance = 1e-5)
  # a start away from the reference
  cold <- .dctmax_draws(ramp_settings(), ramp_draws, 0.1, 14, 0.5, NULL, 45)
  expect_equal(cold[, 1], closed(0.5, 0.1, 14), tolerance = 1e-5)
  # more failing takes a higher temperature, and a slower ramp a lower one
  expect_true(all(out[, 1] < out[, 2] & out[, 2] < out[, 3]))
  slow <- .dctmax_draws(ramp_settings(), ramp_draws, 0.05, 20, 0.5, NULL, 45)
  expect_true(all(slow[, 1] < out[, 2]))
})


test_that(".dctmax_draws agrees with the ramp's integral under the other baselines", {
  # the cumulative hazard by direct numerical integration, for the first draw
  a <- exp(ramp_draws$loga[1]); b <- exp(ramp_draws$logb[1])
  k <- exp(ramp_draws$logshape[1]); s_ll <- exp(ramp_draws$logscale[1]); g <- ramp_draws$gam[1]
  rate <- 0.2; T0 <- 18
  h0 <- list(
    weibull  = \(t) (k / 60) * (t / 60)^(k - 1),          # the scale is the reference time
    llogis   = \(t) (k / s_ll) * (t / s_ll)^(k - 1) / (1 + (t / s_ll)^k),
    gompertz = \(t) exp(g * t))
  by_integral <- function(baseline, p) {
    hazard <- \(t) a * exp(b * (T0 + rate * t - 20)) * h0[[baseline]](t)
    gap <- \(t) stats::integrate(hazard, 0, t, rel.tol = 1e-10)$value + log1p(-p)
    T0 + rate * stats::uniroot(gap, c(1e-6, 200), tol = 1e-10)$root
  }
  for (baseline in names(h0)) {
    out <- .dctmax_draws(ramp_settings(baseline = baseline), ramp_draws, rate, T0,
                         p = c(0.2, 0.5), levels = NULL, T_max = 50)
    expect_equal(out[1, 1], by_integral(baseline, 0.2), tolerance = 1e-4)
    expect_equal(out[1, 2], by_integral(baseline, 0.5), tolerance = 1e-4)
  }
})


test_that(".dctmax_draws gives Inf where the ramp ends first, and handles the zero-type rates", {
  # followed only to 22 degrees, half have not failed
  short <- .dctmax_draws(ramp_settings(), ramp_draws, 0.25, 20, c(0.001, 0.5), NULL, T_max = 22)
  expect_true(all(is.finite(short[, 1])))
  expect_equal(short[, 2], c(Inf, Inf))
  # zeroexp below its reference accrues nothing, so a colder start changes little
  from_ref  <- .dctmax_draws(ramp_settings("zeroexp"), ramp_draws, 0.25, 20, 0.5, NULL, 45)
  from_cold <- .dctmax_draws(ramp_settings("zeroexp"), ramp_draws, 0.25, 15, 0.5, NULL, 45)
  expect_equal(from_cold, from_ref, tolerance = 1e-3)
  expect_error(.dctmax_draws(ramp_settings("power"), ramp_draws, 0.25, 15, 0.5, NULL, 45),
               class = "rTFT_below_Tr")
  # a ramp followed over a single step still gives one value per draw
  one_step <- .dctmax_draws(ramp_settings(), ramp_draws, 0.25, 20, 0.5, NULL, T_max = 20.005)
  expect_equal(dim(one_step), c(2L, 1L))
  expect_equal(one_step[, 1], c(Inf, Inf))
})


test_that("derive_dctmax reports dCTmax for each row of newdata, ramp and p", {
  nd <- data.frame(mass = c(8, 12))
  out <- suppressMessages(derive_dctmax(fit_exp, rate = c(0.1, 0.25), T_start = 20, newdata = nd,
                                        p = c(0.1, 0.5), n_draws = 200))
  expect_named(out, c("mass", "rate", "T_start", "p", "median", "lower", "upper", "prob",
                      "n_draws", "n_unreached"))
  expect_equal(nrow(out), 2L * 2L * 2L)
  expect_equal(out$mass, rep(c(8, 12), each = 4))
  expect_equal(out$rate, rep(rep(c(0.1, 0.25), each = 2), 2))
  expect_true(all(out$T_start == 20))                       # one start for every ramp
  expect_true(all(out$lower < out$median & out$median < out$upper))
  expect_true(all(out$n_unreached == 0))
  one <- function(mass, rate, p) out$median[out$mass == mass & out$rate == rate & out$p == p]
  expect_lt(one(8, 0.1, 0.5), one(8, 0.25, 0.5))            # the slower ramp, the lower
  expect_lt(one(8, 0.25, 0.1), one(8, 0.25, 0.5))
  expect_lt(one(12, 0.25, 0.5), one(8, 0.25, 0.5))          # heavier individuals fail sooner here
  # the median is that of the closed form over the draws
  dr <- .fit_draws(fit_exp, c("loga", "logb"), nd[1, , drop = FALSE],
                   draw_ids = .thin_ids(brms::ndraws(fit_exp), 200))
  a <- exp(dr$loga[, 1]); b <- exp(dr$logb[, 1])
  closed <- 20 + log1p(log(2) * b * 0.25 / (a * exp(b * (20 - fit_exp$rtft$Tr)))) / b
  expect_equal(one(8, 0.25, 0.5), stats::median(closed), tolerance = 1e-5)
  # a start per ramp
  two <- suppressMessages(derive_dctmax(fit_exp, rate = c(0.1, 0.25), T_start = c(18, 20),
                                        newdata = nd[1, , drop = FALSE], n_draws = 100))
  expect_equal(two$T_start, c(18, 20))
  expect_equal(two$p, c(0.5, 0.5))                           # p defaults to one half
})


test_that("derive_dctmax says when the ramp ends first or leaves the fitted temperatures", {
  nd <- data.frame(mass = 10)
  expect_message(short <- derive_dctmax(fit_exp, 0.25, 20, newdata = nd, p = 0.9, T_max = 24,
                                        n_draws = 100),
                 "do not reach", class = "rTFT_dctmax_notes")
  expect_true(all(short$n_unreached > 0))
  # a fast ramp takes dCTmax above the warmest temperature the model was fitted to
  top <- fit_exp$rtft$Ta_range[2]
  expect_message(fast <- derive_dctmax(fit_exp, 5, 20, newdata = nd, p = 0.9, n_draws = 100),
                 "above the temperatures", class = "rTFT_dctmax_notes")
  expect_gt(fast$median, top)
  expect_no_message(derive_dctmax(fit_exp, 0.25, 20, newdata = nd, p = 0.1, n_draws = 100))
})


test_that("derive_dctmax refuses bad arguments", {
  nd <- data.frame(mass = 10)
  expect_error(derive_dctmax(fit_exp, T_start = 20, newdata = nd), class = "rTFT_bad_rate")
  expect_error(derive_dctmax(fit_exp, rate = 0, T_start = 20, newdata = nd), class = "rTFT_bad_rate")
  expect_error(derive_dctmax(fit_exp, rate = -0.3, T_start = 20, newdata = nd), class = "rTFT_bad_rate")
  expect_error(derive_dctmax(fit_exp, rate = 0.3, newdata = nd), class = "rTFT_bad_T_start")
  expect_error(derive_dctmax(fit_exp, rate = c(0.1, 0.2, 0.3), T_start = c(18, 20), newdata = nd),
               class = "rTFT_bad_T_start")
  expect_error(derive_dctmax(fit_exp, 0.3, 20, newdata = nd, p = 1), class = "rTFT_bad_p")
  expect_error(derive_dctmax(fit_exp, 0.3, 20, newdata = nd, T_max = 19), class = "rTFT_bad_T_max")
  expect_error(derive_dctmax(fit_exp, 0.3, 20, newdata = nd, prob = 2), class = "rTFT_bad_prob")
  expect_error(derive_dctmax(fit_exp, 0.3, 20), class = "rTFT_missing_covariates")
  expect_error(derive_dctmax(list(), 0.3, 20), class = "rTFT_not_brmsfit")
})



# ===========================================================================
# 9. A fitted model of static and ramping trials together
# ===========================================================================

# A small fitted model saved as a fixture (made by fixtures/make-hazard-fit-mixed.R):
# three static trials checked every 5 minutes and two ramps watched throughout, with
# the exponential failure rate, the memoryless baseline and no covariate, simulated
# with log a = -5 at a reference of 20 degC and b = 0.4 per degC. 600 draws.
fit_mixed   <- readRDS(test_path("fixtures", "hazard-fit-mixed.rds"))
coefs_mixed <- as.data.frame(posterior::as_draws_df(fit_mixed))

test_that("sCTmax and dCTmax of a fit of both kinds of assay match their closed forms", {
  cf <- coefs_mixed[.thin_ids(600, 50), ]
  a  <- exp(cf$b_loga_Intercept)
  b  <- exp(cf$b_logb_Intercept)
  # sCTmax, as for the fit of one ramp: T = Tr + (log(-log(1 - p) / t_ref) - loga) / b
  out <- derive_tdt_ttl(fit_mixed, p = 0.5, t_ref = 60, n_draws = 50)
  expect_named(out, c("params", "sctmax", "tdt_curves", "ttl_surface"))
  T_by_hand <- 20 + (log(-log(0.5) / 60) - log(a)) / b
  summary_columns <- c("median", "lower", "upper", "prob", "n_draws", "n_unreached")
  expect_equal(out$sctmax[summary_columns], .draw_summary(T_by_hand), tolerance = 1e-7)
  # dCTmax under a ramp that starts at the reference temperature: one half have
  # failed when a (exp(b r t) - 1) / (b r) = log(2)
  ramp <- suppressMessages(
    derive_dctmax(fit_mixed, rate = 0.25, T_start = 20, n_draws = 50))
  closed <- 20 + log1p(log(2) * b * 0.25 / a) / b
  expect_equal(ramp$median, stats::median(closed), tolerance = 1e-5)
  expect_equal(ramp$n_unreached, 0L)
})


test_that("a fit of both kinds of assay covers the values the data were simulated from", {
  out <- derive_tdt_ttl(fit_mixed, p = 0.5, t_ref = 60)
  covers <- function(row, truth) row$lower < truth && truth < row$upper
  expect_true(covers(out$params[out$params$quantity == "a", ], exp(-5)))
  expect_true(covers(out$params[out$params$quantity == "b", ], 0.4))
  expect_true(covers(out$sctmax, 20 + (log(-log(0.5) / 60) + 5) / 0.4))
  # the curves and the landscape are defined over the whole of their grids
  expect_true(all(is.finite(out$tdt_curves$median)))
  expect_true(all(out$ttl_surface$median >= 0 & out$ttl_surface$median <= 1))
})
