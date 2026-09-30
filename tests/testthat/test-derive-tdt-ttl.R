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
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
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
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
  draws <- list(loga = -5, logb = log(0.5), cDO = 0.2, dTDO = 0.03)
  expect_error(.M_draws("exp", draws, Ta = 32, Tr = 30, costressors = list(do_cs)),
               class = "rTFT_bad_levels")
  expect_error(.M_draws("exp", draws, Ta = 32, Tr = 30, costressors = list(do_cs),
                       levels = c(pH = 7)), class = "rTFT_bad_levels")
})


test_that(".M_draws warns about levels the model has no co-stressor for", {
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
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
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
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
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
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
  do_cs <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
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
