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
