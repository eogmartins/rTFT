# make-hazard-fit.R
#
# Fits one small failure-rate model and saves it, so the tests that read posterior
# draws out of a fitted model can run against a real brmsfit without compiling or
# sampling anything. The model is the exponential failure rate with a memoryless
# baseline, interval-censored, with body mass as a covariate on log a, so the tests
# can check that covariate values in newdata reach the parameter they belong to.
#
# The data are simulated from known values under a linear ramp from 20 degC at
# 0.25 degC per minute, checked every 5 minutes for two hours. Under that ramp the
# cumulative hazard is H(t) = a (exp(b r t) - 1) / (b r), so each failure time is
# drawn by inverting H at a unit exponential. The fit is short (two chains of 300
# draws after warm-up): the tests check how draws are extracted, not what they
# estimate.
#
# Rerun this script when a brms update changes what a saved brmsfit holds, then
# rerun the tests. It is run from the package root.
#
# Inputs:
#   R/ (the package, loaded from source)
#
# Outputs:
#   tests/testthat/fixtures/hazard-fit-exp.rds

# Setup ------------------------------------------------------------------------

library(pkgload)

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from tests/testthat/fixtures/", call. = FALSE)
}
load_all(quiet = TRUE)

# Simulate failure times under the ramp ----------------------------------------

set.seed(20260930)
n_fish <- 40
t_max  <- 120
check  <- 5
T_start <- 20
rate    <- 0.25
b_true  <- 0.4

mass   <- round(stats::runif(n_fish, 5, 15), 1)
loga   <- -5 + 0.2 * (mass - 10) 
br     <- b_true * rate
t_fail <- log1p(stats::rexp(n_fish) * br / exp(loga)) / br

ind <- data.frame(
  id      = seq_len(n_fish),
  path_id = "ramp",
  event   = as.integer(t_fail <= t_max),
  t_upper = pmin(ceiling(t_fail / check) * check, t_max),
  mass    = mass
)
path   <- ramp_path("ramp", T_start = T_start, rate = rate, t_max = t_max)
checks <- data.frame(path_id = "ramp", time = seq(0, t_max, by = check))

long <- build_hazard_data(ind, path, "interval", covariates = "mass",
                          Tr = T_start, checks = checks)

# Fit and save -----------------------------------------------------------------

fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass),
                  chains = 2, iter = 600, warmup = 300, seed = 20260930,
                  backend = "cmdstanr", refresh = 0)

stopifnot(brms::ndraws(fit) == 600)
saveRDS(fit, "tests/testthat/fixtures/hazard-fit-exp.rds", compress = "xz")
