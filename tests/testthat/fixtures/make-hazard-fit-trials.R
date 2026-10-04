# make-hazard-fit-trials.R
#
# Fits one small failure-rate model in which individuals share a group-level term,
# and saves it, for the cross-validation tests. Eight trials of six individuals each
# are run on the same linear ramp, and each trial has its own random intercept on
# log a, so the individuals of a trial share that trial's effect. Leaving out one
# individual then still leaves its trial's effect in the fit, which is the case
# where cross-validation must leave out whole trials instead.
#
# The data are simulated as in make-hazard-fit.R: a ramp from 20 degC at 0.25 degC
# per minute, checked every 5 minutes for two hours, with each failure time drawn
# by inverting the cumulative hazard H(t) = a (exp(b r t) - 1) / (b r) at a unit
# exponential. Here log a = -5 plus a trial effect drawn with a standard deviation
# of 0.6. The fit is short (two chains of 300 draws after warm-up): the tests check
# how the cross-validation is set up and summed, not what it estimates.
#
# Rerun this script when a brms update changes what a saved brmsfit holds, then
# rerun the tests. It is run from the package root.
#
# Inputs:
#   R/ (the package, loaded from source)
#
# Outputs:
#   tests/testthat/fixtures/hazard-fit-trials.rds

# Setup ------------------------------------------------------------------------

library(pkgload)

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from tests/testthat/fixtures/", call. = FALSE)
}
load_all(quiet = TRUE)

# Simulate failure times by trial ----------------------------------------------

set.seed(20261002)
n_trials    <- 8
per_trial   <- 6
t_max       <- 120
check       <- 5
T_start     <- 20
rate        <- 0.25
b_true      <- 0.4
trial_sd    <- 0.6

trial  <- rep(sprintf("T%d", seq_len(n_trials)), each = per_trial)
effect <- stats::rnorm(n_trials, 0, trial_sd)
loga   <- -5 + effect[match(trial, unique(trial))]
br     <- b_true * rate
t_fail <- log1p(stats::rexp(length(trial)) * br / exp(loga)) / br

ind <- data.frame(
  id      = seq_along(trial),
  path_id = "ramp",
  event   = as.integer(t_fail <= t_max),
  t_upper = pmin(ceiling(t_fail / check) * check, t_max),
  trial   = trial
)
path   <- ramp_path("ramp", T_start = T_start, rate = rate, t_max = t_max)
checks <- data.frame(path_id = "ramp", time = seq(0, t_max, by = check))

long <- build_hazard_data(ind, path, "interval", covariates = "trial",
                          Tr = T_start, checks = checks)

# Fit and save -----------------------------------------------------------------

fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + (1 | trial)),
                  chains = 2, iter = 600, warmup = 300, seed = 20261002,
                  backend = "cmdstanr", refresh = 0,
                  control = list(adapt_delta = 0.95))

stopifnot(brms::ndraws(fit) == 600)
saveRDS(fit, "tests/testthat/fixtures/hazard-fit-trials.rds", compress = "xz")
