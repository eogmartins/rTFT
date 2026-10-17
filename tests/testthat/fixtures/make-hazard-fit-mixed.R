# make-hazard-fit-mixed.R
#
# Fits one small failure-rate model to a table that holds static and ramping trials
# together, and saves it, so the tests can follow a fit of both kinds of assay
# through everything that is derived from a fitted model. Three static trials, at
# 23, 25 and 27 degC, are checked every 5 minutes for two hours, so their failures
# are known only to the interval between two checks, and an individual that has not
# failed by the last check is censored there. Two ramps from 20 degC, at 0.25 and at
# 0.1 degC per minute, are watched throughout, so their failures are at known times.
# The table is built with observation = "mixed", which gives each row the likelihood
# of its own kind of observation.
#
# The data are simulated from known values, the exponential failure rate with a
# memoryless baseline: the hazard is a exp(b (T - 20)), with log a = -5 and b = 0.4
# per degC. At a constant temperature the failure time is exponential with that
# rate. Under a ramp of rate r the cumulative hazard is
# H(t) = a (exp(b r t) - 1) / (b r), so each failure time is drawn by inverting H at
# a unit exponential, as in make-hazard-fit.R. The fit is short (two chains of 300
# draws after warm-up): the tests check how the fit is handled, not what it
# estimates.
#
# Unlike the other stored fits, the table is not kept at one row per check
# (grid_dT = Inf). A ramp watched throughout has no checks, so each of its
# individuals would be left with a single interval from the start to the failure.
# The intervals are divided so that temperature changes by at most 0.5 degC within
# one (grid_dT = 0.5), coarser than the default of 0.1 to keep the fit small, and
# each failure at a known time gets its short closing interval.
#
# Rerun this script when a brms update changes what a saved brmsfit holds, then
# rerun the tests. It is run from the package root.
#
# Inputs:
#   R/ (the package, loaded from source)
#
# Outputs:
#   tests/testthat/fixtures/hazard-fit-mixed.rds

# Setup ------------------------------------------------------------------------

library(pkgload)

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from tests/testthat/fixtures/", call. = FALSE)
}
load_all(quiet = TRUE)

# Simulate failure times in the static trials and under the ramps --------------

set.seed(20261013)
per_trial <- 10
check     <- 5
T_start   <- 20
loga_true <- -5
b_true    <- 0.4

# Static trials: a constant hazard, seen at the first check after the failure.
static_temp <- c(23, 25, 27)
static_max  <- 120
static <- data.frame(path_id = rep(sprintf("T%d", static_temp), each = per_trial),
                     temp    = rep(static_temp, each = per_trial))
static$t_fail <- stats::rexp(nrow(static)) /
  exp(loga_true + b_true * (static$temp - T_start))

# Ramps: the failure time inverts the cumulative hazard of a linear ramp.
ramp_rate <- c(0.25, 0.1)
ramp_max  <- c(120, 200)
ramp <- data.frame(path_id = rep(sprintf("r%s", ramp_rate), each = per_trial),
                   rate    = rep(ramp_rate, each = per_trial),
                   t_max   = rep(ramp_max, each = per_trial))
br <- b_true * ramp$rate
ramp$t_fail <- log1p(stats::rexp(nrow(ramp)) * br / exp(loga_true)) / br

ind <- rbind(
  data.frame(id       = seq_len(nrow(static)),
             path_id  = static$path_id,
             event    = as.integer(static$t_fail <= static_max),
             obs_type = "interval",
             t_upper  = pmin(ceiling(static$t_fail / check) * check, static_max),
             t_event  = NA_real_),
  data.frame(id       = nrow(static) + seq_len(nrow(ramp)),
             path_id  = ramp$path_id,
             event    = as.integer(ramp$t_fail <= ramp$t_max),
             obs_type = "continuous",
             t_upper  = NA_real_,
             t_event  = pmin(ramp$t_fail, ramp$t_max))
)
path <- rbind(
  data.frame(path_id = rep(sprintf("T%d", static_temp), each = 2),
             time    = rep(c(0, static_max), times = length(static_temp)),
             temp    = rep(static_temp, each = 2)),
  ramp_path(sprintf("r%s", ramp_rate), T_start = T_start, rate = ramp_rate,
            t_max = ramp_max)
)
checks <- data.frame(path_id = rep(sprintf("T%d", static_temp),
                                   each = static_max / check + 1),
                     time    = rep(seq(0, static_max, by = check),
                                   times = length(static_temp)))

long <- suppressMessages(
  build_hazard_data(ind, path, "mixed", Tr = T_start, checks = checks, grid_dT = 0.5)
)

# The table must hold what the tests rely on: failures of both kinds, and at least
# one individual that survived its static trial.
failure <- long$failed == 1
stopifnot(
  "no failure between two checks in the table" = any(failure & long$exact == 0),
  "no failure at a known time in the table"    = any(failure & long$exact == 1),
  "no individual survived a static trial"      = any(ind$event == 0),
  "an individual is missing from the table"    = setequal(long$id, ind$id)
)

# Fit and save -----------------------------------------------------------------

fit <- fit_hazard(long, "exp", observation = "mixed",
                  chains = 2, iter = 600, warmup = 300, seed = 20261013,
                  backend = "cmdstanr", refresh = 0)

stopifnot(brms::ndraws(fit) == 600)
saveRDS(fit, "tests/testthat/fixtures/hazard-fit-mixed.rds", compress = "xz")
