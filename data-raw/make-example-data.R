# make-example-data.R - the example datasets shipped with the package
#
# Simulates the three datasets the help pages, the README and the vignettes use:
# a static assay, a ramping assay and a record of stream temperature. They are
# simulated, not measured, so that they can be shipped freely and so that what a
# model estimates from them can be set against the values they were made from.
#
# Both assays are simulated from one model of a juvenile salmonid (the values are
# plausible for a juvenile rainbow trout acclimated to 14 degC, and stand for no
# particular study). The hazard of losing equilibrium is h0(tau) M(T), with an
# exponential failure rate M(T) = a exp(b (T - 26)) and a Weibull baseline
# h0(tau) = (k / s) (tau / s)^(k - 1), times in minutes. The amplitude a varies
# among trials and rises with body mass. The values are set so that at a constant
# temperature half the fish fail in about 1 h at 27 degC, 10 h at 25 degC and 4 days
# at 23 degC.
#
# Low dissolved oxygen raises the failure rate, and the more so the warmer the water:
# each mg/L of oxygen below 10 adds b_do + d_do (T - 26) to log M(T). The effect is
# 0.3 per mg/L at 26 degC and changes by 0.04 for each degree (0.38 at 28 degC, 0.22
# at 24 degC), so it is a straight line in temperature, the form that the interaction
# of a co-stressor with temperature fits. A straight line crosses zero somewhere, and
# the slope is small so that it does so at 18.5 degC, below the temperatures at which
# the stream record is warm enough for fish to fail.
#
# The static assay holds 12 trials of 10 fish, four at each of three constant
# temperatures (23, 25 and 27 degC), checked every hour for five days. Its failure times
# are therefore known to an interval. The coolest temperature is within reach of a
# summer heat wave, which lets the assay inform a prediction for the field. The ramping
# assay holds 32 trials of 10 fish warmed from 14 degC to 29.5 degC, watched throughout,
# so its failure times are known exactly. Four ramp rates (0.025, 0.1, 0.2 and 0.3 degC
# per min) are crossed with two levels of dissolved oxygen, four trials of each. The
# slow rate is there so that the same data can show what a narrow set of rates (0.1 to
# 0.3, a 3-fold span) and a wide one (0.025 to 0.3, a 12-fold span) allow a model to
# estimate. Either design can be observed either way; the two datasets differ in this
# only so that both ways are shown.
#
# The stream record is 30 days of hourly temperature with a daily cycle and a
# heat wave in the third week that takes the afternoon maximum to about 24 degC.
#
# Rerun this script from the package root to remake the datasets.
#
# Inputs:
#   none
#
# Outputs:
#   data/static_assay.rda
#   data/ramp_assay.rda
#   data/stream_record.rda

# Setup ------------------------------------------------------------------------

library(tibble)

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from data-raw/", call. = FALSE)
}

set.seed(20261006)

# The model the assays are simulated from --------------------------------------

truth <- list(
  T_ref     = 26,      # temperature at which M(T) = a, degC
  log_a     = log(log(2) / (110 / 60)^1.5),   # half fail in 110 min at 26.5 degC ...
  b         = 1.7,     # ... with this thermal sensitivity, per degC
  shape     = 1.5,     # Weibull shape k: the hazard rises with elapsed exposure
  scale     = 60,      # Weibull scale s, min
  sd_trial  = 0.3,     # among-trial standard deviation of log a
  b_mass    = 0.03,    # change in log a per gram above 15 g
  b_do      = 0.3,     # added to log M per mg/L of dissolved oxygen below 10, at T_ref
  d_do      = 0.04,    # change in that effect per degC
  do_normal = 10       # dissolved oxygen in normoxia, mg/L
)
# log a above is the value at 26.5 degC; move it to the reference of 26 degC.
truth$log_a <- truth$log_a - truth$b * 0.5

# The time to failure of one fish under a temperature trajectory: the cumulative
# hazard is accumulated on a fine grid of time and the fish fails where it first
# exceeds a unit exponential draw. Returns Inf for a fish that survives the trial.
fail_time <- function(temp_at, t_max, log_a, deficit = 0, step = 0.01) {
  t    <- seq(step, t_max, by = step)
  dT   <- temp_at(t - step / 2) - truth$T_ref
  rate <- exp(log_a + truth$b * dT + deficit * (truth$b_do + truth$d_do * dT))
  H0   <- (t / truth$scale)^truth$shape
  hazard <- cumsum(rate * diff(c(0, H0)))
  failed <- which(hazard >= stats::rexp(1))[1]
  if (is.na(failed)) Inf else t[failed]
}

# One trial's fish: their masses, log a for each from the trial's own departure and
# the fish's mass, and the trial's oxygen deficit.
trial_fish <- function(n, do) {
  mass <- round(stats::runif(n, 5, 30), 1)
  list(mass = mass,
       log_a = truth$log_a + stats::rnorm(1, 0, truth$sd_trial) +
         truth$b_mass * (mass - 15),
       deficit = truth$do_normal - do)
}

# The static assay -------------------------------------------------------------

static_design <- tibble::tibble(
  trial = sprintf("S%02d", 1:12),
  temp  = rep(c(23, 25, 27), each = 4))
static_end    <- 5 * 24 * 60
static_checks <- seq(0, static_end, by = 60)

static_assay <- do.call(rbind, lapply(seq_len(nrow(static_design)), \(i) {
  fish <- trial_fish(10, truth$do_normal)
  t_fail <- vapply(fish$log_a, \(la)
    fail_time(\(t) rep(static_design$temp[i], length(t)), static_end, la), numeric(1))
  tibble::tibble(
    trial   = static_design$trial[i],
    path_id = sprintf("T%s", static_design$temp[i]),
    temp    = static_design$temp[i],
    mass    = fish$mass,
    event   = as.integer(is.finite(t_fail)),
    # the check at which the fish was found to have failed, or the last check
    t_upper = vapply(t_fail, \(t)
      if (is.finite(t)) min(static_checks[static_checks >= t]) else static_end,
      numeric(1)))
}))
static_assay <- tibble::add_column(static_assay, id = seq_len(nrow(static_assay)),
                                   .before = 1)

# The ramping assay ------------------------------------------------------------

ramp_start <- 14
ramp_stop  <- 29.5
ramp_design <- tibble::tibble(
  trial = sprintf("R%02d", 1:12),
  rate  = rep(c(0.1, 0.2, 0.3), each = 4),
  do    = rep(rep(c(10, 5), each = 2), times = 3))

# The fish of the trials in a design, one row per fish.
simulate_ramps <- function(design) {
  do.call(rbind, lapply(seq_len(nrow(design)), \(i) {
    rate  <- design$rate[i]
    t_max <- (ramp_stop - ramp_start) / rate
    fish  <- trial_fish(10, design$do[i])
    t_fail <- vapply(fish$log_a, \(la)
      fail_time(\(t) ramp_start + rate * t, t_max, la, fish$deficit), numeric(1))
    tibble::tibble(
      trial   = design$trial[i],
      path_id = sprintf("r%s", rate),
      rate    = rate,
      do      = design$do[i],
      mass    = fish$mass,
      event   = as.integer(is.finite(t_fail)),
      # the time of failure to the nearest tenth of a minute, or the end of the ramp
      t_event = round(pmin(t_fail, t_max), 1))
  }))
}

ramp_first <- simulate_ramps(ramp_design)

# The stream record ------------------------------------------------------------

hours <- 0:(30 * 24)
day   <- hours / 24
stream_record <- tibble::tibble(
  time = hours * 60,
  day  = pmin(floor(day) + 1, 30),
  temp = round(
    17 + 1.5 * sin(pi * day / 30) +                    # the season: warmest mid-record
      2.2 * sin(2 * pi * (hours %% 24 - 9) / 24) +     # the day: warmest mid-afternoon
      3.5 * exp(-(day - 18)^2 / (2 * 2^2)) +           # a heat wave around day 18
      stats::rnorm(length(hours), 0, 0.15), 2))

# The further ramping trials --------------------------------------------------

# Two more trials of each combination of the three faster rates and the two levels
# of oxygen, four trials at the slow rate in normoxia and, last, four at the slow
# rate in hypoxia.
ramp_design_more <- tibble::tibble(
  trial = sprintf("R%02d", 13:32),
  rate  = c(rep(c(0.1, 0.2, 0.3), each = 4), rep(0.025, 8)),
  do    = c(rep(rep(c(10, 5), each = 2), times = 3), rep(10, 4), rep(5, 4)))

ramp_assay <- rbind(ramp_first, simulate_ramps(ramp_design_more))
ramp_assay <- tibble::add_column(ramp_assay, id = seq_len(nrow(ramp_assay)), .before = 1)

# Checks on what was simulated -------------------------------------------------

stopifnot(
  nrow(static_assay) == 120, nrow(ramp_assay) == 320,
  # four trials of ten fish of each rate at each level of oxygen
  all(table(ramp_assay$rate, ramp_assay$do) == 40),
  # low oxygen raises the failure rate at every temperature of the stream record
  # that is warm enough to matter, and at every temperature at which a fish failed
  truth$T_ref - truth$b_do / truth$d_do < 19,
  !anyNA(static_assay), !anyNA(ramp_assay), !anyNA(stream_record),
  all(static_assay$t_upper %in% static_checks),
  # at the coolest static temperature a good share of the fish fail and a good share
  # survive the five days, and at the two warmer ones nearly all fail
  mean(static_assay$event[static_assay$temp == 23]) > 0.25,
  mean(static_assay$event[static_assay$temp == 23]) < 0.75,
  mean(static_assay$event[static_assay$temp > 23]) > 0.9,
  # the heat wave reaches the coolest static temperature and stays below the next
  max(stream_record$temp) > 23, max(stream_record$temp) < 25,
  # nearly every fish fails before a ramp ends
  mean(ramp_assay$event) > 0.9,
  all(diff(stream_record$time) == 60))

# Save -------------------------------------------------------------------------

attr(static_assay, "truth") <- truth
attr(ramp_assay, "truth")   <- truth

save(static_assay,  file = "data/static_assay.rda",  compress = "xz", version = 2)
save(ramp_assay,    file = "data/ramp_assay.rda",    compress = "xz", version = 2)
save(stream_record, file = "data/stream_record.rda", compress = "xz", version = 2)
