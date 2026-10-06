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
# among trials, rises with body mass and rises under low dissolved oxygen. The
# values are set so that at a constant temperature half the fish fail in about
# 20 min at 28 degC, 110 min at 26.5 degC and 10 h at 25 degC.
#
# The static assay holds 12 trials of 10 fish, four at each of three constant
# temperatures, checked every 5 min for two hours and every 30 min after that, to
# 24 h. Its failure times are therefore known to an interval. The ramping assay
# holds 12 trials of 10 fish warmed from 14 degC to 29.5 degC: three ramp rates by two
# levels of dissolved oxygen, two trials of each, watched throughout, so its failure
# times are known exactly. Either design can be observed either way; the two
# datasets differ in this only so that both ways are shown.
#
# The stream record is 30 days of hourly temperature with a daily cycle and a
# warm spell in the third week.
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
  b_do      = 0.3,     # change in log a per mg/L of dissolved oxygen below 10
  do_normal = 10       # dissolved oxygen in normoxia, mg/L
)
# log a above is the value at 26.5 degC; move it to the reference of 26 degC.
truth$log_a <- truth$log_a - truth$b * 0.5

# The time to failure of one fish under a temperature trajectory: the cumulative
# hazard is accumulated on a fine grid of time and the fish fails where it first
# exceeds a unit exponential draw. Returns Inf for a fish that survives the trial.
fail_time <- function(temp_at, t_max, log_a, step = 0.01) {
  t    <- seq(step, t_max, by = step)
  rate <- exp(log_a + truth$b * (temp_at(t - step / 2) - truth$T_ref))
  H0   <- (t / truth$scale)^truth$shape
  hazard <- cumsum(rate * diff(c(0, H0)))
  failed <- which(hazard >= stats::rexp(1))[1]
  if (is.na(failed)) Inf else t[failed]
}

# One trial's fish: their masses, and log a for each from the trial's own
# departure, the fish's mass and the trial's dissolved oxygen.
trial_fish <- function(n, do) {
  mass <- round(stats::runif(n, 5, 30), 1)
  list(mass = mass,
       log_a = truth$log_a + stats::rnorm(1, 0, truth$sd_trial) +
         truth$b_mass * (mass - 15) + truth$b_do * (truth$do_normal - do))
}

# The static assay -------------------------------------------------------------

static_design <- tibble::tibble(
  trial = sprintf("S%02d", 1:12),
  temp  = rep(c(25, 26.5, 28), each = 4))
static_end    <- 24 * 60
static_checks <- c(seq(0, 120, by = 5), seq(150, static_end, by = 30))

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

ramp_assay <- do.call(rbind, lapply(seq_len(nrow(ramp_design)), \(i) {
  rate  <- ramp_design$rate[i]
  t_max <- (ramp_stop - ramp_start) / rate
  fish  <- trial_fish(10, ramp_design$do[i])
  t_fail <- vapply(fish$log_a, \(la)
    fail_time(\(t) ramp_start + rate * t, t_max, la), numeric(1))
  tibble::tibble(
    trial   = ramp_design$trial[i],
    path_id = sprintf("r%s", rate),
    rate    = rate,
    do      = ramp_design$do[i],
    mass    = fish$mass,
    event   = as.integer(is.finite(t_fail)),
    # the time of failure to the nearest tenth of a minute, or the end of the ramp
    t_event = round(pmin(t_fail, t_max), 1))
}))
ramp_assay <- tibble::add_column(ramp_assay, id = seq_len(nrow(ramp_assay)), .before = 1)

# The stream record ------------------------------------------------------------

hours <- 0:(30 * 24)
day   <- hours / 24
stream_record <- tibble::tibble(
  time = hours * 60,
  day  = pmin(floor(day) + 1, 30),
  temp = round(
    17 + 1.5 * sin(pi * day / 30) +                    # the season: warmest mid-record
      2.2 * sin(2 * pi * (hours %% 24 - 9) / 24) +     # the day: warmest mid-afternoon
      5.5 * exp(-(day - 18)^2 / (2 * 2^2)) +           # a warm spell around day 18
      stats::rnorm(length(hours), 0, 0.15), 2))

# Checks on what was simulated -------------------------------------------------

stopifnot(
  nrow(static_assay) == 120, nrow(ramp_assay) == 120,
  !anyNA(static_assay), !anyNA(ramp_assay), !anyNA(stream_record),
  all(static_assay$t_upper %in% static_checks),
  # some fish, and not most, survive the coolest static trials
  mean(static_assay$event[static_assay$temp == 25]) > 0.5,
  mean(static_assay$event) < 1,
  # nearly every fish fails before a ramp ends
  mean(ramp_assay$event) > 0.9,
  all(diff(stream_record$time) == 60))

# Save -------------------------------------------------------------------------

attr(static_assay, "truth") <- truth
attr(ramp_assay, "truth")   <- truth

save(static_assay,  file = "data/static_assay.rda",  compress = "xz", version = 2)
save(ramp_assay,    file = "data/ramp_assay.rda",    compress = "xz", version = 2)
save(stream_record, file = "data/stream_record.rda", compress = "xz", version = 2)
