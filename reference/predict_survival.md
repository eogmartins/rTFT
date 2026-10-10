# Survival along a record of temperature, from a fitted failure-rate model

Predicts, from a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md),
the survival of an individual exposed to a measured record of
temperature: the probability of surviving to each time of the record and
of failing by its end, each with its posterior median and credible
interval, at the covariate values in each row of `newdata`. Random
effects are set to zero, so the prediction is for a typical group.

## Usage

``` r
predict_survival(
  fit,
  record,
  newdata = NULL,
  levels = NULL,
  by = NULL,
  prob = 0.95,
  n_draws = 1000,
  grid_dT = 0.1
)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- record:

  A data frame with the columns `time` and `temp`, the record of
  temperature through time, in the time unit and temperature scale of
  the data the model was fitted to, with times increasing and no gaps.
  It may also hold the column of a co-stressor and the column named in
  `by`.

- newdata:

  A data frame with one row per set of covariate values to predict at,
  holding every covariate in the parameters' formulas. `NULL` (the
  default) suits a model with no covariates.

- levels:

  A named numeric vector giving the level of each co-stressor that
  `record` does not follow, on its measurement scale, e.g. `c(DO = 6)`.

- by:

  A character string, the name of a column of `record` giving the period
  each reading belongs to (e.g., a day). When given, the result has each
  period's share of the hazard accrued over the record. An interval
  between two readings belongs to the period of the first.

- prob:

  A single number between 0 and 1, the coverage of the credible
  intervals. Defaults to 0.95.

- n_draws:

  A single number, the number of posterior draws to use, evenly spaced
  through the posterior. Defaults to 1000. The draws are held for every
  step along the record, so a long record with many draws needs memory
  in proportion to their product.

- grid_dT:

  A single positive number, the largest change in temperature within a
  step along the record, in degrees. Between two readings the
  temperature is taken to change linearly, and the hazard is accrued in
  steps over which it changes by at most this much, since the failure
  rate at the midpoint of a step stands for the whole of it (see *The
  time grid* in
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)).
  Defaults to 0.1; `Inf` takes one step between consecutive readings.

## Value

A list of tibbles, each with the columns of `newdata`:

- `survival`: for each time of the record, the posterior `mean`,
  `median`, `lower` and `upper` of the probability of surviving to it (1
  at the first time).

- `mortality`: the same for the probability of failing by the end of the
  record.

- `risk` (when `by` is given): each period's share of the hazard accrued
  over the record. The shares of a posterior draw add to one; `n_draws`
  is the number of draws with any hazard, which the summary rests on.

- `coverage`: the lowest and highest temperatures at which failures were
  observed in the data (`Ta_failed_min`, `Ta_failed_max`); the share of
  the record's duration `below`, `within` and `above` them; the share of
  the hazard accrued along the record at those temperatures
  (`risk_below`, `risk_within`, `risk_above`); the longest exposure in
  the data (`tau_max`); and the share of the record's duration that
  comes after it (`beyond`) with the share of the hazard accrued there
  (`risk_beyond`). Each share of the hazard is the median over the
  posterior draws with any hazard.

- `record`: the times and temperatures of the record, with any
  co-stressor it follows, kept so that
  [`plot_survival()`](https://eogmartins.github.io/rTFT/reference/plot_survival.md)
  can draw them above the survival they led to.

## Details

Along the record the cumulative hazard grows, over each step, by M(T)
(H0(tau_1) - H0(tau_0)): the failure rate at the temperature in the
middle of the step, times the growth of the cumulative baseline hazard
over it. A step is the interval between two readings, or a part of it
where the temperature changes by more than `grid_dT` between them.
Survival to a time is exp(-cumulative hazard up to it). With the
memoryless baseline the growth of H0 is the length of the interval, so
the cumulative hazard is the failure rate summed over the record, and it
does not matter where the record begins. With the other baselines the
hazard depends on elapsed exposure, which is counted from the first time
in the record.

## Outside the temperatures and exposures of the assay

An assay is run at temperatures that cause failure within its duration,
and a field record may lie below them. At the temperatures at which
failures were observed the data hold the failure-rate function in place,
so functions that fit equally well predict alike. Below them nothing
does, and the prediction depends on which function was chosen: functions
that fit the assay equally well can differ by orders of magnitude there.
The temperatures that count are those at which individuals failed, and
not all those they were exposed to. In a ramping assay every individual
starts at a temperature at which none fails, and the time spent there
shows only that failures were rare over that time.

The `coverage` table gives the share of the record below, within and
above the temperatures at which failures were observed, by duration and
by the risk accrued. The second is the one to read: a record may spend
most of its time at cool temperatures where little hazard accrues, and
take nearly all of its risk from a few warm hours that the assay did
cover.

A record may also last longer than any exposure in the assay. With a
baseline other than the memoryless one, the prediction then rests on how
the baseline extends beyond the exposures the data hold, and with the
memoryless one on the assumption that the hazard does not change with
the time under exposure, which the data support only over those
exposures. The `coverage` table also gives the longest exposure in the
data and the share of the record, by duration and by the risk accrued,
that comes after it.

A message is given when any of the record lies outside the temperatures
or after the longest exposure. Comparing the predictions of several
failure-rate functions and baselines, or averaging them with
[`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md),
shows how much the prediction depends on those choices.

The zero-type failure rates (`zeroexp`, `zeropower`) are zero at the
reference temperature `Tr` and are taken as zero below it. The `power`
function is not defined below `Tr`, and a record that goes below it is
refused.

## The mean and the median

Each posterior draw gives a probability of surviving, and the tables
summarise those probabilities by their mean and by their median. The
mean is the probability of surviving with the uncertainty in the
parameters taken into account (and, from
[`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md),
the uncertainty about the model, as far as the weights measure it): it
is the fraction of individuals expected to survive in a typical group
(random effects at zero). The median is the value that half of the
posterior lies below. The two are close when the posterior is
concentrated, and differ when it is skewed or has two separate groups of
values, as it can when models that disagree are averaged; the credible
interval should then be reported with either.

## No recovery

The hazard accrued along the record is only ever added to, so the
prediction assumes that an individual does not recover when the
temperature returns to a range in which repair is expected to exceed
thermal damage. Along a record that moves in and out of stressful
temperatures, such as one with a daily cycle, the predicted mortality is
then too high if individuals recover in between, and the more so the
more cycles there are.

## Co-stressors

For a model with co-stressors, a co-stressor whose column is in `record`
is followed through time, taking in each interval the value halfway
between the two readings. One that is not in `record` is held at the
level given in `levels`.

## Examples

``` r
if (FALSE) { # \dontrun{
# a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
# degrees over 60 minutes and checked every 10, in two trials (see fit_hazard())
path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))
ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
                     trial = rep(c("a", "b"), each = 10),
                     t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
                            checks = checks)
fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))

pred <- predict_survival(fit, record, newdata = data.frame(mass = 10), by = "day")
pred$mortality
pred$risk
} # }
```
