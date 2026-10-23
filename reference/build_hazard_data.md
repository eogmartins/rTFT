# Build the row-per-interval table a failure-rate model is fitted to

Reshapes one row per individual, together with a temperature trajectory,
into one row per individual and at-risk interval. An interval-censored
failure gives the failure its whole censoring interval; a known-time one
truncates the final at-risk interval at the recorded time. The
cumulative hazard is integrated numerically over the at-risk intervals,
with the failure-check times serving as the integration nodes. For a
known-time failure, the fit integrates on the checks (its final interval
ending exactly at the failure time) and no finer sub-grid is required
but can be provided.

## Usage

``` r
build_hazard_data(
  ind,
  path,
  observation = c("interval", "continuous", "mixed"),
  covariates = character(0),
  Tr = NULL,
  tau_ref = NULL,
  checks = NULL,
  grid_dt = NULL,
  grid_dT = 0.1,
  costressors = list()
)
```

## Arguments

- ind:

  A data frame with one row per individual. Required: `id`, `path_id`
  (key into `path`), `event` (1 = failed, 0 = right-censored), plus the
  covariate columns referenced in the model. Timing depends on how the
  failure was observed: interval-censored individuals have `t_upper`
  (the check at which failure was found, or a survivor's last observed
  check); known-time ones have `t_event` (the known time of failure or
  censoring). For `observation = "mixed"` the data also needs an
  `obs_type` column per individual (`"interval"` or `"continuous"`), and
  each type supplies its own timing column. Every individual is taken to
  be at risk from time 0 of its path: the table has no entry time, so an
  individual that joins an exposure after it began cannot be
  represented.

- path:

  A data frame giving the temperature trajectory for each path. Required
  columns: `path_id`, `time`, `temp`. The function interpolates these
  points linearly in time, so a minimal two-point path (start and end)
  encodes a linear ramp, and each path has one reading at each time: a
  sudden change is written as two readings a moment apart (e.g., at 10
  and 10.01). Interpolation uses the full series, so a temperature data
  logger whose records are finer than the checks still supplies the true
  temperature at each interval midpoint (never a check-to-check straight
  line); the checks control only the integration grid, not the
  temperature. Temperatures outside the recorded time range are held
  flat at the nearest endpoint rather than returned as `NA`, because the
  at-risk grid (the checks, the exact failure time, and the `t = 0`
  start) need not sit strictly inside the logged times: a check or
  failure time can fall just past the last reading, or the first
  interval can begin before the first reading, and an `NA` temperature
  there would propagate into `Ta` and the linear predictor.

- observation:

  How each failure time was observed. One of `"interval"` (every failure
  known only to fall within a censoring interval), `"continuous"` (every
  failure observed at a known time), or `"mixed"` (both kinds in one
  table, each individual's kind given by its `obs_type` column). This
  determines which response family the table is fitted with, and whether
  the `exact` flag is added; the response column is `failed` in every
  mode.

- covariates:

  A character vector of column names in `ind` to carry through onto
  every row of the output table, so they are available in the model
  formula as fixed-effect predictors or as grouping factors for random
  effects. A trial identifier is the common case: include it here, then
  place a trial-level random intercept on a parameter in
  [build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md),
  e.g. `dpar_formulas = list(loga = ~ 1 + (1 | trial))`. The individual
  identifier `id` is always kept, so an individual-level grouping term
  needs no entry here.

- Tr:

  The reference temperature, written to a constant `Tr` column so the
  functions can form `Ta - Tr`. It is a single model reference, not the
  per-trial starting temperature (see
  [M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md));
  the per-trial starts live in the temperature trajectories (`path`).
  For the forms that take a logarithm or fractional power of `Ta - Tr`
  (`zeroexp`, `power`, `zeropower`) it is the temperature at which the
  rate is zero (or, for `power`, falls to its baseline), so it must be
  given, at or below the coldest temperature in the data. Otherwise it
  defaults (`NULL`) to the median temperature at which failures were
  observed, which places `loga` where the data are and keeps it from
  being tied to the thermal sensitivity b; a reference far from the
  failures makes the two hard to tell apart. It must be set to a common
  value across datasets when comparing `loga`, whose meaning is anchored
  at `Tr`.

- tau_ref:

  A single positive number, the reference time, in the time unit of the
  data, written to a constant `tau_ref` column. The scale of the
  `weibull` baseline is fixed at it (see
  [h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md)),
  and the default priors of the parameters that depend on the time unit
  are scaled by it (see
  [build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)).
  It defaults (`NULL`) to the median time at which failures were
  observed, which places `loga` where the data are, as the default `Tr`
  does. A reference time much longer than the times at which individuals
  fail makes `loga` an extrapolation in time that depends heavily on the
  shape of the baseline. Like `Tr`, set it to a common value across
  datasets when comparing `loga` under a `weibull` baseline, whose
  meaning is anchored at it.

- checks:

  An optional data frame with `path_id`, `time` giving the failure-check
  schedule per path (i.e., the times at which each individual was
  inspected for failure). Defaults to the distinct `path` times (logger
  reading = check). Time 0, the start of the exposure, is added to every
  schedule when it is not listed. An interval-censored `t_upper` must be
  one of its path's checks: a time between two checks, or past the last
  one, is refused. As noted above, these times double as the integration
  nodes for the cumulative hazard, so their spacing sets the
  midpoint-rule resolution. For a known-time failure, a regular grid can
  be used instead (see `grid_dt`).

- grid_dt:

  An optional numeric giving a fixed time step (in the units of
  `path$time`). If set, a known-time failure is followed on a regular
  grid of this step, from `t = 0` to the failure time, instead of on
  `checks`. Defaults to `NULL` (the checks, or the readings of `path`).
  Either way the intervals are divided further where `grid_dT` requires
  it.

- grid_dT:

  A single positive number, the largest change in temperature allowed
  within an interval at risk, in degrees. An interval over which
  temperature changes by more is divided into equal parts (see *The time
  grid*). Defaults to 0.1; `Inf` keeps the intervals as given.

- costressors:

  An optional list of
  [costressor](https://eogmartins.github.io/rTFT/reference/costressor.md)
  objects declaring non-temperature stressors (dissolved oxygen and the
  like). Each adds a standardised column `z_<name>`, read from its own
  trajectory (interpolated to the interval midpoints) or, for a stressor
  held fixed within a trial, from a column of `ind`. Defaults to none.

## Value

A long-format data frame, one row per individual and at-risk interval,
holding the identifiers `id` and `path_id`; the model columns `Ta`,
`Tr`, `tau`, `log_tau`, `tau_ref`, `exposure`, `log_exposure`; the
carried covariates; one `z_<name>` column per co-stressor; the 0/1
response column `failed`; and, for mixed observation, the integer
`exact` flag the family branches on.

## The time grid

Within each interval at risk the hazard is taken as constant, at its
value at the midpoint of the interval. That is accurate while
temperature changes little within the interval, because the failure rate
rises steeply with temperature. With a thermal sensitivity of b = 1.7
per degree, for example, the rate is exp(1.7 x 0.1) = 1.19 times higher
at the end of an interval that spans 0.1 degrees than at its start, and
exp(1.7 x 0.5) = 2.3 times higher across one that spans half a degree;
the hazard accrued over the interval, taken at its midpoint, is then
understated by about (b dT)^2 / 24, that is by 0.1% over a tenth of a
degree and by 3% over half a degree, an error that grows with the square
of the change in temperature. So any interval over which temperature
changes by more than `grid_dT` is divided into equal parts, each with
its own row; a part that still changes by more, where the path bends
within the interval (a logged record, or a ramp held and then resumed),
is divided again; and a message says how many rows resulted.

The error left by the default step grows with the square of b. In
simulations of large ramping assays (eight assays of 3,000 individuals
at three ramp rates, fitted by maximum likelihood), the default of 0.1
degrees raised b and the Weibull shape by about 0.1% at b = 1.7 and by
about 1% at b = 4.7, and sCTmax by 0.007 degrees at b = 4.7; dCTmax was
not affected. Those are small against the uncertainty of an assay of a
hundred or so individuals. Where the fitted b is above about 3, a step
of 0.05 makes them about four times smaller, at the cost of twice the
rows and about twice the time to fit.

Two cases are treated apart. A failure at a known time needs the hazard
at the moment it happened, and the midpoint of the interval it fell in
is earlier and, in a ramp, cooler; so that interval is closed with a
very short one ending at the failure. This is the larger of the two
errors: without the closing interval, a ramp followed in intervals of a
minute underestimated b = 1.7 by about 4% in a simulation of eight large
assays, and with it b was recovered whatever the length of the
intervals. A failure known only to the interval between two checks keeps
that interval whole, since where in it the individual failed is not
known, and its contribution is the probability of failing anywhere in
it. Its hazard is still taken at its midpoint, and on a ramp checked at
long intervals that is the largest error of all, because the grid cannot
reduce it: in simulations of ramps of 0.1 to 0.3 degrees per minute
checked every 3 minutes (failure intervals of 0.3 to 0.9 degrees), b was
3% low and sCTmax 0.07 degrees low at b = 1.7, and 20% and 0.23 degrees
low at b = 4.7, whatever the step. A warning is given when temperature
changes by more than half a degree within such an interval. On a ramp,
check often enough that temperature changes by a few tenths of a degree
at most between checks, the fewer the higher b, or record the time of
each failure. A static assay is not affected, since temperature does not
change within its intervals.

When neither `checks` nor `grid_dt` is given for failures at known
times, the intervals are those between the readings of `path`, divided
as above, and into at least 20 per individual, so that a trajectory
given by its two ends is still followed closely.

## Examples

``` r
# A single linear ramp encoded as two points: 20 -> 35 degrees over 60 minutes,
# with failure checked for every minute.
path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))

# Interval-censored: one failure found at the 30-minute check, one survivor.
ind <- data.frame(id = 1:2, path_id = "A", event = c(1, 0), t_upper = c(30, 60))
build_hazard_data(ind, path, "interval", checks = checks)
#> ℹ The intervals at risk were divided further, into 268 rows from 90, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
#> # A tibble: 268 × 10
#>       id path_id   tau exposure    Ta failed    Tr log_tau tau_ref log_exposure
#>    <int> <chr>   <dbl>    <dbl> <dbl>  <int> <dbl>   <dbl>   <dbl>        <dbl>
#>  1     1 A       0.167    0.333  20.0      0  27.4  -1.79       30        -1.10
#>  2     1 A       0.5      0.333  20.1      0  27.4  -0.693      30        -1.10
#>  3     1 A       0.833    0.333  20.2      0  27.4  -0.182      30        -1.10
#>  4     1 A       1.17     0.333  20.3      0  27.4   0.154      30        -1.10
#>  5     1 A       1.5      0.333  20.4      0  27.4   0.405      30        -1.10
#>  6     1 A       1.83     0.333  20.5      0  27.4   0.606      30        -1.10
#>  7     1 A       2.17     0.333  20.5      0  27.4   0.773      30        -1.10
#>  8     1 A       2.5      0.333  20.6      0  27.4   0.916      30        -1.10
#>  9     1 A       2.83     0.333  20.7      0  27.4   1.04       30        -1.10
#> 10     1 A       3.17     0.333  20.8      0  27.4   1.15       30        -1.10
#> # ℹ 258 more rows

# Known-time: the failure time itself, integrated on the same checks.
ind_c <- data.frame(id = 1, path_id = "A", event = 1, t_event = 42)
build_hazard_data(ind_c, path, "continuous", checks = checks)
#> ℹ The intervals at risk were divided further, into 127 rows from 42, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
#> # A tibble: 127 × 10
#>       id path_id   tau exposure    Ta failed    Tr log_tau tau_ref log_exposure
#>    <dbl> <chr>   <dbl>    <dbl> <dbl>  <int> <dbl>   <dbl>   <dbl>        <dbl>
#>  1     1 A       0.167    0.333  20.0      0  30.5  -1.79       42        -1.10
#>  2     1 A       0.5      0.333  20.1      0  30.5  -0.693      42        -1.10
#>  3     1 A       0.833    0.333  20.2      0  30.5  -0.182      42        -1.10
#>  4     1 A       1.17     0.333  20.3      0  30.5   0.154      42        -1.10
#>  5     1 A       1.5      0.333  20.4      0  30.5   0.405      42        -1.10
#>  6     1 A       1.83     0.333  20.5      0  30.5   0.606      42        -1.10
#>  7     1 A       2.17     0.333  20.5      0  30.5   0.773      42        -1.10
#>  8     1 A       2.5      0.333  20.6      0  30.5   0.916      42        -1.10
#>  9     1 A       2.83     0.333  20.7      0  30.5   1.04       42        -1.10
#> 10     1 A       3.17     0.333  20.8      0  30.5   1.15       42        -1.10
#> # ℹ 117 more rows

# A regular integration grid instead of the checks (known-time only);
# the last interval still closes exactly on the failure time.
build_hazard_data(ind_c, path, "continuous", grid_dt = 5)
#> ℹ The intervals at risk were divided further, into 110 rows from 9, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
#> # A tibble: 110 × 10
#>       id path_id   tau exposure    Ta failed    Tr log_tau tau_ref log_exposure
#>    <dbl> <chr>   <dbl>    <dbl> <dbl>  <int> <dbl>   <dbl>   <dbl>        <dbl>
#>  1     1 A       0.192    0.385  20.0      0  30.5 -1.65        42       -0.956
#>  2     1 A       0.577    0.385  20.1      0  30.5 -0.550       42       -0.956
#>  3     1 A       0.962    0.385  20.2      0  30.5 -0.0392      42       -0.956
#>  4     1 A       1.35     0.385  20.3      0  30.5  0.297       42       -0.956
#>  5     1 A       1.73     0.385  20.4      0  30.5  0.549       42       -0.956
#>  6     1 A       2.12     0.385  20.5      0  30.5  0.749       42       -0.956
#>  7     1 A       2.5      0.385  20.6      0  30.5  0.916       42       -0.956
#>  8     1 A       2.88     0.385  20.7      0  30.5  1.06        42       -0.956
#>  9     1 A       3.27     0.385  20.8      0  30.5  1.18        42       -0.956
#> 10     1 A       3.65     0.385  20.9      0  30.5  1.30        42       -0.956
#> # ℹ 100 more rows

# Mixed: each individual carries its kind in obs_type; the unused timing
# column is left NA.
ind_m <- data.frame(id = 1:2, path_id = "A", event = 1,
                    obs_type = c("interval", "continuous"),
                    t_upper = c(30, NA), t_event = c(NA, 42))
build_hazard_data(ind_m, path, "mixed", checks = checks)
#> ℹ The intervals at risk were divided further, into 215 rows from 72, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
#> # A tibble: 215 × 11
#>       id path_id   tau exposure    Ta failed exact    Tr log_tau tau_ref
#>    <int> <chr>   <dbl>    <dbl> <dbl>  <int> <int> <dbl>   <dbl>   <dbl>
#>  1     1 A       0.167    0.333  20.0      0     0  28.9  -1.79       36
#>  2     1 A       0.5      0.333  20.1      0     0  28.9  -0.693      36
#>  3     1 A       0.833    0.333  20.2      0     0  28.9  -0.182      36
#>  4     1 A       1.17     0.333  20.3      0     0  28.9   0.154      36
#>  5     1 A       1.5      0.333  20.4      0     0  28.9   0.405      36
#>  6     1 A       1.83     0.333  20.5      0     0  28.9   0.606      36
#>  7     1 A       2.17     0.333  20.5      0     0  28.9   0.773      36
#>  8     1 A       2.5      0.333  20.6      0     0  28.9   0.916      36
#>  9     1 A       2.83     0.333  20.7      0     0  28.9   1.04       36
#> 10     1 A       3.17     0.333  20.8      0     0  28.9   1.15       36
#> # ℹ 205 more rows
#> # ℹ 1 more variable: log_exposure <dbl>

# A covariate carried through to every row, and a shared reference Tr.
ind$mass <- c(12.1, 9.8)
build_hazard_data(ind, path, "interval", covariates = "mass", Tr = 15, checks = checks)
#> ℹ The intervals at risk were divided further, into 268 rows from 90, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
#> # A tibble: 268 × 11
#>       id  mass path_id   tau exposure    Ta failed    Tr log_tau tau_ref
#>    <int> <dbl> <chr>   <dbl>    <dbl> <dbl>  <int> <dbl>   <dbl>   <dbl>
#>  1     1  12.1 A       0.167    0.333  20.0      0    15  -1.79       30
#>  2     1  12.1 A       0.5      0.333  20.1      0    15  -0.693      30
#>  3     1  12.1 A       0.833    0.333  20.2      0    15  -0.182      30
#>  4     1  12.1 A       1.17     0.333  20.3      0    15   0.154      30
#>  5     1  12.1 A       1.5      0.333  20.4      0    15   0.405      30
#>  6     1  12.1 A       1.83     0.333  20.5      0    15   0.606      30
#>  7     1  12.1 A       2.17     0.333  20.5      0    15   0.773      30
#>  8     1  12.1 A       2.5      0.333  20.6      0    15   0.916      30
#>  9     1  12.1 A       2.83     0.333  20.7      0    15   1.04       30
#> 10     1  12.1 A       3.17     0.333  20.8      0    15   1.15       30
#> # ℹ 258 more rows
#> # ℹ 1 more variable: log_exposure <dbl>
```
