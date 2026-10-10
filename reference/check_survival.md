# Check a fitted model against the survival observed in each group

Compares, for each group of individuals, the proportion observed to
survive through time with the proportion surviving in datasets simulated
from the fitted model: a posterior predictive check. Where the model
describes the data, the observed proportions lie within the band of the
simulated ones.

## Usage

``` r
check_survival(fit, by = NULL, prob = 0.95, n_draws = 200)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- by:

  An optional character string naming a variable of the fitted model
  with one value per individual (a grouping such as a trial, or a
  covariate with few values), to split the individuals of each
  trajectory by.

- prob:

  A single number between 0 and 1, the coverage of the band of simulated
  proportions. Defaults to 0.95.

- n_draws:

  A single number, the number of posterior draws to simulate from,
  evenly spaced through the posterior. Defaults to 200.

## Value

A tibble with a row for each group and time: `path_id`, the `by`
variable when given, `n` (the number of individuals in the group),
`time`, `Ta` (the temperature of the trajectory in the interval ending
at that time), `observed` (the Kaplan-Meier estimate of survival), and
`median`, `lower` and `upper` of the simulated proportion surviving.
Each group starts with survival 1 at time 0.

## Details

A group is the set of individuals that followed one temperature
trajectory (one `path_id` of the data), or, with `by`, a subset of them
that share a value of another variable, such as a trial. Individuals on
different trajectories are never pooled, since at the same elapsed time
they are at different temperatures.

The observed proportion is the Kaplan-Meier estimate, which allows for
individuals no longer observed; a failure known only to an interval is
placed at the end of the interval. For the simulated proportion, each
posterior draw is used to simulate every individual from the start of
its trajectory to the last time any individual on that trajectory was
observed, with its own covariates and the group-level effects estimated
for its groups. The band therefore includes the variation expected among
samples of this size, which is wide for a group of ten. The simulation
is random: call [`set.seed()`](https://rdrr.io/r/base/Random.html) first
for the same result again.

## See also

[`plot_check()`](https://eogmartins.github.io/rTFT/reference/plot_check.md)
to draw it,
[`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
for residuals by individual.

## Examples

``` r
if (FALSE) { # \dontrun{
chk <- check_survival(fit)
plot_check(chk, time_unit = "min")

# each trial on its own
plot_check(check_survival(fit, by = "trial"))
} # }
```
