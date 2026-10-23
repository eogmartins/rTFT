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

The smaller the group, the less the check can see. With fewer than about
five individuals the band runs from 0 to 1 at many of the times, and an
individual falls outside it only if it failed far earlier, or lasted far
longer, than the model predicts; a model that is wrong in the same
direction for every individual can then pass: in a simulation on the
design of the stored exponential fit, with the hazard e times too high
for every individual, the check caught 15% of groups of one and 64% of
groups of five, against every group of 40. The function warns when a
group holds fewer than five individuals, as every group does when each
individual has a temperature trajectory of its own (its own logged
temperatures, in the laboratory or in the field) or when `by` splits the
groups finely. The residuals of
[`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
pool every individual whatever the groups, and are then the check to
judge the model by. In groups of any size, the band is read at every
time of the group, so even under a model that describes the data a group
leaves it at some time now and then.

## See also

[`plot_check()`](https://eogmartins.github.io/rTFT/reference/plot_check.md)
to draw it,
[`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
for residuals by individual.

## Examples

``` r
if (FALSE) { # \dontrun{
# a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
# degrees over 60 minutes and checked every minute, in two trials (see fit_hazard())
path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
                     trial = rep(c("a", "b"), each = 10),
                     t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
                            checks = checks)
fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))

chk <- check_survival(fit)
plot_check(chk, time_unit = "min")

# each trial on its own
plot_check(check_survival(fit, by = "trial"))
} # }
```
