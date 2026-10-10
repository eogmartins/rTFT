# Residuals of a fitted model, one per individual

Computes a randomised residual for each individual, in the manner of the
simulation-based residuals of the DHARMa package: the position of what
was observed within what the fitted model predicts for that individual.
When the model describes the data, the residuals are uniform between 0
and 1.

## Usage

``` r
hazard_residuals(fit, n_draws = 400)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- n_draws:

  A single number, the number of posterior draws to average the
  predicted survival over, evenly spaced through the posterior. Defaults
  to 400.

## Value

A tibble with one row per individual: `id`; `path_id`; the variables of
the model that have one value per individual (covariates and groupings);
`time`, the time at which it failed or was last observed; `Ta`, the
temperature at that time; `event`, 1 for a failure and 0 for a survivor;
and `residual`.

## Details

The residual is the fraction of individuals that the model predicts to
have failed by the time this one did, 1 - S(t), with the individual's
own covariates and the group-level effects estimated for its groups,
averaged over the posterior draws. A residual near 0 is an individual
that failed earlier than the model predicts, and one near 1 an
individual that lasted longer. What is observed is not always a single
time, so the residual is drawn at random within what was observed:

- a failure at a known time has the residual 1 - S(t);

- a failure known only to an interval has a residual drawn uniformly
  between that fraction at the start of the interval and at its end;

- a survivor, known only to have outlasted its last observation, has a
  residual drawn uniformly between that fraction at its last observation
  and 1.

The draws are random: call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for the same
residuals again. One residual per individual is used, and not one per
interval, because the intervals of an individual are slices of one
exposure and are not separate observations.

The residuals are to be examined against what was set before the outcome
was known: a covariate, a treatment, a trial. They should not be
examined against the time of failure, or the temperature at failure in a
ramp. An individual that failed late has, by that fact alone, a high
predicted fraction failed at its time of failure, so the residuals rise
with the time of failure under any model, right or wrong.

## References

Hartig, F. (2024). DHARMa: Residual diagnostics for hierarchical
(multi-level / mixed) regression models. R package version 0.4.7.
[doi:10.32614/CRAN.package.DHARMa](https://doi.org/10.32614/CRAN.package.DHARMa)

## See also

[`plot_residuals()`](https://eogmartins.github.io/rTFT/reference/plot_residuals.md)
to draw them,
[`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md)
for a check by group.

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

set.seed(1)
res <- hazard_residuals(fit)
plot_residuals(res)                      # against the uniform distribution
plot_residuals(res, against = "trial")   # by trial
} # }
```
