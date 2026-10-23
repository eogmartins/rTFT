# The fitted failure rate over temperature, for each posterior draw

Evaluates the failure rate M(T) of a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
at each temperature in `Ta`, for every posterior draw, at the covariate
values in each row of `newdata`. The rate is computed from the same
expression the model was fitted with, and the failure-rate function,
reference temperature, softening steepness and co-stressors are read
from the fit, so none of them is entered by hand. Random effects are set
to zero, so the rate is that of a typical group.

## Usage

``` r
M_draws(fit, Ta, newdata = NULL, levels = NULL, n_draws = NULL)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- Ta:

  A numeric vector of temperatures at which to evaluate the rate.

- newdata:

  A data frame of covariate values, one row per setting at which to
  evaluate the rate. It must hold every covariate in the model's
  parameter formulas; other columns the model needs are filled in. May
  be `NULL` for a model without covariates.

- levels:

  A named numeric vector giving the level of each co-stressor, on its
  measurement scale, e.g. `c(DO = 6)` for 6 mg/L of dissolved oxygen.
  Required when the model has co-stressors.

- n_draws:

  An optional single number, the number of posterior draws to keep,
  evenly spaced through the posterior. `NULL` (the default) keeps them
  all.

## Value

A tibble with one row per draw, temperature and row of `newdata`: the
columns of `newdata`, `.draw` (the draw's number in the fit, so draws
from different calls can be matched), `Ta` and the failure rate `M`.

## Details

The draws are returned as they are, not summarised, for plotting the
rate with its uncertainty or for further calculation.

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

rate <- M_draws(fit, Ta = seq(25, 40, by = 0.5),
                newdata = data.frame(mass = c(8, 12)))
} # }
```
