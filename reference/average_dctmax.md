# Dynamic critical thermal maximum averaged over several fitted models

Derives the dynamic critical thermal maximum dCTmax under a ramp as
[`derive_dctmax()`](https://eogmartins.github.io/rTFT/reference/derive_dctmax.md)
does, averaged over several models fitted to the same data. Each model's
posterior draws of dCTmax are pooled in proportion to its weight, and
the pooled draws are summarised, so the median and credible interval
reflect both the uncertainty within each model and the uncertainty about
which model predicts best.

## Usage

``` r
average_dctmax(
  fits,
  weights,
  rate,
  T_start,
  newdata = NULL,
  p = 0.5,
  levels = NULL,
  prob = 0.95,
  n_draws = 1000,
  T_max = NULL
)
```

## Arguments

- fits:

  A named list of `brmsfit`s, models fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  to the same data.

- weights:

  The weights from
  [`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md),
  whose model names are names of `fits`. Only the models it kept are
  averaged.

- rate:

  A numeric vector of ramp rates, in degrees per unit of time of the
  data the model was fitted to (e.g. degrees C per minute), all
  positive.

- T_start:

  A numeric vector of starting temperatures, one per rate or a single
  one for all.

- newdata:

  A data frame of covariate values, one row per setting at which to
  derive. It must hold every covariate in the model's parameter
  formulas. May be `NULL` for a model without covariates.

- p:

  A numeric vector of the fractions of individuals failing at which
  sCTmax and the thermal death time curves are derived, each strictly
  between 0 and 1.

- levels:

  A named numeric vector giving the level of each co-stressor, on its
  measurement scale, e.g. `c(DO = 6)`. Required when the model has
  co-stressors.

- prob:

  A single number strictly between 0 and 1, the coverage of the credible
  intervals.

- n_draws:

  A single number, the number of posterior draws in the pool. Each model
  contributes its share of them, evenly spaced through its posterior.
  Where a model holds fewer draws than its share, the pool is made
  smaller so that the shares stay in proportion to the weights, and a
  message says so.

- T_max:

  A single number, the temperature up to which the ramp is followed.
  `NULL` (the default) takes 15 degrees above the warmest temperature
  the models were fitted to.

## Value

A list of two tibbles: `dctmax`, in the layout of
[`derive_dctmax()`](https://eogmartins.github.io/rTFT/reference/derive_dctmax.md)'s
result, and `weights`, the table of weights with the number of draws
each model contributed.

## Details

The models usually agree closely on dCTmax under a ramp like those they
were fitted to, and less under ramps much slower or faster, where the
failure-rate and baseline functions are extrapolated.

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
fit_mem     <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
fit_weibull <- fit_hazard(long, "exp", baseline = "weibull",
                          dpar_formulas = list(loga = ~ 1 + mass))
fits        <- list(memoryless = fit_mem, weibull = fit_weibull)

w    <- hazard_weights(lapply(fits, loo_hazard))
average_dctmax(fits, w, rate = c(0.1, 0.3), T_start = 14,
               newdata = data.frame(mass = 10))$dctmax
} # }
```
