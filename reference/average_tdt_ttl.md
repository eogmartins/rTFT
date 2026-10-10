# Thermal tolerance averaged over several fitted models

Derives the static critical thermal maximum sCTmax, the thermal death
time curves and the thermal tolerance landscape as
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
does, averaged over several models fitted to the same data. Each model's
posterior draws of a quantity are pooled in proportion to its weight,
and the pooled draws are summarised, so the median and credible interval
reflect both the uncertainty within each model and the uncertainty about
which model holds.

## Usage

``` r
average_tdt_ttl(
  fits,
  weights,
  newdata = NULL,
  p = c(0.1, 0.5, 0.9),
  t_ref,
  Ta = NULL,
  tau = NULL,
  levels = NULL,
  prob = 0.95,
  n_draws = 1000
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

- newdata:

  A data frame of covariate values, one row per setting at which to
  derive. It must hold every covariate in the model's parameter
  formulas. May be `NULL` for a model without covariates.

- p:

  A numeric vector of the fractions of individuals failing at which
  sCTmax and the thermal death time curves are derived, each strictly
  between 0 and 1.

- t_ref:

  A single positive number, the reference exposure time for sCTmax, in
  the time unit of the data the model was fitted to (for example `60`
  for one hour when times are in minutes). It has no default, because
  the fit does not record that unit.

- Ta:

  An optional numeric vector of the constant temperatures for the
  thermal death time curves and the tolerance landscape. `NULL` (the
  default) takes 50 evenly spaced across the temperatures the data
  cover.

- tau:

  An optional numeric vector of the elapsed exposure times for the
  tolerance landscape, in the time unit of the data, none negative.
  `NULL` (the default) takes 50 spaced evenly on a log scale across the
  times the data cover. A time of zero is allowed, and gives a survival
  of 1, so that a survival curve can start there.

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

## Value

A list of four tibbles: `sctmax`, `tdt_curves` and `ttl_surface` as from
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md),
and `weights`, the table of weights with the number of draws each model
contributed.

## Details

The parameters are not averaged. The parameters of different
failure-rate functions are different quantities under the same names, so
a pooled a or b would mean nothing; they are reported for each model by
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md).

The result has the layout of
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)'s,
without its `params` table, so
[`plot_tdt()`](https://eogmartins.github.io/rTFT/reference/plot_tdt.md),
[`plot_ttl()`](https://eogmartins.github.io/rTFT/reference/plot_ttl.md)
and
[`plot_ttl_3d()`](https://eogmartins.github.io/rTFT/reference/plot_ttl_3d.md)
draw it as they draw a single model's.

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
avg  <- average_tdt_ttl(fits, w, newdata = data.frame(mass = 10), t_ref = 60)
plot_tdt(avg)
} # }
```
