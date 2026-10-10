# Survival along a record of temperature, averaged over several fitted models

Predicts survival along a record as
[`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
does, averaged over several models fitted to the same data. Each model's
posterior draws are pooled in proportion to its weight and the pool is
summarised, so the credible intervals reflect both the uncertainty
within each model and the uncertainty about which model holds. This
matters most where the record lies outside the temperatures of the
assay, where failure-rate functions that fit the data equally well can
predict very differently (see
[`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)).

## Usage

``` r
average_survival(
  fits,
  weights,
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

- fits:

  A named list of `brmsfit`s, models fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  to the same data.

- weights:

  The weights from
  [`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md),
  whose model names are names of `fits`. Only the models it kept are
  averaged.

- record:

  A data frame with the columns `time` and `temp`, the record of
  temperature through time, in the time unit and temperature scale of
  the data the model was fitted to, with times increasing and no gaps.
  It may also hold the column of a co-stressor and the column named in
  `by`.

- newdata:

  A data frame of covariate values, one row per setting at which to
  derive. It must hold every covariate in the model's parameter
  formulas. May be `NULL` for a model without covariates.

- levels:

  A named numeric vector giving the level of each co-stressor, on its
  measurement scale, e.g. `c(DO = 6)`. Required when the model has
  co-stressors.

- by:

  A character string, the name of a column of `record` giving the period
  each reading belongs to (e.g., a day). When given, the result has each
  period's share of the hazard accrued over the record. An interval
  between two readings belongs to the period of the first.

- prob:

  A single number strictly between 0 and 1, the coverage of the credible
  intervals.

- n_draws:

  A single number, the number of posterior draws in the pool. Each model
  contributes its share of them, evenly spaced through its posterior.
  Where a model holds fewer draws than its share, the pool is made
  smaller so that the shares stay in proportion to the weights, and a
  message says so.

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

The list
[`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
returns, with `weights`, the table of weights with the number of draws
each model contributed. In `coverage`, the shares of the risk are those
of the pooled draws.

## Details

The weights say how well each model predicts failures within the assay.
They do not say which model extrapolates better below it, which the data
cannot tell.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- list("exp" = fit_exp, "zeroexp" = fit_zeroexp)
w    <- hazard_weights(lapply(fits, loo_hazard))
average_survival(fits, w, record, newdata = data.frame(mass = 10))$mortality
} # }
```
