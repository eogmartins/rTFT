# Cross-validation by refitting, leaving out whole individuals or groups

Estimates how well a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
predicts data it was not fitted to, exactly, by refitting it with each
unit left out in turn and predicting the unit from the refit. It answers
the same question as
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md),
with the same choice of unit, but without its approximation: use it
where
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
warns that PSIS cannot be trusted, which is typical when whole trials
with a strong trial effect are left out.

## Usage

``` r
kfold_hazard(fit, unit = NULL, K = NULL, ...)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- unit:

  An optional character string, the unit to leave out, as in
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md).

- K:

  An optional single whole number, the number of folds, from 2 to the
  number of units. `NULL` (the default) gives one fold per unit, leaving
  out each unit on its own.

- ...:

  Further arguments for the refits, passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  through [`stats::update()`](https://rdrr.io/r/stats/update.html), such
  as `chains`, `iter` or `cores`.

## Value

A `kfold` object in the layout of
[`loo::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html),
with one pointwise value per unit, which
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
accepts, and extra elements `unit` (its name), `K`, and `folds` (the
fold of each unit).

## Details

The units are split into `K` folds, and the model is refitted once per
fold without that fold's units. A unit's rows are summed within each
draw of the refit, so its value is its joint log predictive density,
averaged over the draws. Where the unit's level of a group-level term
was left out entirely (a whole trial), its effect is drawn from the
fitted population distribution, so the trial is predicted as a new trial
would be.

The refits take time: one model fit per fold. With `K` smaller than the
number of units, the units are assigned to folds at random; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first to get the same
folds again.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + (1 | trial)))
cv  <- kfold_hazard(fit, cores = 4)       # one refit per trial
} # }
```
