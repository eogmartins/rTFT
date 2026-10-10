# Weights for averaging failure-rate models

Weights several models fitted to the same data by how well each predicts
data it was not fitted to, from their cross-validation results. The
models can be any set: any failure-rate functions with any baselines.
Which models are averaged is chosen by which results are given.

## Usage

``` r
hazard_weights(cv, method = c("stacking", "pseudobma"), min_weight = 1)
```

## Arguments

- cv:

  A named list of cross-validation results, one per model, from
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  or
  [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md).
  All must have the same unit and the same units left out. The names
  identify the models.

- method:

  A character string, `"stacking"` (the default) or `"pseudobma"`.

- min_weight:

  A single number above 0 and at most 1, the cumulative weight the
  models kept must reach. The default, 1, keeps every model.

## Value

A tibble with one row per model, highest weight first: `model`; `elpd`,
its summed log predictive density; `weight_all`, its weight among all
the models, and `cumulative`, the running total; `kept`, whether it is
kept; and `weight`, its weight among the models kept, which is the
weight to average with (0 for a model not kept). The unit, the number of
units and the method are attached as attributes.

## Details

With `method = "stacking"` the weights are those of the combination of
models that predicts best. It is the method Yao et al. (2018) recommend,
and the default here. Stacking can give a model a weight of zero. This
happens when two models make nearly the same predictions: including both
adds nothing to the combination, so one takes the weight and the other
gets none, even though it would predict almost as well on its own. With
`method = "pseudobma"` each model is weighted by its own predictive
density, allowing for its uncertainty over the units by the Bayesian
bootstrap; these weights are smoother, and steadier when the units are
few. The Bayesian bootstrap is random: call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for the same
weights again.

With `min_weight` below 1, only the highest-weighted models are kept, as
many as it takes for their weights to add up to at least `min_weight`,
and the weights of those kept are then computed again among themselves.
Recomputing matters for stacking, whose weights are found jointly for a
set of models and are not simply rescaled when a model leaves the set.

The weights rest on one value per unit left out, so they are less
certain the fewer units there are, as when the unit is one of a handful
of trials. A message says how many units there were.

Results from
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
and
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
can be mixed, as when the leave-one-out approximation is unreliable for
one model and that model alone is cross-validated by refitting. The two
agree in what they estimate when each fold leaves out a single unit.
With fewer folds than units, each refit sees less of the data than a
leave-one-out fit does, so its predictions score slightly lower and its
model is at a small disadvantage. A message names the models this
applies to.

## References

Yao, Y., Vehtari, A., Simpson, D. and Gelman, A. (2018). Using stacking
to average Bayesian predictive distributions. *Bayesian Analysis* 13,
917-1007.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- list("exp" = fit_exp, "zeroexp" = fit_zeroexp)
cv   <- lapply(fits, loo_hazard)
hazard_weights(cv)
hazard_weights(cv, min_weight = 0.95)
} # }
```
