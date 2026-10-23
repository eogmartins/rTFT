# Fit, check and cross-validate the models of a table

Fits each model of a table built with
[`build_hazard_models()`](https://eogmartins.github.io/rTFT/reference/build_hazard_models.md)
to the same data with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md),
checks the sampling of each fit by one rule, and cross-validates each
with
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md),
leaving out the same unit. With `kfold = TRUE` each is also
cross-validated by refitting, with
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
and the same folds for every model. What it returns is what
[`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
and the functions that average models read.

## Usage

``` r
fit_hazard_models(
  models,
  long,
  observation = c("interval", "continuous", "mixed"),
  unit = NULL,
  kfold = FALSE,
  K = NULL,
  rhat = 1.01,
  ess = 400,
  divergent = 0.001,
  ...
)
```

## Arguments

- models:

  A data frame, the table of models from
  [`build_hazard_models()`](https://eogmartins.github.io/rTFT/reference/build_hazard_models.md).

- long:

  A data frame, the table of intervals at risk from
  [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
  to which every model is fitted.

- observation:

  A character string, how the failure times were observed, as in
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- unit:

  An optional character string, the unit to leave out in
  cross-validation, as in
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md).
  `NULL` (the default) chooses it as
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  does.

- kfold:

  A single logical value. `TRUE` also cross-validates each model by
  refitting, which takes one further fit per fold and model. Defaults to
  `FALSE`.

- K:

  An optional single whole number, the number of folds when `kfold` is
  `TRUE`, as in
  [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md).
  `NULL` (the default) gives one fold per unit. With fewer folds than
  units the units are assigned to folds at random, once for all the
  models; call [`set.seed()`](https://rdrr.io/r/base/Random.html) first
  to get the same folds, and the same estimates by refitting, again (see
  [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)).

- rhat:

  A single number above 1, the limit below which the Rhat of every
  parameter must lie. Defaults to 1.01.

- ess:

  A single positive number, the limit above which the bulk and the tail
  effective sample size of every parameter must lie. Defaults to 400.

- divergent:

  A single number from 0 to 1, the largest proportion of the draws that
  may come from divergent transitions. Defaults to 0.001.

- ...:

  Further arguments for every fit, passed to
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  and through it to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
  such as `chains`, `iter`, `cores`, `seed` and `backend`, or the priors
  that
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  takes. Those that
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  takes are given to the refits of `kfold = TRUE` as well.

## Value

A list with four elements. `fits` is a named list of the fitted models
(`brmsfit`), without those that failed. `loo` is a named list of their
cross-validation results from
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md).
`kfold` is a named list of their results from
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md),
or `NULL` when `kfold` is `FALSE`. `summary` is a data frame with one
row per model of the table: its name, failure-rate function and
baseline; `status` (`"fitted"` or `"failed"`) and, for a model that
failed, the reason in `message`; the number of `chains` and `draws`; the
largest `rhat`, the smallest `ess_bulk` and `ess_tail` and the number of
`divergent` draws; `rule_met`, whether the sampling rule is met;
`rhat_group` and `ess_group`, the largest Rhat and the smallest bulk or
tail effective sample size over the group-level effects, outside the
rule (`NA` for a model without them); `elpd_loo` with its standard error
`se_loo`, and `high_k`, the number of units for which importance
sampling cannot be trusted; `elpd_kfold` and `se_kfold` when `kfold` is
`TRUE`; and `seconds`, the time the fit and its cross-validation took.

## Details

The unit left out is settled before any model is fitted, as
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
chooses it, from the group-level terms of the models. Models that are
compared must be cross-validated with the same unit, so the models of
the table must share it; where their group-level terms would give
different units, give `unit`.

A model whose fit fails does not stop the others. It is left out of the
fits returned and is reported, with the reason, in the summary and in a
warning. A fit that lost one of its chains is treated as failed: the
diagnostics of the chains that are left would make it look complete.

The sampling rule is applied to the population-level coefficients and to
the standard deviations and correlations of the group-level effects: the
rank-normalised split-Rhat of each below `rhat`, the bulk and the tail
effective sample size of each above `ess`, and at most the proportion
`divergent` of the draws from divergent transitions. The three check
different things: Rhat whether the chains agree, the effective sample
sizes how precise the draws are, and the divergent transitions whether
the sampler could reach every region; the rule is not called
convergence, which is only the first of them. A fit that does not meet
it is kept, since it is often the one to look at, and is marked in the
summary and named in a warning. Whether to refit it, drop it or keep it
is left to the user; zero divergent transitions and a rule that is met
are not, by themselves, evidence that a model was sampled well. The
limit of 400 on the effective sample sizes is 100 per chain for the four
chains of the default (Vehtari et al. 2021); with more chains, raise it
in proportion. The tolerance on divergent transitions exists so that a
stray one does not fail a fit that is otherwise well sampled; a
divergent transition is still evidence of a region the sampler could not
enter, so where any remain, look at where those draws fall against the
group-level standard deviations before reading the fit, and set
`divergent = 0` to refuse any.

The group-level effects themselves, such as the effect of each trial,
are not part of the rule: a model with many levels would take its worst
case over many parameters and exceed the limits by chance. They are
reported beside it, as the largest Rhat and the smallest effective
sample size over them, because the cross-validation by importance
sampling,
[`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md)
and
[`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
use each trial's fitted effect. The derived quantities do not: they are
computed at the level of the population.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved Rhat
for assessing convergence of MCMC. *Bayesian Analysis* 16, 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

[`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
to weight the models, and
[`average_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/average_tdt_ttl.md),
[`average_dctmax()`](https://eogmartins.github.io/rTFT/reference/average_dctmax.md)
and
[`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md)
to average them.

## Examples

``` r
if (FALSE) { # \dontrun{
# a small ramping trial: 20 individuals warmed from 20 to 35 degrees over 60
# minutes and checked every minute, in two trials (see fit_hazard())
path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
ind    <- data.frame(id = 1:20, path_id = "A", event = 1,
                     trial = rep(c("a", "b"), each = 10),
                     t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
long   <- build_hazard_data(ind, path, "interval", covariates = "trial",
                            checks = checks)

# the exponential failure rate with two baselines, each with a random intercept
# per trial on log a
models <- build_hazard_models("exp", c("mem", "weibull"),
                              dpar_formulas = list(loga = ~ 1 + (1 | trial)))
out <- fit_hazard_models(models, long, chains = 4, iter = 2000, cores = 4, seed = 1)
out$summary

# weight the models by how well each predicts a trial it was not fitted to, and
# average what is derived from them
weights <- hazard_weights(out$loo)
average_tdt_ttl(out$fits, weights, t_ref = 60)$sctmax
} # }
```
