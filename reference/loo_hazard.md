# Leave-one-out cross-validation that leaves out whole individuals or groups

Estimates how well a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
predicts data it was not fitted to, by Pareto-smoothed importance
sampling (PSIS-LOO), leaving out one unit at a time: an individual, or a
group of individuals.

## Usage

``` r
loo_hazard(fit, unit = NULL)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- unit:

  An optional character string, the unit to leave out: `"id"` for the
  individual, or one of the groupings of the model's group-level terms,
  such as `"trial"`, or `"trial:tank"` for the tanks of
  `(1 | trial/tank)`. `NULL` (the default) chooses as described above.

## Value

A `loo` object (see
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)) with one
pointwise value per unit, and an extra element `unit` naming the unit.

## Details

The model is fitted one row per individual and interval at risk, so
[`brms::loo()`](https://mc-stan.org/loo/reference/loo.html) would leave
out one interval, and ask how well the model predicts the rest of an
individual's exposure from its beginning, which is nearly always easy.
Here the pointwise log-likelihood is summed within each unit, which
given the parameters is that unit's log-likelihood exactly, and the unit
is left out whole.

Where individuals share a group-level effect in the model, such as a
random intercept per trial, leaving out one individual still leaves its
trial's effect, estimated from the other individuals of the trial, in
the fit, and the prediction looks better than it would for a new trial.
So by default the unit is the coarsest grouping variable whose levels
individuals share, with a message saying so, and otherwise the
individual. Nested groupings, such as tanks within trials, are handled
by this: leaving out a trial takes its tanks with it. Asking for a finer
unit gives a warning that the coarser grouping's effect stays in the
fit, which is appropriate when the data to be predicted come from the
same levels of it, and optimistic when they come from new ones.

Crossed groupings, such as families spread across trials, have no
default: no single unit removes the effects of both, and leaving out
trials answers how well the model predicts new trials of known families,
while leaving out families answers how well it predicts new families.
The unit must then be given.

Leaving out a large unit can change the posterior too much for PSIS to
stand in for refitting. Where the Pareto k diagnostic of any unit is
above the threshold for the number of draws, a warning says so, and
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
gives an estimate by refitting, without the approximation.

Models compared with
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
must be cross-validated with the same unit.

## References

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing* 27, 1413-1432.

Vehtari, A., Simpson, D., Gelman, A., Yao, Y. and Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research* 25(72), 1-58.

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

fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + (1 | trial)))
loo_hazard(fit)                 # leaves out one trial at a time
} # }
```
