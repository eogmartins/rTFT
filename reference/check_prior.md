# Check whether the priors of a fitted model pull its estimates

Measures how much each estimate of a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
depends on its priors, block by block, by power-scaling (Kallioinen et
al. 2024): each block of priors is raised to a power slightly above and
below 1, the posterior is reweighted to match, and the change in each
estimate's distribution is measured. An estimate that moves when one
block is scaled is pulled by that block. A narrow prior can pull an
estimate with no sign in the convergence diagnostics. That is most
likely in a design whose groups fail far from the reference temperature
or far from each other. log a is the log failure rate at the reference
temperature, so a group that fails some degrees away has its log a
extrapolated from its failures over that distance, a change of about b
times the distance. A prior narrow beside that change draws log a back
towards the prior's centre, and the group's other parameters move to fit
its failures.

## Usage

``` r
check_prior(fit, newdata = NULL, t_ref = NULL, p = 0.5, levels = NULL)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- newdata:

  A data frame of covariate values, one row per setting at which to
  check the reported parameters (and sCTmax, with `t_ref`), e.g. one row
  per species. It must hold every covariate in the model's parameter
  formulas, as in
  [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md),
  and the quantities are computed with the random effects at zero, as
  there. `NULL` (the default) checks the coefficients only.

- t_ref:

  A single positive number, the reference exposure time for sCTmax, in
  the time unit of the data the model was fitted to. `NULL` (the
  default) leaves sCTmax out; it needs `newdata`.

- p:

  A single number strictly between 0 and 1, the fraction failing within
  `t_ref` that defines sCTmax.

- levels:

  A named numeric vector giving the level of each co-stressor, on its
  measurement scale, at which the parameters and sCTmax are checked;
  required with `newdata` when the model has co-stressors, as in
  [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md),
  and ignored, with a warning, without `newdata`.

## Value

A list of two tibbles:

- `sensitivity`, one row per block and quantity: `block`; `quantity`, a
  coefficient as brms names its draws (`b_loga_Intercept`) or, for each
  row of `newdata`, a reported parameter (`a`, `b`, `shape`, ...) or
  `sCTmax`; the columns of `newdata`, NA for a coefficient; `prior` and
  `likelihood`, the sensitivity to scaling the block and to scaling the
  likelihood; `diagnosis`, priorsense's reading of the two, in its own
  words; and `k_prior` and `k_likelihood`, the Pareto k of the scaling
  of the block and of the likelihood (see *Reading the sensitivity*).

- `priors`, one row per population-level coefficient: `parameter`,
  `coefficient`, `prior`; `prior_centre` and `prior_scale` (the mean and
  standard deviation of a normal, the location and scale of a Student-t,
  NA for a flat prior); the posterior `median` and `sd`; `distance`, the
  median's distance from the prior's centre in units of its scale; and
  `sd_ratio`, the posterior standard deviation over the prior's scale.

## Details

The model is not refitted. The check needs the log prior of each draw,
which brms keeps, and the log-likelihood of each draw, which it computes
once over all draws (about half a minute for 4000 draws of a model
fitted to 14,000 rows).

## Blocks of priors

The priors are scaled in blocks: one per parameter of the model (`loga`,
`logb`, `Tc`, the baseline's parameters), holding its intercept and its
slopes; one per co-stressor, named `c<name>` after its main effect,
holding the main effect and the interaction with temperature; `sd`,
holding every prior on a group-level term (standard deviations and
correlations); and `all`, every prior together. The priors of the
population-level coefficients are read from the fit and evaluated at
each draw, and `sd` is the rest of the log prior brms kept. Without
group-level terms that rest must be zero, which checks that the priors
were read as brms used them. Normal and Student-t priors (and flat ones)
on the coefficients can be read; any other is refused, as are bounds on
the coefficients and terms other than linear and group-level ones
(splines, monotonic effects). A block whose priors are all flat has
nothing to scale and is left out, with a message.

## Reading the sensitivity

priorsense reads the two sensitivities together, against a threshold of
0.05. When both are at or above it, the prior and the data disagree, and
it reports a "potential prior-likelihood conflict". When only the
sensitivity to the prior is above it, the prior rather than the data
sets the estimate, and it reports a "potential strong prior / weak
likelihood"; every other row reads "-". In a model with group-level
terms, scaling the likelihood is like changing the number of individuals
within each group, not the number of groups, so a quantity that the
groups inform, such as their standard deviation, can show little
likelihood sensitivity without the likelihood being weak (Kallioinen et
al. 2024, section 2.5). The threshold marks a sensitivity that can be
detected, not one that matters. In simulations of a design with three
species whose failures lay some degrees apart, a sensitivity of the
Weibull shape of 0.18 to 0.28 to the block `loga` came with a bias in
the shape of about two and a half posterior standard deviations, while
sensitivities up to 0.09 came with biases under a third of one; no fit
fell between 0.09 and 0.18. Nor is the size of a sensitivity the size of
a bias: in the same simulations sCTmax reached sensitivities of 0.06 to
0.12 while its bias was at most 0.08 degrees. The other blocks never
exceeded 0.03.

Each scaling is an importance-sampling estimate, and the Pareto k of its
weights says whether they can be trusted: `k_prior` is the larger k of
the two scalings of the block and `k_likelihood` that of the
likelihood's two. A k above priorsense's threshold (min(1 - 1 /
log10(S), 0.7) for S draws; 0.7 from 4,000 draws) warns, and the
sensitivities of that scaling are a sign at most. A block whose prior
lies far from the posterior is where the weights degrade, and that is
the case the check is for, so the warning is to be read rather than
silenced.

The test of a flagged block is a refit with that block's priors widened
(through `intercept_priors` and `slope_prior` of
[`build_hazard_prior()`](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)),
compared with the first fit. The table `priors` shows which coefficients
lie far out in their prior or have a posterior almost as wide as it.

## References

Kallioinen, N., Paananen, T., Burkner, P.-C. and Vehtari, A. (2024).
Detecting and diagnosing prior and likelihood sensitivity with
power-scaling. *Statistics and Computing* 34, 57.

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

chk <- check_prior(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
chk$sensitivity
chk$priors
} # }
```
