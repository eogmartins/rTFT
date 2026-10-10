# Assemble the brms formula and family for a failure-rate hazard model

Builds a nonlinear brms formula whose linear predictor is the log
failure rate log M(T), plus (unless the baseline is memoryless) the log
baseline hazard log h0(tau), plus the log exposure of each at-risk
interval. The failure rate and baseline are named by their keys in
[M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md)
and
[h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md),
and their log-scale expression strings
([brms_log_M](https://eogmartins.github.io/rTFT/reference/brms_log_M.md),
[brms_log_h0](https://eogmartins.github.io/rTFT/reference/brms_log_h0.md))
are pasted straight into the predictor, so the formula fitted here is
the same expression
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
assumes and that reporting later evaluates.

## Usage

``` r
build_hazard_spec(
  m_form,
  baseline = "mem",
  dpar_formulas = list(),
  observation = c("interval", "continuous", "mixed"),
  k_sw = NULL,
  costressors = list()
)
```

## Arguments

- m_form:

  A character string denoting a key of
  [M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md)
  (e.g. `"exp"` or `"exp_threshold"`).

- baseline:

  A character string denoting a key of
  [h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md):
  `"mem"` (default), `"weibull"`, `"llogis"` or `"gompertz"`.

- dpar_formulas:

  A named list of one-sided formulas, one per parameter you want
  covariates on. A parameter the functions need but you omit defaults to
  `~ 1` (intercept only); a name not among the parameters is ignored
  with a warning. For example,
  `list(loga = ~ 1 + acclim + (1 | trial), logb = ~ 1 + (1 | trial))`.

- observation:

  A character string denoting how each failure time was observed, as in
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md):
  `"interval"` (the
  [cloglog_bern_family](https://eogmartins.github.io/rTFT/reference/cloglog_bern_family.md)),
  `"continuous"` (a Poisson piecewise-exponential), or `"mixed"` (the
  [hazard_mixed_family](https://eogmartins.github.io/rTFT/reference/hazard_mixed_family.md),
  which reads the per-row `exact` flag through `vint()`).

- k_sw:

  An optional positive numeric representing the softening steepness of a
  threshold function's gate, in inverse temperature units (the
  transition band is roughly its reciprocal). This is an advanced
  override: leave `NULL` (the default) to use the value in
  [M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md)
  (50), which suits most fits. A larger value is a sharper switch over a
  narrower band, but too large returns the sampler to the hard step
  geometry, where the critical temperature `Tc` mixes poorly (low
  effective sample size, sometimes divergences). This is a
  sampling-geometry effect, separate from identifiability: `Tc` is
  poorly identified whenever the observed temperatures do not bracket
  it, whatever `k_sw` is. Ignored, with a warning, for a function that
  has no threshold.

- costressors:

  An optional list of
  [costressor](https://eogmartins.github.io/rTFT/reference/costressor.md)
  objects. Each appends its terms to log M(T) (a main effect `c<name>`,
  and unless switched off a temperature interaction `dT<name>`), so the
  failure rate becomes a surface M(T, Z). Their coefficients join the
  parameters getting a sub-formula, so a `dpar_formulas` entry can put
  covariates on one. Defaults to none.

## Value

A named list with the brms `formula`, the `family`, and the
`observation` mode.

## Details

Each parameter the chosen functions need (`loga`, `logb`, `Tc`, and so
on) gets its own sub-formula, an intercept only (`~ 1`) by default. Pass
a one-sided formula in `dpar_formulas` to put covariates on a parameter;
a formula named for a parameter the chosen functions do not use is
ignored with a warning.

## Examples

``` r
# An exponential failure rate and memoryless baseline, with interval-censored
# observation (the default).
build_hazard_spec("exp")
#> $formula
#> failed ~ loga + exp(logb) * (Ta - Tr) + log_exposure 
#> loga ~ 1
#> logb ~ 1
#> 
#> $family
#> 
#> Custom family: cloglog_bern 
#> Link function: identity 
#> Parameters: mu 
#> 
#> 
#> $observation
#> [1] "interval"
#> 

# An exponential-threshold failure rate with a Weibull baseline hazard in elapsed
# exposure time.
build_hazard_spec("exp_threshold", baseline = "weibull")
#> $formula
#> failed ~ logshape - log(tau_ref) + (exp(logshape) - 1) * (log_tau - log(tau_ref)) + loga + exp(logb) * (Ta - Tc) - log1p(exp(-50 * (Ta - Tc))) + log_exposure 
#> loga ~ 1
#> logb ~ 1
#> Tc ~ 1
#> logshape ~ 1
#> 
#> $family
#> 
#> Custom family: cloglog_bern 
#> Link function: identity 
#> Parameters: mu 
#> 
#> 
#> $observation
#> [1] "interval"
#> 

# A different structure per parameter: acclim as a fixed covariate on the
# amplitude loga, and a trial-level random intercept on both loga and the
# temperature sensitivity logb (`acclim` and `trial` are carried through by
# build_hazard_data()).
build_hazard_spec("exp",
  dpar_formulas = list(loga = ~ 1 + acclim + (1 | trial), logb = ~ 1 + (1 | trial)))
#> $formula
#> failed ~ loga + exp(logb) * (Ta - Tr) + log_exposure 
#> loga ~ 1 + acclim + (1 | trial)
#> logb ~ 1 + (1 | trial)
#> 
#> $family
#> 
#> Custom family: cloglog_bern 
#> Link function: identity 
#> Parameters: mu 
#> 
#> 
#> $observation
#> [1] "interval"
#> 

# Mixed interval-censored and known-time failures: the response holds the per-row
# `exact` flag through vint().
build_hazard_spec("exp", observation = "mixed")
#> $formula
#> failed | vint(exact) ~ loga + exp(logb) * (Ta - Tr) + log_exposure 
#> loga ~ 1
#> logb ~ 1
#> 
#> $family
#> 
#> Custom family: hazard_mixed 
#> Link function: identity 
#> Parameters: mu 
#> 
#> 
#> $observation
#> [1] "mixed"
#> 

# A crossed temperature-by-oxygen design: dissolved oxygen enters log M as a main
# effect (cDO) and a temperature interaction (dTDO).
do <- costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
                 interaction = TRUE)
build_hazard_spec("exp", costressors = list(do))
#> $formula
#> failed ~ loga + exp(logb) * (Ta - Tr) + cDO * z_DO + dTDO * (Ta - Tr) * z_DO + log_exposure 
#> loga ~ 1
#> logb ~ 1
#> cDO ~ 1
#> dTDO ~ 1
#> 
#> $family
#> 
#> Custom family: cloglog_bern 
#> Link function: identity 
#> Parameters: mu 
#> 
#> 
#> $observation
#> [1] "interval"
#> 

# Advanced: a gentler threshold switch than the catalogue default.
build_hazard_spec("exp_threshold", k_sw = 20)
#> $formula
#> failed ~ loga + exp(logb) * (Ta - Tc) - log1p(exp(-20 * (Ta - Tc))) + log_exposure 
#> loga ~ 1
#> logb ~ 1
#> Tc ~ 1
#> 
#> $family
#> 
#> Custom family: cloglog_bern 
#> Link function: identity 
#> Parameters: mu 
#> 
#> 
#> $observation
#> [1] "interval"
#> 
```
