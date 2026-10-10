# Weakly-informative default priors for a failure-rate hazard model

Supplies a sensible starting prior for the model
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md)
assembles, keyed to each parameter's intercept and unaffected by
whatever covariates the `dpar_formulas` add. For every parameter the
chosen functions need, it sets a weakly-informative prior on the
intercept (the parameter's value at the reference temperature and start
of exposure), a generic Normal(0, 1) on any covariate slopes, and an
Exponential(2) on the standard deviation of any random effect. The slope
and random-effect priors are added only when the parameter actually has
those terms, so brms raises no "no parameter" warning. Override the
centre of any parameter through `intercept_priors`, or the shared slope
and random-effect standard-deviation priors through `slope_prior` and
`sd_prior`; or replace the whole set with your own `prior` when fitting
with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

## Usage

``` r
build_hazard_prior(
  m_form,
  baseline = "mem",
  dpar_formulas = list(),
  intercept_priors = list(),
  slope_prior = "normal(0, 1)",
  sd_prior = "exponential(2)",
  costressors = list(),
  long
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

  A named list of one-sided formulas, one per parameter, as in
  [build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md).
  A parameter's fixed and random terms decide which slope and
  standard-deviation priors are added. A formula without an intercept
  (`~ 0 + acclim`) gives the parameter one coefficient per level in
  place of the intercept, and each of them takes the intercept's prior.

- intercept_priors:

  A named list of prior strings, one per parameter, that replace the
  default intercept (centre) for those parameters, e.g.
  `list(Tc = "normal(26, 4)")`. A parameter left out keeps its default;
  a name that is not a parameter of the chosen functions is ignored with
  a warning.

- slope_prior:

  A single prior string for every covariate slope (a population-level
  effect beyond the intercept), applied to each parameter that has one.
  Defaults to `"normal(0, 1)"`.

- sd_prior:

  A single prior string for every random-effect standard deviation,
  applied to each parameter that has one. Defaults to
  `"exponential(2)"`.

- costressors:

  An optional list of
  [costressor](https://eogmartins.github.io/rTFT/reference/costressor.md)
  objects, as in
  [build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md).
  Each adds an intercept prior for its main effect and, when the
  temperature interaction is kept, its interaction coefficient. Defaults
  to none.

- long:

  A data frame, the row-per-interval table from
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  that the prior is for. The defaults read from it the reference time
  `tau_ref`, the reference temperature `Tr` and the temperatures at
  which failures were observed. It has no default.
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  passes its own table when it builds the default prior.

## Value

A `brmsprior`.

## Details

The defaults take three things from the data in `long`, so that they say
the same thing whatever the organism, the time unit and the length of
the exposures: the reference time `tau_ref`, the reference temperature
`Tr`, and the median temperature at which failures were observed, T_f.
They are weakly-informative and meant to be changed to suit your study:

- `loga` – log amplitude. Stated where the failures are: at T_f, the
  cumulative hazard reached over `tau_ref` is centred at exp(-1), about
  0.37, with a standard deviation of 2 on the log scale (about 0.007 to
  20). With the memoryless and Gompertz baselines a is a rate per unit
  time, so that centre is -1 - log(tau_ref); with the Weibull and
  log-logistic it is -1. The prior is then carried from T_f back to the
  temperature a refers to (see below).

- `logb ~ normal(-0.7, 1)` – log temperature sensitivity

- `logc` – log baseline term of the power function (see below)

- `Tc ~ normal(T_f, 3)` – critical temperature, centred where failures
  were observed

- `logshape ~ normal(0, 1)` – log Weibull / log-logistic shape

- `logscale ~ normal(log(tau_ref), 1)` – log log-logistic scale

- `gam ~ normal(0, 3 / tau_ref)` – Gompertz rate: the baseline hazard
  changes by a factor of exp(gamma tau_ref) across the reference time,
  within about 400-fold either way with 95% probability. It is centred
  at zero, the memoryless case, so unless the data show a trend the fit
  stays near a flat baseline.

with `normal(0, 1)` on any covariate slope and `exponential(2)` on any
random-effect standard deviation. Each co-stressor in `costressors` adds
a `normal(0, 1)` on its main effect `c<name>` and a tighter
`normal(0, 0.5)` on its temperature interaction `dT<name>`, both
overridable through `intercept_priors`.

The constants in these priors (the centre of exp(-1) and the standard
deviation of 2 for the amplitude, the centre of -0.7 for log b, the rate
of 2 for the group standard deviations, the 0.5 for an interaction) are
judgements of what is plausible for a thermal failure assay, not values
estimated from data, and the centres that are read from `long` make the
prior depend on the data in a weak way, as a standard deviation of 2 on
the log scale is wide. The prior on a covariate slope is per unit of the
covariate, so a covariate measured in large units (body mass in grams,
say) should be centred and scaled, or the slope prior widened, before
`normal(0, 1)` is taken as weak. With few groups, the prior on a group
standard deviation can matter for predictions of new groups (which
cross-validation by group scores, see
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md));
a sensitivity check against `exponential(1)` is cheap insurance.

## The prior on a when the reference is far from the failures

For the rates anchored at the reference temperature (`exp`, `zeroexp`,
`zeropower`, `power`), a is the rate at `Tr`, or its scale there, and
`Tr` may lie a distance d = T_f - Tr from where failures were observed.
The log rate at T_f is then loga plus a term in b and d (b d for `exp`),
so loga is the log rate at T_f less that term, and it inherits the
uncertainty in b across the distance:

- its centre is the centre at T_f less the term at the prior median of
  b;

- its variance is 2^2 plus the square of (the change in the term per
  unit of b, times the prior standard deviation of b).

For `exp` this is `normal(centre - 0.5 d, sqrt(4 + (1.07 d)^2))`: the
usual prior when the reference is at the failures (d = 0), and a much
wider one 14 degrees below them. A prior of standard deviation 2 on loga
at a distant reference would in effect be a tight prior on b, and biases
it. The threshold functions are anchored at `Tc`, not `Tr`, and take the
centre at T_f unchanged. The calculation uses the default prior on
`logb`; if you change that prior and keep a distant reference, set
`loga` to match.

## The baseline term c of the power function

In M(T) = c + a (T - Tr)^b, c is the rate that remains at the reference
temperature: failure that happens whatever the warming. Its default is
centred so that about 1% of individuals would fail over `tau_ref` for
reasons unrelated to temperature, `normal(log(0.01) - log(tau_ref), 2)`
with the memoryless or Gompertz baseline and `normal(log(0.01), 2)`
otherwise, which covers background failure from about 0.02% to 40%. The
1% is a judgement, not a measured value. Where controls were held at the
reference temperature, set it from them: if a fraction q of controls
failed over a time t, c is about -log(1 - q) / t, so for instance
`intercept_priors = list(logc = "normal(log(-log(1 - 0.02) / 120), 1)")`
for 1 in 50 failing over 120 minutes.

## Examples

``` r
# A ramp from 20 to 35 degrees over 60 minutes, with two failures.
path   <- ramp_path("A", T_start = 20, rate = 0.25, t_max = 60)
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))
ind    <- data.frame(id = 1:3, path_id = "A", event = c(1, 1, 0),
                     t_upper = c(30, 40, 60), acclim = c(10, 15, 15))
long   <- build_hazard_data(ind, path, "interval", covariates = "acclim",
                            checks = checks)
#> ℹ The intervals at risk were divided further, into 277 rows from 13, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
#> Warning: ! For 2 individuals, temperature changes by up to 2.5 degrees within the
#>   interval in which failure was found.
#> ℹ That interval cannot be divided, since where in it the individual failed is
#>   not known, and the hazard is taken at its midpoint, so the estimates may be
#>   somewhat biased.

# Default weakly-informative priors for an exponential-threshold + Weibull model.
build_hazard_prior("exp_threshold", baseline = "weibull", long = long)
#>            prior class      coef group resp dpar    nlpar   lb   ub tag source
#>    normal(-1, 2)     b Intercept                     loga <NA> <NA>       user
#>  normal(-0.7, 1)     b Intercept                     logb <NA> <NA>       user
#>  normal(27.5, 3)     b Intercept                       Tc <NA> <NA>       user
#>     normal(0, 1)     b Intercept                 logshape <NA> <NA>       user

# Override a couple of priors (a cooler species, a tighter amplitude) and narrow
# the covariate-slope prior; every other parameter keeps its default.
build_hazard_prior("exp_threshold", baseline = "weibull",
  dpar_formulas = list(loga = ~ 1 + acclim),
  intercept_priors = list(Tc = "normal(26, 4)", loga = "normal(-6, 1.5)"),
  slope_prior = "normal(0, 0.5)", long = long)
#>            prior class      coef group resp dpar    nlpar   lb   ub tag source
#>  normal(-6, 1.5)     b Intercept                     loga <NA> <NA>       user
#>   normal(0, 0.5)     b                               loga <NA> <NA>       user
#>  normal(-0.7, 1)     b Intercept                     logb <NA> <NA>       user
#>    normal(26, 4)     b Intercept                       Tc <NA> <NA>       user
#>     normal(0, 1)     b Intercept                 logshape <NA> <NA>       user
```
