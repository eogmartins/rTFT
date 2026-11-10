# Weakly-informative default priors for a failure-rate hazard model

Supplies a sensible starting prior for the model
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md)
assembles, keyed to each parameter's intercept and unaffected by
whatever covariates the `dpar_formulas` add. For every parameter the
chosen functions need, it sets a weakly-informative prior on the
intercept (the parameter's value at the reference temperature and start
of exposure), a generic Normal(0, 1) on any covariate slopes, and an
Exponential(2) on the standard deviation of any random effect; where the
formula of log a has terms that mark groups, its priors are scaled to
the groups instead (see *Groups whose failures lie far apart*). The
slope and random-effect priors are added only when the parameter
actually has those terms, so brms raises no "no parameter" warning.
Override the centre of any parameter through `intercept_priors`, the
slope priors through `slope_prior` (one for every parameter, or one per
parameter) and the random-effect standard-deviation priors through
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
  effect beyond the intercept), applied to each parameter that has one;
  or a named list of prior strings, one per parameter, e.g.
  `list(loga = "normal(0, 5)")`, the parameters left out taking
  `"normal(0, 1)"`. A name that is not a parameter of the chosen
  functions is ignored with a warning, as is one whose formula has no
  slope for it to reach (no covariate, or no intercept, when each
  coefficient takes the intercept prior). The group terms of log a take
  the wider of this and the prior scaled to the groups (see *Groups
  whose failures lie far apart*). Defaults to `"normal(0, 1)"`.

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

with `normal(0, 1)` on any covariate slope (except the group terms of
log a, see below) and `exponential(2)` on any random-effect standard
deviation. Each co-stressor in `costressors` adds a `normal(0, 1)` on
its main effect `c<name>` and a tighter `normal(0, 0.5)` on its
temperature interaction `dT<name>`, both overridable through
`intercept_priors`.

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
The log rate at T_f is then log a plus a term in b and d (b d for
`exp`), so log a is the log rate at T_f less that term, and it inherits
the uncertainty in b across the distance:

- its centre is the centre at T_f less the term at the prior median of
  b;

- its variance is 2^2 plus the square of (the change in the term per
  unit of b, times the prior standard deviation of b).

For `exp` with the default prior on `logb` this is
`normal(centre - 0.5 d, sqrt(4 + (1.07 d)^2))`: the usual prior when the
reference is at the failures (d = 0), and a much wider one 14 degrees
below them. A prior of standard deviation 2 on log a at a distant
reference would in effect be a tight prior on b, and biases it. The
calculation uses the prior on `logb` in force: the one given in
`intercept_priors` when it is a normal, normal(mu, s), whose b has
median exp(mu) and standard deviation sqrt((exp(s^2) - 1) exp(2 mu +
s^2)); the default for any other.

The threshold functions are anchored at `Tc`, not `Tr`, so for them d is
the distance from the centre of the prior on `Tc` in force to the
failures. The default prior on `Tc` is centred at T_f, and so the centre
of log a stays at T_f; so it does for a prior on `Tc` centred above the
failures, where the rate is held near zero. However, `Tc` may lie
anywhere in its prior, and well above `Tc` the data fix only log a less
the function's term at the distance from `Tc` to the failures (log a - b
`Tc` for `exp_threshold`). A prior of standard deviation 2 on log a
would then hold `Tc` close to the failures with a confidence the data do
not give: in simulations with `Tc` 4 degrees below the failures it
placed `Tc` about 1.2 degrees too high, with intervals that never
covered the truth. Its variance therefore adds (slope b_hi)^2, where
slope is the change per unit of b of the function's term at the reach of
the prior on `Tc` in force, half the width of its central 95% interval
(1.96 times the standard deviation of a normal, 3 by default, or the
97.5% quantile of a Student-t times its scale; the default reach for any
other prior), and b_hi the 97.5% quantile of the prior on b (about 3.5
by default): with the defaults, about 21 for `exp_threshold`, 22 for
`zeroexp_threshold` and 6.6 for `zeropower_threshold`. `const_threshold`
has no b and keeps 2. Where the data do place `Tc`, the wider prior
changes nothing; where they do not, the intervals on `Tc` and log a are
as wide as the data leave them, and the fit may take longer.

## Groups whose failures lie far apart

The intercept prior on log a is centred at the failure temperature
pooled over all individuals. Where log a differs between groups
(species, say, or acclimation temperatures), a group's coefficient on
log a is about -b times its distance in sCTmax from the reference group,
so that with a steep b a tight prior holds it back. The result is a
biased b and baseline shape while sCTmax looks right: in simulations of
a three-species design with a Weibull shape of 4, the default priors of
version 0.2.1 biased the shape of the reference species by about -0.7.

Since version 0.2.2 the defaults on log a are therefore scaled to the
groups, which removed that bias in the same simulations. The rule needs
a b to extrapolate with, so `const` and `const_threshold` keep the
priors of the section above for every term. A term of the formula of log
a marks groups when each of its variables is a factor, a character or
logical column, or a numeric column with two values, alone or in an
interaction. Where there is such a term, with b_hi the 97.5% quantile of
the prior on b (exp(-0.7 + 1.96), about 3.5, by default; read from
`intercept_priors` when the prior on `logb` given there is a normal):

- the intercept keeps its centre and takes standard deviation sqrt(2^2 +
  (d_g b_hi)^2), where d_g is the largest distance from `Tr` of any
  group's median failure temperature, a group being a level of one
  variable or a combination of levels of an interaction (for the
  exponential function; for the other functions anchored at `Tr`, the
  change per unit of b of the function's term at that distance, as in
  the section above). For a threshold function, anchored at `Tc`, the
  distance from `Tr` plays no part and the intercept takes the standard
  deviation of the section above;

- each group term's coefficients take `normal(0, s)`, with s = sqrt(1 +
  (b_hi max(D, 1))^2), where D is the range in degrees of the median
  failure temperatures over the term's groups (its levels, or the
  combinations of levels of an interaction). The floor of 1 degree is
  there because a median mixes design with effect: a treatment given
  only in the ramps can lower tolerance while its median failure
  temperature barely moves. s is a prior on the difference between
  groups, which is what the coefficient of a factor or of a 0/1 column
  is; the coefficient of a numeric column with two other values is per
  unit of the column (per percentage point of oxygen coded 30 and 100),
  so s is divided by the gap between the two values, and across an
  interaction by the product of such gaps. The model and its
  coefficients stay in the units of the data. A group none of whose
  individuals failed has no median failure temperature and is left out
  of the rule, so its coefficient takes the width the other groups set
  among themselves; where such a group may be much more tolerant than
  the others, widen its prior through `slope_prior`;

- a continuous covariate keeps the slope prior.

A group term takes the wider of this prior and the `slope_prior` that
applies to log a, so that setting `slope_prior` for the other terms
never narrows it; a `slope_prior` whose scale cannot be read (other than
a normal, a Student-t or the flat prior `""`) gives way to it. A prior
named for `loga` in `intercept_priors` replaces the rule for the
intercept. A group term given a narrower prior than the rule's must be
set through a full `prior` passed to
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).
To widen further:

    build_hazard_prior(..., intercept_priors = list(loga = "normal(-0.8, 20)"),
                       slope_prior = list(loga = "normal(0, 20)"))

Whether the priors of a fitted model pull its estimates, with these
defaults or any others, can be checked without refitting with
[`check_prior()`](https://eogmartins.github.io/rTFT/reference/check_prior.md).

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
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
ind    <- data.frame(id = 1:3, path_id = "A", event = c(1, 1, 0),
                     t_upper = c(30, 40, 60), acclim = c(10, 12.5, 15))
long   <- build_hazard_data(ind, path, "interval", covariates = "acclim",
                            checks = checks)
#> ℹ The intervals at risk were divided further, into 386 rows from 130, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.

# Default weakly-informative priors for an exponential-threshold + Weibull model.
build_hazard_prior("exp_threshold", baseline = "weibull", long = long)
#>              prior class      coef group resp dpar    nlpar   lb   ub tag
#>  normal(-1, 20.82)     b Intercept                     loga <NA> <NA>    
#>    normal(-0.7, 1)     b Intercept                     logb <NA> <NA>    
#>   normal(28.62, 3)     b Intercept                       Tc <NA> <NA>    
#>       normal(0, 1)     b Intercept                 logshape <NA> <NA>    
#>  source
#>    user
#>    user
#>    user
#>    user

# Override a couple of priors (a cooler species, a tighter amplitude) and narrow
# the slope prior of log a; every other parameter keeps its default.
build_hazard_prior("exp_threshold", baseline = "weibull",
  dpar_formulas = list(loga = ~ 1 + acclim),
  intercept_priors = list(Tc = "normal(26, 4)", loga = "normal(-6, 1.5)"),
  slope_prior = list(loga = "normal(0, 0.5)"), long = long)
#>            prior class      coef group resp dpar    nlpar   lb   ub tag source
#>  normal(-6, 1.5)     b Intercept                     loga <NA> <NA>       user
#>   normal(0, 0.5)     b                               loga <NA> <NA>       user
#>  normal(-0.7, 1)     b Intercept                     logb <NA> <NA>       user
#>    normal(26, 4)     b Intercept                       Tc <NA> <NA>       user
#>     normal(0, 1)     b Intercept                 logshape <NA> <NA>       user

# Two acclimation groups: log a has a group term, and its priors are scaled to
# how far apart the groups fail.
ind$acclim <- c(10, 15, 15)
long <- build_hazard_data(ind, path, "interval", covariates = "acclim",
                          checks = checks)
#> ℹ The intervals at risk were divided further, into 386 rows from 130, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
build_hazard_prior("exp", dpar_formulas = list(loga = ~ 1 + acclim), long = long)
#>                  prior class      coef group resp dpar nlpar   lb   ub tag
#>  normal(-4.555, 4.839)     b Intercept                  loga <NA> <NA>    
#>           normal(0, 1)     b                            loga <NA> <NA>    
#>        normal(-0.7, 1)     b Intercept                  logb <NA> <NA>    
#>       normal(0, 1.774)     b    acclim                  loga <NA> <NA>    
#>  source
#>    user
#>    user
#>    user
#>    user
```
