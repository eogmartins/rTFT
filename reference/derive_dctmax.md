# Dynamic critical thermal maximum under a ramp, from a fitted model

Derives, from a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md),
the dynamic critical thermal maximum dCTmax of Jorgensen et al. (2021):
the temperature at which a fraction `p` of individuals has failed while
warmed from `T_start` at a constant `rate`. It is the critical thermal
maximum a ramping trial would record under that ramp, with its posterior
median and credible interval, at the covariate values in each row of
`newdata`. Random effects are set to zero, so the results are those of a
typical group.

## Usage

``` r
derive_dctmax(
  fit,
  rate,
  T_start,
  newdata = NULL,
  p = 0.5,
  levels = NULL,
  prob = 0.95,
  n_draws = 1000,
  T_max = NULL
)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- rate:

  A numeric vector of ramp rates, in degrees per unit of time of the
  data the model was fitted to (e.g. degrees C per minute), all
  positive.

- T_start:

  A numeric vector of starting temperatures, one per rate or a single
  one for all.

- newdata:

  A data frame with one row per set of covariate values to report at,
  holding every covariate in the parameters' formulas. `NULL` (the
  default) suits a model with no covariates.

- p:

  A numeric vector of fractions of individuals failing, each strictly
  between 0 and 1. Defaults to 0.5.

- levels:

  A named numeric vector giving the level of each co-stressor, on its
  measurement scale, e.g. `c(DO = 6)`. Required when the model has
  co-stressors.

- prob:

  A single number between 0 and 1, the coverage of the credible
  intervals. Defaults to 0.95.

- n_draws:

  A single number, the number of posterior draws to use, evenly spaced
  through the posterior. Defaults to 1000.

- T_max:

  A single number, the temperature up to which the ramp is followed.
  `NULL` (the default) takes 15 degrees above the warmest temperature
  the model was fitted to.

## Value

A tibble with a row for each row of `newdata`, ramp and `p`: the columns
of `newdata`, `rate`, `T_start`, `p`, and the posterior `median`,
`lower` and `upper` of dCTmax, with `prob`, `n_draws` and `n_unreached`.

## Details

Unlike the static critical thermal maximum from
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md),
dCTmax depends on the ramp. Individuals accrue hazard throughout the
ramp, so a slower one gives them more time to fail before a given
temperature is reached, and its dCTmax is lower. Giving several rates
shows how much. The ramp need not be one that was run: the model was
fitted to failure times, from a static or a ramping assay, and is
applied here to any ramp.

Under the ramp T(t) = T_start + rate t, the cumulative hazard is the
integral of M(T(t)) h0(t) from the start of the ramp, and dCTmax solves,
for each posterior draw, cumulative hazard = -log(1 - p). With a
baseline other than the memoryless one, the hazard also depends on the
time elapsed since the ramp began.

The ramp is followed up to `T_max`. A posterior draw whose fraction
failed has not reached `p` by then counts as infinite, as in
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md),
and is counted in `n_unreached`. A message is given when a dCTmax lies
above the temperatures the model was fitted to, where the data do not
constrain the failure-rate function.

## References

Jorgensen, L. B., Malte, H., Orsted, M., Klahn, N. A. and Overgaard, J.
(2021). A unifying model to estimate thermal tolerance limits in
ectotherms across static, dynamic and fluctuating exposures to thermal
stress. *Scientific Reports* 11, 12840.

## See also

[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
for the static critical thermal maximum, and
[`average_dctmax()`](https://eogmartins.github.io/rTFT/reference/average_dctmax.md)
to average dCTmax over several models.

## Examples

``` r
if (FALSE) { # \dontrun{
# the CTmax expected under three ramps from 14 degrees
derive_dctmax(fit, rate = c(0.1, 0.2, 0.3), T_start = 14,
              newdata = data.frame(mass = 15))
} # }
```
