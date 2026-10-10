# Thermal tolerance derived from a fitted failure-rate model

Derives, from a model fitted with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md),
the failure-rate parameters on their natural scale, the static critical
thermal maximum sCTmax, the thermal death time curves and the thermal
tolerance landscape, each with its posterior median and credible
interval, at the covariate values in each row of `newdata`. The model's
settings (failure-rate and baseline functions, reference temperature,
softening steepness, co-stressors) are read from the fit. Random effects
are set to zero, so the results are those of a typical group.

## Usage

``` r
derive_tdt_ttl(
  fit,
  newdata = NULL,
  p = c(0.1, 0.5, 0.9),
  t_ref,
  Ta = NULL,
  tau = NULL,
  levels = NULL,
  prob = 0.95,
  n_draws = 1000
)
```

## Arguments

- fit:

  A `brmsfit`, a model fitted with
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- newdata:

  A data frame of covariate values, one row per setting at which to
  derive. It must hold every covariate in the model's parameter
  formulas. May be `NULL` for a model without covariates.

- p:

  A numeric vector of the fractions of individuals failing at which
  sCTmax and the thermal death time curves are derived, each strictly
  between 0 and 1.

- t_ref:

  A single positive number, the reference exposure time for sCTmax, in
  the time unit of the data the model was fitted to (for example `60`
  for one hour when times are in minutes). It has no default, because
  the fit does not record that unit.

- Ta:

  An optional numeric vector of the constant temperatures for the
  thermal death time curves and the tolerance landscape. `NULL` (the
  default) takes 50 evenly spaced across the temperatures the data
  cover.

- tau:

  An optional numeric vector of the elapsed exposure times for the
  tolerance landscape, in the time unit of the data, none negative.
  `NULL` (the default) takes 50 spaced evenly on a log scale across the
  times the data cover. A time of zero is allowed, and gives a survival
  of 1, so that a survival curve can start there.

- levels:

  A named numeric vector giving the level of each co-stressor, on its
  measurement scale, e.g. `c(DO = 6)`. Required when the model has
  co-stressors.

- prob:

  A single number strictly between 0 and 1, the coverage of the credible
  intervals.

- n_draws:

  An optional single number, the number of posterior draws to use,
  evenly spaced through the posterior. `NULL` uses them all.

## Value

A list of four tibbles, each starting with the columns of `newdata` and
ending with the posterior summary `median`, `lower`, `upper`, `prob`,
`n_draws` and `n_unreached`:

- `params`: one row per parameter, named in `quantity`;

- `sctmax`: one row per value of `p`, with `p` and `t_ref`;

- `tdt_curves`: one row per value of `p` and temperature, with `p` and
  `Ta`; the summary is of the time until that fraction has failed;

- `ttl_surface`: one row per temperature and elapsed time, with `Ta` and
  `tau`; the summary is of survival.

## Details

sCTmax is the constant temperature at which a fraction `p` of
individuals fail within the reference exposure time `t_ref`, the static
critical thermal maximum of Jorgensen et al. (2021), who report it for
one hour. Under a constant temperature, survival is exp(-M(T)
H0(t_ref)), so for each posterior draw sCTmax solves M(T) = -log(1 - p)
/ H0(t_ref). It is not the critical thermal maximum a ramping experiment
records (their dynamic CTmax).

The parameters reported are a, and where the failure-rate function has
them b, c and Tc; the baseline's shape (Weibull, log-logistic), its
scale (log-logistic; the Weibull's is fixed at the reference time, see
[h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md))
or its rate constant gamma (Gompertz); and each co-stressor's
coefficients. With a co-stressor that interacts with temperature,
`b_at_levels` is the thermal sensitivity b + dT z at the chosen levels.
The thermal sensitivity z, the temperature change that shortens the time
to failure tenfold, is reported only where the thermal death time curve
is a straight line in log time against temperature, so that one number
describes it: the exponential failure rate with the memoryless or
Weibull baseline, where z = k log(10) / b (k = 1 for the memoryless
baseline, the Weibull shape otherwise). Elsewhere the curve bends and no
single z stands for it.

The thermal death time curve for a fraction `p` gives, at each constant
temperature in `Ta`, the time until that fraction has failed: t_p =
H0^-1(-log(1 - p) / M(T)). The thermal tolerance landscape gives
survival exp(-M(T) H0(tau)) over the grid of temperatures `Ta` and
elapsed exposure times `tau`. "Thermal death time" is the literature's
name for the curve, whatever the failure endpoint.

## Draws that never reach the quantity

A posterior draw for which a quantity is never reached (an sCTmax the
rate cannot attain within 1000 degrees of the reference temperature, or
at all; a time to failure that never comes, under a Gompertz baseline
whose cumulative hazard levels off) counts as infinite rather than being
dropped, so the summaries describe the whole posterior. `n_unreached`
counts these draws, and a message lists the quantities that have any.
With coverage `prob`, let u = (1 - prob) / 2 (0.025 for a 95% interval):

- fewer than a share u of the draws unreached: both bounds are finite;
  report the summary as it stands.

- between u and half: the outer bound is infinite; report it as
  unbounded (for example "above 40 degC", "not reached"), and consider
  whether the question lies beyond the temperatures or exposure times
  the data cover.

- half or more: the median is infinite; the quantity is not estimable
  there. Ask at a lower `p` or a shorter `t_ref`.

## References

Jorgensen, L. B., Malte, H., Orsted, M., Klahn, N. A. and Overgaard, J.
(2021). A unifying model to estimate thermal tolerance limits in
ectotherms across static, dynamic and fluctuating exposures to thermal
stress. *Scientific Reports* 11, 12840.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), p = 0.5,
                      t_ref = 60)
tol$sctmax
} # }
```
