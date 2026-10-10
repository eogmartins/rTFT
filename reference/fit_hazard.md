# Fit a failure-rate hazard model with brms

Assembles the formula and family with
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md),
the default prior with
[build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)
(unless you supply your own), selects the Stan code the custom families
need, and hands it all to
[brms::brm](https://paulbuerkner.com/brms/reference/brm.html). Anything
in `...` is forwarded to `brm` (chains, iterations, cores, backend,
initial values, and so on).

## Usage

``` r
fit_hazard(
  long,
  m_form,
  baseline = "mem",
  dpar_formulas = list(),
  observation = c("interval", "continuous", "mixed"),
  prior = NULL,
  intercept_priors = list(),
  slope_prior = NULL,
  sd_prior = NULL,
  costressors = list(),
  k_sw = NULL,
  ...
)
```

## Arguments

- long:

  A data frame, the row-per-interval table from
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md).

- m_form, baseline, dpar_formulas, observation:

  As in
  [build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md).

- prior:

  An optional `brmsprior` for full manual control. When `NULL` (the
  default) a weakly-informative prior is built with
  [build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md).
  A supplied prior is used as given, in place of that default rather
  than merged with it, so any coefficient it leaves out falls back to
  brms's own default, not the weakly-informative one. To adjust only
  some priors and keep the rest, use `intercept_priors` / `slope_prior`
  / `sd_prior` instead of this.

- intercept_priors, slope_prior, sd_prior:

  As in
  [build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md),
  forwarded when the default prior is built (ignored if you supply
  `prior`). Each changes only the priors it names and leaves the rest at
  their defaults, so a user who only calls `fit_hazard` can retune
  without building a prior by hand.

- costressors:

  An optional list of
  [costressor](https://eogmartins.github.io/rTFT/reference/costressor.md)
  objects, as in
  [build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md),
  forwarded to both the formula and the default prior. The same list
  must have been passed to
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
  so `long` carries the `z_<name>` columns the formula reads. Defaults
  to none.

- k_sw:

  An optional softening steepness for a threshold function's gate, as in
  [build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md).
  An advanced override: leave `NULL` (the default) to use the value in
  [M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md).
  A fit made with a custom value must be reported with the same one.

- ...:

  Further arguments forwarded to
  [brms::brm](https://paulbuerkner.com/brms/reference/brm.html). Any
  argument it accepts (e.g. `chains`, `iter`, `cores`, `backend`,
  `init`), except `formula`, `data`, `family`, `prior` and `stanvars`,
  which `fit_hazard` sets itself. See
  [brms::brm](https://paulbuerkner.com/brms/reference/brm.html) for the
  full set.

## Value

A `brmsfit`, carrying in `$rtft` the settings it was fitted with: the
failure-rate and baseline functions, the observation mode, the reference
temperature `Tr`, the softening steepness, the co-stressors, the ranges
of temperature and elapsed exposure time the data cover, and the
individual each row of `long` belongs to. The reporting functions read
them from there, so they are never re-entered by hand.

## Details

A missing value in any variable the model uses is an error. brms would
drop the rows concerned without saying so, leaving some individuals with
only part of their exposure, so remove those individuals or fill in the
values first.

Before fitting, the temperature is checked against the reference `Tr`
for the failure rates that take a logarithm or fractional power of
`Ta - Tr` (`zeroexp`, `power`, `zeropower`): these are undefined where
the temperature reaches `Tr` (`zeroexp`, `zeropower`) or drops below it
(`power`), so a run that would hit that is stopped with a clear error
rather than a cryptic sampler failure.

## Examples

``` r
if (FALSE) { # \dontrun{
path <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
ind  <- data.frame(id = 1:20, path_id = "A", event = 1,
                   t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
long <- build_hazard_data(ind, path, "interval")
fit  <- fit_hazard(long, "exp", chains = 4, iter = 2000, cores = 4)
} # }
```
