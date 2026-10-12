# Declare a co-stressor for a failure-rate hazard model

Describes one non-temperature environmental variable (Z) the individual
is exposed to and that co-determines the failure rate: dissolved oxygen,
salinity, pH, a toxicant concentration, or, for a terrestrial ectotherm,
humidity. The failure rate then becomes a surface M(T, Z). The
co-stressor enters log M(T) as its own term and, when asked for with
`interaction = TRUE`, an interaction with temperature, so a crossed
temperature-by-stressor design (a ramp run at two or more stressor
levels) can estimate how the stressor changes the thermal response.
Build one for each co-stressor and pass them together, as a list, to
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md),
[build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)
and
[fit_hazard](https://eogmartins.github.io/rTFT/reference/fit_hazard.md),
so all four read the same declaration.

## Usage

``` r
costressor(
  name,
  column,
  ref,
  scale,
  direction,
  interaction = FALSE,
  trajectory = NULL,
  transform = NULL
)
```

## Arguments

- name:

  A character string, a short name for the co-stressor, used to name its
  column (`z_<name>`) and its parameters (`c<name>`, `dT<name>`).
  Letters and digits only, starting with a letter, because brms allows
  no dot or underscore in a parameter name. For example `"DO"` gives the
  column `z_DO` and the parameters `cDO` and `dTDO`.

- column:

  A character string naming the column that holds the raw measurement
  `X`: a column of `ind` (one value per individual) for a stressor held
  fixed within a trial, or a column of `trajectory` for one that changes
  over the trial.

- ref:

  The benign, no-additional-stress reference level of the raw
  measurement (air-saturation for dissolved oxygen, zero for a toxicant,
  the control level for salinity or pH). The failure rate reduces to the
  pure M(T) here.

- scale:

  A single positive number, the unit a coefficient is expressed per (per
  mg/L or per kPa for oxygen, per unit concentration for a toxicant).

- direction:

  Either `1` or `-1`, set so that a positive `z` always means more
  stress: `-1` for dissolved oxygen (less oxygen is more stress), `1`
  for a toxicant (more is more stress).

- interaction:

  A single logical value. When `FALSE` (the default) the co-stressor
  enters through its main effect alone, changing the failure rate by the
  same proportion at every temperature. When `TRUE` it also enters
  through a temperature interaction `dT<name>`, so it can change the
  temperature slope of log M. For a multiplicative alternative that
  rescales thermal sensitivity instead, see ‘Two kinds of temperature
  interaction’ below, and for what an interaction implies outside the
  temperatures of the data, see ‘An interaction outside the temperatures
  of the data’.

- trajectory:

  An optional data frame giving the co-stressor over time, with columns
  `path_id`, `time` and the `column` named above, in the shape of the
  `path` argument of
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md).
  Supply it for a stressor that changes during a trial; it is
  interpolated to the same interval midpoints as temperature. Leave
  `NULL` for a stressor held fixed within a trial, which is then read
  from `ind`.

- transform:

  An optional function of the standardised deviation (X - ref) / scale
  that returns the column `z`. Use it for a stressor whose harm is not
  signed-linear, for example a squared deviation for a U-shaped response
  to salinity or pH, where a departure either way is harmful. When
  `NULL` (the default) the signed-linear deficit direction \* (X - ref)
  / scale is used.

## Value

A `costressor` object: a named list holding the declaration, for the
`costressors` argument of
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md),
[build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)
and
[fit_hazard](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

## Details

The raw measurement `X` becomes a standardised column z = direction \*
(X - ref) / scale, which is zero at the benign reference. There the
failure rate reduces to the pure thermal M(T): the amplitude `loga` and
critical temperature `Tc` keep their unstressed meaning, and a model
with the co-stressor held at its reference is the same as one without
it. This follows from centring the column at a chosen reference rather
than at the sample mean.

The three pieces the standardisation needs (a reference, a scale, and a
direction) are all required, as none has a safe default. A domain-chosen
scale (per mg/L of oxygen, per unit concentration) is preferred to a
data standard deviation, which is unstable with only two or three
treatment levels. As an example, for dissolved oxygen with a reference
at air-saturation (ref = 8.5 mg/L), a scale of 1 mg/L (scale = 1), and
less oxygen meaning more stress (direction = -1), z = -1 \* (X - 8.5) /
1 = 8.5 - X, the oxygen deficit is: z = 0 at saturation and z = 2.5 at 6
mg/L (below saturation). The direction re-orients every stressor onto a
common axis on which a larger z always means more stress, so the main
effect `c<name>` has the same sign meaning for every stressor: a
positive value is a hazard that rises under stress.

This signed-linear column assumes the stressor harm runs one way. For a
stressor that is harmful in either direction (salinity or pH, where a
departure above or below the optimum both add stress), supply a
`transform` function so z measures distance from the reference rather
than a signed deficit: a squared standardised deviation, transform =
function(u) u^2, makes z = ((X - ref) / scale)^2, zero at the optimum
and rising on both sides. The `direction` then no longer applies, since
z is symmetric.

The co-stressor contributes one or two parameters, named from `name`.
The main effect `c<name>` is the change in the log failure rate per
`scale` of stressor beyond the reference, and the interaction `dT<name>`
(with `interaction = TRUE`) is how that change varies with temperature.
These are the names brms estimates and that
[build_hazard_prior](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)
sets priors on. `c<name>` is a continuous dose term: with three or more
stressor levels it is a per-`scale` slope through them, but with only
two it cannot be told apart from a per-treatment shift of the amplitude
`loga`. A main effect on its own is in fact just that shift: it is
equivalent to a covariate on the amplitude, `loga ~ 1 + z` in
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md)'s
`dpar_formulas`. The additive `dT<name>` term is not: no `dpar_formulas`
entry reproduces it, so it is the part `costressor()` supplies. A
temperature interaction of a different shape can still be had through
`dpar_formulas` (the multiplicative form in ‘Two kinds of temperature
interaction’ below); what is unique to `costressor()` is this particular
additive term.

Because `c<name>` and `dT<name>` are ordinary parameters,
`dpar_formulas` in
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md)
can put covariates on either, so the stressor's effect can itself vary:
`cDO = ~ 1 + acclim` lets the oxygen sensitivity depend on acclimation
(a stressor-by-covariate interaction), and the same on `dT<name>` lets
the temperature interaction vary. These are higher-order terms, so they
need a design that varies temperature, the stressor and the covariate
together enough to identify them; the formula is built either way, but a
thin design shows up as poor mixing rather than a clean fit.

## Categorical stressor treatments

`costressor()` expects a stressor measured on a numeric scale, so `ref`,
`scale` and `direction` mean something. Quantitative treatment levels
belong here even though they are "treatments": dissolved oxygen dosed at
8.5, 5 and 2 mg/L is numeric, so pass the measured level as `column`. A
genuinely unordered categorical stressor, with no scale on which a
reference and a per-unit step make sense, does not. Model it through
`dpar_formulas` in
[build_hazard_spec](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md)
instead: a factor on the amplitude `loga` gives a per-treatment level
shift, the same as a stressor main effect. To let the treatment effect
vary with temperature, put the factor on the parameter that governs the
temperature response (`logb` for the exponential forms, `Tc` for the
threshold forms). That is a genuine temperature-by-treatment
interaction, though a different parameterisation from the co-stressor's
`dT<name>` term, as it scales the temperature slope rather than adding
to it (see ‘Two kinds of temperature interaction’ below). A two-level
treatment can instead be encoded 0/1 and passed to `costressor()`, which
then gives the co-stressor's own main and interaction terms.

## An interaction outside the temperatures of the data

The interaction is off by default because of how it behaves away from
the temperatures at which individuals failed. It makes the effect of the
stressor a straight line in temperature, c + dT (Ta - Tr), which is
informed only over those temperatures: a few degrees, in a ramping
assay. Extended below them the line reaches zero at Tr - c / dT and
changes sign beyond, so a model in which the stressor is harmful where
individuals failed may have it protective at cooler temperatures, such
as those of a field record. Use the interaction to describe how the
stressor's effect changes within the temperatures of the assay, and
compare what the models with and without it imply at the temperatures a
prediction is needed for (see the vignette *Temperature and a
co-stressor*).

## Two kinds of temperature interaction

The `dT<name>` interaction is additive on the temperature slope: it adds
`dT<name> * (Ta - Tr) * z` to log M, so the per-degree effect shifts by
a constant amount per unit of stressor. This is the standard
temperature-by-stressor interaction, it works the same on every M(T)
function, and it lets the slope change sign, so a stressor can in
principle flatten or reverse the thermal effect. The sign of `dT<name>`
then reads as the synergy or antagonism of temperature and the stressor
(on the proportional-hazards scale): positive steepens the thermal
slope, negative flattens it, and a large enough negative value reverses
it (an antagonistic crossover where warming becomes protective under
severe harm caused by the co-stressor). The alternative is to let the
stressor rescale thermal sensitivity by a proportion, which keeps it
positive: set `interaction = FALSE` (so the co-stressor adds only its
main effect and the `z_<name>` column) and put that column on the
sensitivity parameter through `dpar_formulas` instead:
`logb ~ 1 + z_<name>` for the exponential forms (for the threshold forms
a factor on `Tc` shifts the threshold location, a different effect).
Choose the additive form for a conventional interaction and the freedom
for the slope to change sign (e.g., an antagonistic crossover), the
multiplicative form when the biology is a proportional change in
sensitivity that should stay positive (it can attenuate the thermal
slope but not reverse it).

## Examples

``` r
# Dissolved oxygen held at a fixed level per trial: benign at air-saturation
# (about 8.5 mg/L), expressed per mg/L, with less oxygen meaning more stress. It
# enters through its main effect alone.
costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1)
#> $name
#> [1] "DO"
#> 
#> $column
#> [1] "do"
#> 
#> $ref
#> [1] 8.5
#> 
#> $scale
#> [1] 1
#> 
#> $direction
#> [1] -1
#> 
#> $interaction
#> [1] FALSE
#> 
#> $trajectory
#> NULL
#> 
#> $transform
#> NULL
#> 
#> attr(,"class")
#> [1] "rTFT_costressor"

# The same with its interaction with temperature.
costressor("DO", column = "do", ref = 8.5, scale = 1, direction = -1,
           interaction = TRUE)
#> $name
#> [1] "DO"
#> 
#> $column
#> [1] "do"
#> 
#> $ref
#> [1] 8.5
#> 
#> $scale
#> [1] 1
#> 
#> $direction
#> [1] -1
#> 
#> $interaction
#> [1] TRUE
#> 
#> $trajectory
#> NULL
#> 
#> $transform
#> NULL
#> 
#> attr(,"class")
#> [1] "rTFT_costressor"

# A toxicant that changes over the trial, supplied as its own trajectory.
tox_path <- data.frame(path_id = "A", time = c(0, 60), conc = c(0, 4))
costressor("Tox", column = "conc", ref = 0, scale = 1, direction = 1,
           trajectory = tox_path)
#> $name
#> [1] "Tox"
#> 
#> $column
#> [1] "conc"
#> 
#> $ref
#> [1] 0
#> 
#> $scale
#> [1] 1
#> 
#> $direction
#> [1] 1
#> 
#> $interaction
#> [1] FALSE
#> 
#> $trajectory
#>   path_id time conc
#> 1       A    0    0
#> 2       A   60    4
#> 
#> $transform
#> NULL
#> 
#> attr(,"class")
#> [1] "rTFT_costressor"

# Salinity, harmful above or below an optimum: a squared standardised deviation
# makes z measure distance from the optimum (about 15 PSU, per 5 PSU) in either
# direction. With a transform supplied, `direction` no longer applies.
costressor("Sal", column = "sal", ref = 15, scale = 5, direction = 1,
           transform = function(u) u^2)
#> $name
#> [1] "Sal"
#> 
#> $column
#> [1] "sal"
#> 
#> $ref
#> [1] 15
#> 
#> $scale
#> [1] 5
#> 
#> $direction
#> [1] 1
#> 
#> $interaction
#> [1] FALSE
#> 
#> $trajectory
#> NULL
#> 
#> $transform
#> function (u) 
#> u^2
#> <environment: 0x55724c12ea88>
#> 
#> attr(,"class")
#> [1] "rTFT_costressor"
```
