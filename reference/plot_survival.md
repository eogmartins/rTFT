# Plot survival along a record of temperature

Draws the survival predicted by
[`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
or
[`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md):
the probability of surviving to each time of the record, as the
posterior median with its credible interval. The result is a ggplot
object, to be refined by the user with further ggplot2 layers, scales,
labels and themes.

## Usage

``` r
plot_survival(
  x,
  temperature = FALSE,
  facet = TRUE,
  colour_by = NULL,
  interval = is.null(colour_by),
  time_unit = NULL
)
```

## Arguments

- x:

  A list, the output of
  [`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
  or
  [`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md).

- temperature:

  A single logical value, whether to draw the record of temperature, and
  of any co-stressor that changes along it, in panels above the survival
  curve. Defaults to `FALSE`.

- facet:

  A single logical value or the character string `"all"`. `TRUE` (the
  default) gives a panel per setting of `newdata`, cut down to quartile
  values when there are many; `"all"` a panel for every setting; `FALSE`
  no panels, with a colour per setting.

- colour_by:

  An optional character string naming a numeric covariate of `newdata`,
  to colour the curves by.

- interval:

  A single logical value, whether to draw the credible interval.
  Defaults to `TRUE`, and to `FALSE` with `colour_by`.

- time_unit:

  An optional character string giving the time unit for the axis label,
  e.g. `"min"`. The fit does not record the unit of its times, so none
  is shown unless given.

## Value

A ggplot object.

## Details

By default each setting of `newdata` has its own panel, following the
same rule as in
[`plot_tdt()`](https://eogmartins.github.io/rTFT/reference/plot_tdt.md).
There are two ways to draw the settings together in one panel. With
`facet = FALSE` each setting has its own colour, from teal to deep
purple in the order of `newdata`, and an entry in the legend, which
suits a few settings or a covariate that is not numeric; with more than
eight settings the colours are too close to tell apart, and a message
says to use panels or `colour_by`. With `colour_by`, the curves of every
value of one numeric covariate are coloured on a continuous scale from
teal for its lowest value to deep purple for its highest, and the
intervals are left out by default, since they overlap.

With `temperature = TRUE` the record of temperature is drawn above the
survival curve on the same time axis, so that a drop in survival can be
read against the temperatures that caused it. A co-stressor the record
follows through time is drawn in a panel of its own between the two,
named by its column, when it changes along the record; one that stays
constant is left out.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- predict_survival(fit, record, newdata = data.frame(mass = c(8, 12)))
plot_survival(pred, time_unit = "min")
plot_survival(pred, temperature = TRUE, time_unit = "min")

# both masses in one panel, the record drawn once
plot_survival(pred, temperature = TRUE, facet = FALSE, time_unit = "min")
} # }
```
