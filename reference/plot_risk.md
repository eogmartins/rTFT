# Plot each period's share of the risk along a record

Draws the `risk` table of
[`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
or
[`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md):
the share of the hazard accrued over the record that falls in each
period, as a bar at the posterior median with its credible interval. It
shows where along the record the risk is concentrated. The result is a
ggplot object, which can be refined by the user with further ggplot2
layers, scales, labels and themes.

## Usage

``` r
plot_risk(x, facet = TRUE, interval = TRUE)
```

## Arguments

- x:

  A list, the output of
  [`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
  or
  [`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md)
  called with `by`.

- facet:

  A single logical value or the character string `"all"`. `TRUE` (the
  default) gives a panel per setting of `newdata`, cut down to quartile
  values when there are many; `"all"` a panel for every setting; `FALSE`
  no panels, with the settings side by side in each period.

- interval:

  A single logical value, whether to draw the credible interval.
  Defaults to `TRUE`.

## Value

A ggplot object.

## Details

The shares are of whatever risk there is, so they say where it falls and
not how large it is; read them with the `mortality` table. A period has
no bar where no posterior draw has any hazard along the record.

By default each setting of `newdata` has its own panel, following the
same rule as in
[`plot_tdt()`](https://eogmartins.github.io/rTFT/reference/plot_tdt.md).
With `facet = FALSE` the settings are drawn side by side within each
period, each in its own colour, from teal to deep purple in the order of
`newdata`; this needs eight settings or fewer.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- predict_survival(fit, record, newdata = data.frame(mass = c(8, 12)),
                         by = "day")
plot_risk(pred)
plot_risk(pred, facet = FALSE)
} # }
```
