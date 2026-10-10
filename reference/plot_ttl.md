# Plot the thermal tolerance landscape

Draws the thermal tolerance landscape from
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md):
survival over constant temperature and elapsed exposure time. The result
is a ggplot object, to be refined by the user with further ggplot2
layers, scales, labels and themes.

## Usage

``` r
plot_ttl(
  x,
  type = c("surface", "curves"),
  facet = TRUE,
  Ta = NULL,
  survival = c(0.1, 0.5, 0.9),
  labels = TRUE,
  time_unit = NULL
)
```

## Arguments

- x:

  A list, the output of
  [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md).

- type:

  A character string, either `"surface"` (the default) for the landscape
  as a map, or `"curves"` for survival curves at a few temperatures.

- facet:

  A single logical value or the character string `"all"`. `TRUE` (the
  default) gives a panel per setting of `newdata`, cut down to quartile
  values when there are many; `"all"` a panel for every setting; `FALSE`
  no panels.

- Ta:

  An optional numeric vector of the temperatures to draw for
  `type = "curves"`, each replaced by the nearest temperature of the
  landscape's grid. `NULL` (the default) takes the grid temperatures
  nearest its minimum, quartiles and maximum.

- survival:

  A numeric vector of survival levels between 0 and 1, at which to draw
  the contour lines for `type = "surface"`.

- labels:

  A single logical value, whether to label each contour line of
  `type = "surface"` with its survival level.

- time_unit:

  An optional character string giving the time unit for the axis label,
  e.g. `"min"`. The fit does not record the unit of its times, so none
  is shown unless given.

## Value

A ggplot object.

## Details

With `type = "surface"` the whole landscape is drawn as a map, elapsed
time across and temperature up, coloured by the posterior median
survival from red (most have failed) through pale yellow (a survival of
one half) to blue (most survive), with contour lines at the survival
levels in `survival`, each labelled with its level. The credible
interval is not shown on it. With `type = "curves"` survival is drawn
against elapsed time at a few temperatures, each as the posterior median
with its credible interval.

Elapsed time is on a log scale, unless the times include zero. The
panels follow the same rule as in
[`plot_tdt()`](https://eogmartins.github.io/rTFT/reference/plot_tdt.md),
except that the map needs a panel for each setting of `newdata`, so
`facet = FALSE` is refused for it when there are several.

## Examples

``` r
if (FALSE) { # \dontrun{
tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
plot_ttl(tol, time_unit = "min")
plot_ttl(tol, type = "curves", Ta = c(22, 25, 28), time_unit = "min")
} # }
```
