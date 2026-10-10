# Plot the thermal death time curves

Draws the thermal death time curves from
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md):
the time until a fraction `p` of individuals has failed, against
constant temperature, as the posterior median with its credible
interval. The result is a ggplot object, which be refined by the user
with further ggplot2 layers, scales, labels and themes.

## Usage

``` r
plot_tdt(
  x,
  facet = TRUE,
  colour_by = NULL,
  interval = is.null(colour_by),
  time_unit = NULL
)
```

## Arguments

- x:

  A list, the output of
  [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md).

- facet:

  A single logical value or the character string `"all"`. `TRUE` (the
  default) gives a panel per setting of `newdata`, cut down to quartile
  values when there are many; `"all"` a panel for every setting; `FALSE`
  no panels.

- colour_by:

  An optional character string naming a numeric covariate of `newdata`,
  to colour the curves by instead of by `p`.

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

By default each fraction `p` has its own colour, from blue for the
fewest failed to red for the most, and each setting of `newdata` its own
panel. Each setting has its own curve, so settings are shown or left
out, never averaged: with more than six, a numeric covariate with more
than five values is cut down to the values nearest its minimum,
quartiles and maximum, and a message says which were kept.

With `colour_by`, the curves of every value of one numeric covariate are
drawn together, coloured on a continuous scale from teal for its lowest
value to deep purple for its highest, with a panel per fraction `p`.
This shows how the curve shifts with the covariate. The intervals are
then left out by default, since they overlap. Where the covariate is
itself a temperature, such as an acclimation temperature, a cool-to-warm
scale may read better, and can be added to the result:
`+ ggplot2::scale_colour_gradient(low = "#2C7BB6", high = "#D7191C")`.

Time is on a log scale, on which an exponential failure rate gives
straight lines. Where a share of the posterior never reaches the
fraction, the interval runs to the top of the panel, and where the
median itself is never reached the line stops.

## Examples

``` r
if (FALSE) { # \dontrun{
tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
plot_tdt(tol, time_unit = "min")

# how the curves shift with body mass
tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = seq(5, 15, by = 0.5)),
                      t_ref = 60)
plot_tdt(tol, colour_by = "mass", time_unit = "min")
} # }
```
