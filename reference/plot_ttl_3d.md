# Plot the thermal tolerance landscape as a surface to rotate

Draws the thermal tolerance landscape from
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
in three dimensions: elapsed exposure time and constant temperature
across the floor, and the posterior median survival as the height. The
figure is interactive, so it can be rotated, zoomed and read by
hovering, in a viewer, a browser or an HTML document.

## Usage

``` r
plot_ttl_3d(x, row = NULL, interval = FALSE, grid = FALSE, time_unit = NULL)
```

## Arguments

- x:

  A list, the output of
  [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md).

- row:

  An optional single whole number, the setting of `newdata` to draw, as
  its row number among the settings. `NULL` (the default) draws the
  first, with a message when there are several.

- interval:

  A single logical value, whether to add the lower and upper bounds of
  the credible interval as two faint surfaces around the median.

- grid:

  A single logical value, or one or two whole numbers. `TRUE` draws a
  grid on the median surface, ten lines at constant temperature and ten
  at constant elapsed time. A number sets how many lines each way, and
  two numbers the lines at constant temperature and at constant time,
  e.g. `c(5, 15)`. `FALSE` (the default) draws none.

- time_unit:

  An optional character string giving the time unit for the axis label,
  e.g. `"min"`.

## Value

A plotly object.

## Details

It needs the plotly package, and returns a plotly object rather than a
ggplot, so it is refined with plotly's own functions (for example
`plotly::layout()`), not with ggplot2 layers. A three-dimensional
surface hides part of itself in a static figure; for print,
[`plot_ttl()`](https://eogmartins.github.io/rTFT/reference/plot_ttl.md)
shows the same landscape as a map.

One setting of `newdata` is drawn at a time. Elapsed time is on a log
scale, unless the times include zero.

## Examples

``` r
if (FALSE) { # \dontrun{
# a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
# degrees over 60 minutes and checked every 10, in two trials (see fit_hazard())
path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))
ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
                     trial = rep(c("a", "b"), each = 10),
                     t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
                            checks = checks)
fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))

tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
plot_ttl_3d(tol, row = 2, time_unit = "min")
} # }
```
