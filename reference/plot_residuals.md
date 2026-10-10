# Plot the residuals of a fitted model

Draws the output of
[`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md).
With `against = NULL`, the residuals are set against the values expected
of a uniform distribution, which is what they follow when the model
describes the data: the points then lie along the diagonal. With
`against` naming a column, the residuals are drawn against it, to look
for a pattern the model has missed: points for a numeric variable, such
as body mass, and boxes for a grouping, such as the trial or the
treatment. A horizontal line marks 0.5, the middle of the residuals
under the model; residuals below it are individuals that failed earlier
than the model predicts. The time of failure is not a variable to draw
them against (see
[`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)).
The result is a ggplot object, which can be refined by the user with
further ggplot2 layers, scales, labels and themes.

## Usage

``` r
plot_residuals(x, against = NULL)
```

## Arguments

- x:

  A tibble, the output of
  [`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md).

- against:

  An optional character string naming a column of `x` to draw the
  residuals against. `NULL` (the default) compares them with the uniform
  distribution.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- hazard_residuals(fit)
plot_residuals(res)
plot_residuals(res, against = "mass")
plot_residuals(res, against = "trial")
} # }
```
