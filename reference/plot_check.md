# Plot the check of a fitted model against observed survival

Draws the output of
[`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md):
for each group, the proportion of individuals surviving in datasets
simulated from the fitted model, as its median and band, with the
proportion observed to survive as a stepped line over it. Where the
model describes the data, the stepped line stays within the band. The
result is a ggplot object, which can be refined by the user with further
ggplot2 layers, scales, labels and themes.

## Usage

``` r
plot_check(x, time_unit = NULL)
```

## Arguments

- x:

  A tibble, the output of
  [`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md).

- time_unit:

  An optional character string giving the time unit for the axis label,
  e.g. `"min"`. The fit does not record the unit of its times, so none
  is shown unless given.

## Value

A ggplot object, with a panel per group.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_check(check_survival(fit), time_unit = "min")
} # }
```
