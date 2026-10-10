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

plot_check(check_survival(fit), time_unit = "min")
} # }
```
