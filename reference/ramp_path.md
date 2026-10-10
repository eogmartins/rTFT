# Build a linear ramp temperature trajectory

A convenience function for the `path` argument of
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
useful when there is no logged temperature series: when temperature
follows a controlled linear ramp, `temp(t) = T_start + rate * t`, the
nominal ramp stands in for a logged trace. It is vectorised over its
arguments: give several ramps at once and it stacks them into one table
keyed by `path_id`.

## Usage

``` r
ramp_path(path_id, T_start, rate, t_max, by = NULL)
```

## Arguments

- path_id, T_start, rate, t_max:

  The ramp specification, one value per ramp: `path_id` the identifier,
  `T_start` the start temperature, `rate` in degrees per time unit, and
  `t_max` the duration. The four are recycled together, so equal-length
  vectors describe several ramps and a length-one argument is reused
  across all of them.

- by:

  An optional numeric giving the sampling step in time. When `NULL` (the
  default) the path is just its two endpoints, which the linear
  interpolation in
  [build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  reproduces exactly. A value such as `by = 5` creates a row every `by`
  time units, which then serves as the check schedule in interval mode
  or a ready-made integration grid in continuous mode.

## Value

A data frame with `path_id`, `time` and `temp`, ready to pass to
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
as its `path`.

## Examples

``` r
# One ramp: 20 -> 35 degrees over 60 minutes, as a two-point path.
ramp_path("A", T_start = 20, rate = 0.25, t_max = 60)
#> # A tibble: 2 × 3
#>   path_id  time  temp
#>   <chr>   <dbl> <dbl>
#> 1 A           0    20
#> 2 A          60    35

# Two ramps at different rates, each sampled every 10 minutes (T_start recycled).
ramp_path(c("slow", "fast"), T_start = 20, rate = c(0.2, 0.4),
          t_max = c(75, 40), by = 10)
#> # A tibble: 14 × 3
#>    path_id  time  temp
#>    <chr>   <dbl> <dbl>
#>  1 slow        0    20
#>  2 slow       10    22
#>  3 slow       20    24
#>  4 slow       30    26
#>  5 slow       40    28
#>  6 slow       50    30
#>  7 slow       60    32
#>  8 slow       70    34
#>  9 slow       75    35
#> 10 fast        0    20
#> 11 fast       10    24
#> 12 fast       20    28
#> 13 fast       30    32
#> 14 fast       40    36
```
