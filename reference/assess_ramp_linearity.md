# Assess how linear a logged temperature trace is

A pre-fitting diagnostic with two uses. As quality control on the ramp,
it fits `temp ~ time` per path and reports the realised start and rate,
the fit (`r_squared`), and, given `nominal`, the `rate_bias` against the
setpoint, so you can see whether the controller held a linear ramp at
the intended rate. And as a check on whether a dense or noisy log can be
simplified for fitting: a finely sampled temperature series turns every
reading into an at-risk-interval boundary in
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
so a linear one can be replaced by a clean two-point ramp from
[ramp_path](https://eogmartins.github.io/rTFT/reference/ramp_path.md)
(exact under the linear interpolation, free of sensor jitter, and far
fewer rows). A faithful series is otherwise best fitted directly; this
only judges that faithfulness.

## Usage

``` r
assess_ramp_linearity(path, nominal = NULL)
```

## Arguments

- path:

  A data frame with `path_id`, `time` and `temp` – the logged
  temperature trace(s).

- nominal:

  An optional named numeric of setpoint rates (names matching
  `path_id`). When given, adds `rate_nominal` and `rate_bias` (realised
  minus nominal), the controller's rate error.

## Value

A tibble, one row per `path_id`, with `n` (the number of readings), the
realised `T_start` and `rate`, `r_squared`, `rmse`, `max_abs_dev`, and
`max_dev_hot` (the worst absolute deviation in the hot quartile); plus
`rate_nominal` and `rate_bias` when `nominal` is supplied.

## Details

As a rule of thumb, a series is linear enough to simplify when its
deviations (and especially `max_dev_hot`, the worst deviation in the
upper temperature quartile, where the failure rate M(T) is steepest) are
small relative to `log(10) / b`, the number of degrees for a tenfold
change in an exponential M(T): half that distance already multiplies
M(T) by `sqrt(10)`, about threefold.

That distance is not always the z that
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
reports. z is the number of degrees for a tenfold change in the time to
failure, `k log(10) / b`, with k the shape of a Weibull baseline. The
two are the same under the memoryless baseline, where k = 1, and differ
by the factor k under a Weibull one. Judge the deviations against
`log(10) / b`, not against a reported z.

## Examples

``` r
# A logged temperature series close to a 0.3 degree/min ramp, with a small
# mid-run wobble.
path <- data.frame(path_id = "r1", time = 0:20,
                   temp = 20 + 0.3 * (0:20) + c(rep(0, 8), 0.4, 0.5, 0.3, rep(0, 10)))
assess_ramp_linearity(path)
#> # A tibble: 1 × 8
#>   path_id     n T_start  rate r_squared  rmse max_abs_dev max_dev_hot
#>   <chr>   <int>   <dbl> <dbl>     <dbl> <dbl>       <dbl>       <dbl>
#> 1 r1         21    20.1 0.298     0.994 0.143       0.441      0.0487

# With the setpoint rate, to see the controller's rate bias.
assess_ramp_linearity(path, nominal = c(r1 = 0.3))
#> # A tibble: 1 × 10
#>   path_id     n T_start  rate r_squared  rmse max_abs_dev max_dev_hot
#>   <chr>   <int>   <dbl> <dbl>     <dbl> <dbl>       <dbl>       <dbl>
#> 1 r1         21    20.1 0.298     0.994 0.143       0.441      0.0487
#> # ℹ 2 more variables: rate_nominal <dbl>, rate_bias <dbl>
```
