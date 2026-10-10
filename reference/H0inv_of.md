# The inverse of the cumulative baseline hazard

Returns the function that gives the elapsed time at which the cumulative
baseline hazard reaches a given value, which is what turns a failure
probability into a thermal death time.

## Usage

``` r
H0inv_of(baseline)
```

## Arguments

- baseline:

  A character string, a key of
  [h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md):
  `"mem"`, `"weibull"`, `"llogis"` or `"gompertz"`.

## Value

A function of `(x, shape, scale)` giving the elapsed time at which the
cumulative baseline hazard equals `x`. `shape` and `scale` mean the same
as in [`H0_of()`](https://eogmartins.github.io/rTFT/reference/H0_of.md):

- `"mem"`: both are ignored.

- `"weibull"` and `"llogis"`: `shape` is k and `scale` is s.

- `"gompertz"`: `shape` is the rate constant gamma and `scale` is
  ignored.

## Details

Under the Gompertz baseline with a negative rate constant the cumulative
hazard levels off below 1/\|gamma\|, so a target above that level is
never reached. The time returned there is `NaN`, rather than a large
number standing in for "never". A target exactly at the level is reached
only after infinite time, and returns `Inf`.

## Examples

``` r
# Invert a Weibull baseline (shape 2, scale 60): the time at which H0 reaches 0.25.
H0inv_of("weibull")(x = 0.25, shape = 2, scale = 60)   # 30
#> [1] 30
```
