# The cumulative baseline hazard of a given shape

Returns the function H0(tau) for the chosen baseline. Under the
memoryless shape it is elapsed exposure itself; the Weibull,
log-logistic and Gompertz shapes let the hazard rise or fall with how
long the individual has already been exposed. Shape and scale arrive as
per-draw vectors, and elapsed time may be a vector or a matrix whose
rows are draws.

## Usage

``` r
H0_of(baseline)
```

## Arguments

- baseline:

  A character string, a key of
  [h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md):
  `"mem"`, `"weibull"`, `"llogis"` or `"gompertz"`.

## Value

A function of `(tau, shape, scale)` giving the cumulative baseline
hazard. What `shape` and `scale` mean depends on the baseline:

- `"mem"`: both are ignored.

- `"weibull"` and `"llogis"`: `shape` is k and `scale` is s.

- `"gompertz"`: `shape` is the rate constant gamma and `scale` is
  ignored.

## Details

For the Gompertz shape the first argument is the rate constant gamma and
the second is ignored. A negative gamma is a valid value and means a
hazard that falls the longer the individual has been exposed. The
cumulative hazard then levels off: it rises towards 1/\|gamma\| but
never reaches it, instead of growing without limit. Survival therefore
never drops to zero, and the model predicts that some individuals never
fail (a cure fraction, in survival-analysis terms).
[`H0inv_of()`](https://eogmartins.github.io/rTFT/reference/H0inv_of.md)
gives the time by which a chosen fraction of individuals has failed;
when that fraction is never reached it returns `NaN` for the time,
rather than a large number standing in for "never".

## Examples

``` r
# The memoryless baseline is elapsed exposure itself.
H0_of("mem")(tau = c(10, 30, 60), shape = 1, scale = 1)
#> [1] 10 30 60

# A Weibull baseline with shape 2 and scale 60: H0(tau) = (tau / 60)^2.
H0_of("weibull")(tau = 30, shape = 2, scale = 60)
#> [1] 0.25
```
