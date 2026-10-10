# Complementary log-log Bernoulli response, on the identity link

The response distribution for interval-censored failures. brms's own
`bernoulli(link = "cloglog")` forms the failure probability first, and
that probability rounds to exactly one once the linear predictor passes
about 3.6, after which every survived interval returns minus infinity.
This is not a rare edge case in thermal failure data: near CTmax the
failure rate grows very large, so the cumulative hazard over an interval
(and with it the linear predictor) routinely passes that threshold, and
the fit becomes unusable. This family keeps the predictor on the
identity link and forms the log density from it directly, which is
numerically stable across the whole range.

## Usage

``` r
cloglog_bern_family()
```

## Value

A `customfamily` for use with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).
