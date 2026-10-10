# Response distribution for a table mixing interval-censored and known-time failures

Some failures are known only to have fallen within an interval, others
at a known time. This family fits both kinds in one model. A row where
the individual survived its interval contributes the same amount either
way, so only the row recording the failure differs: an interval-censored
failure contributes the probability of failing somewhere inside the
interval, whereas a failure seen at a known time contributes the density
at that recorded time. A per-row flag, passed to brms through `vint()`,
tells the model which of the two to use, so a single fit can contain
individuals whose failures were recorded as either kind.

## Usage

``` r
hazard_mixed_family()
```

## Value

A `customfamily` for use with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).
