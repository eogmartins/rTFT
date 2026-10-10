# Draw from the posterior predictive of a `hazard_mixed` observation

Draw from the posterior predictive of a `hazard_mixed` observation

## Usage

``` r
posterior_predict_hazard_mixed(i, prep, ...)
```

## Arguments

- i:

  Observation index.

- prep:

  A `brmsprep` object.

- ...:

  Unused.

## Value

For observation `i`, a vector of zeros and ones, one per draw; brms
binds these columns into the draws-by-observations matrix.
