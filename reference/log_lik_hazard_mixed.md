# Pointwise log-likelihood of a `hazard_mixed` observation

The R counterpart of the `hazard_mixed_lpmf` in Stan, reading the
per-row known-time flag from `prep$data$vint1`. It returns the same log
density in R. brms needs this function after fitting to compute the
pointwise log-likelihood (used by LOO and WAIC), as the Stan function
cannot be called from R.

## Usage

``` r
log_lik_hazard_mixed(i, prep)
```

## Arguments

- i:

  Observation index.

- prep:

  A `brmsprep` object.

## Value

For observation `i`, a vector of log-likelihood values, one per draw;
brms binds these columns into the draws-by-observations matrix.
