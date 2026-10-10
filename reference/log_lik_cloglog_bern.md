# Pointwise log-likelihood of a `cloglog_bern` observation

The R counterpart of the `cloglog_bern_lpmf` in Stan. It returns the
same log density in R. brms needs this function after fitting to compute
the pointwise log-likelihood (used by LOO and WAIC), as the Stan
function cannot be called from R.

## Usage

``` r
log_lik_cloglog_bern(i, prep)
```

## Arguments

- i:

  Observation index.

- prep:

  A `brmsprep` object.

## Value

For observation `i`, a vector of log-likelihood values, one per draw;
brms binds these columns into the draws-by-observations matrix.
