# log M(T) for each failure rate function, as a brms expression string

Each entry is the expression string of log M(T) for one function of
[M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md),
which brms uses to build the model formula. The strings are written in
terms of the parameter names and the data columns `Ta` (the row's
temperature: the interval midpoint, or the temperature at the moment of
a known-time failure) and `Tr` (the reference temperature of the model;
see
[M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md)).
The threshold functions add a small constant inside a logarithm so it
stays finite at the critical temperature. The constant is written
`10^-12` and not `1e-12`, because brms spaces out the minus sign of a
number in scientific notation when it writes the Stan program
(`1e - 12`), which Stan does not accept. The temperature column is named
`Ta` rather than Kingsolver and Umbanhowar's `T`, because `T` is an
alias for `TRUE` in R and a column of that name would silently shadow it
in the data-handling code. The threshold functions also reference
`k_sw`, the softening steepness, taken from the `k_sw` field of
[M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md)
when the model formula is built by brms. The strings are the single
shared definition that fitting, simulation and reporting all evaluate.

## Usage

``` r
brms_log_M
```

## Format

A named list, one expression string per failure rate function of
[M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md).
