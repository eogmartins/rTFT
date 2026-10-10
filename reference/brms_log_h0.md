# log h0(tau) for each baseline function, as a brms expression string

Each entry is the expression string of log h0(tau) for one function of
[h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md),
which brms uses to build the model formula. The strings are written in
terms of the parameter names and the data columns `log_tau` (the log of
the row's elapsed exposure time), for `weibull` `tau_ref` (the reference
time its scale is fixed at) and, for `gompertz`, `tau` itself (the
elapsed exposure time: the midpoint for an interval row, or the time of
a known-time failure). `mem` is the empty string, because a memoryless
baseline contributes nothing: h0(tau) = 1 and the hazard is M(T) alone.

## Usage

``` r
brms_log_h0
```

## Format

A named list, one expression string per baseline function of
[h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md).

## Details

Each string is the baseline **hazard** on the log scale, not the
density. For `llogis` the two differ by one `log1p` term: the hazard has
a single `- log1p(...)`, the density would have two.
