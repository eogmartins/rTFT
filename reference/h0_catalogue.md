# Baseline hazard functions h0(tau)

This catalogue contains the four baseline hazard functions h0(tau), the
factor in elapsed exposure time that multiplies the failure rate M(T).
`mem` is the memoryless case, in which the hazard does not depend on how
long the individual has already been exposed; `weibull`, `llogis` and
`gompertz` allow it to rise or fall with elapsed exposure.

## Usage

``` r
h0_catalogue
```

## Format

A named list, one entry per baseline function.

## Details

`mem` is the baseline canonical to Kingsolver and Umbanhowar (2018), and
is the default throughout. Their framework has no
baseline-hazard-in-time component at all: the only hazard is M(T)
evaluated along the temperature trajectory T(t). The other three are
extensions of their framework.

`gompertz`, h0(tau) = exp(gamma \* tau), is the only one of the four
baselines whose hazard increases exponentially with elapsed exposure
time, which is also the shape an exponential M(T) exhibits as
temperature increases through a ramp. Because the two share that shape,
a single ramp rate cannot separate them, as it only sees the combined
acceleration, not how much comes from the increasing temperature and how
much from time under exposure. Two or more ramp rates can separate them,
because the temperature part steepens with the ramp rate while the
elapsed time part does not. The same confound occurs between `weibull`,
whose hazard is a power of elapsed time, and a power M(T) (`zeropower`
or `zeropower_threshold`), since a power of temperature along a ramp is
also a power of elapsed time. So the two are confounded under a single
ramp and, as for `gompertz`, are separated by two or more ramp rates.
`llogis`, whose hazard rises then falls, matches none of the
failure-rate functions of temperature, which all rise monotonically, so
it is not confounded with any of them.

A Gompertz baseline hazard is a different object from the Gompertz
distribution of failure times that Kingsolver and Umbanhowar's Eqn 6
derives. The distribution is an outcome of the model: an exponential
M(T) traversed along a linear ramp, under a memoryless baseline,
produces Gompertz-distributed failure times. The baseline hazard
presented here is an assumption about how risk depends on elapsed
exposure time, independent of temperature.

Each entry gives the integer the Stan code dispatches on (`stan_id`) and
the parameters of each baseline function (`dpars`), as in
[M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md).
`mem` has no parameters. `weibull` has a shape parameter and `llogis`
shape and scale parameters, on the log scale (`logshape`, `logscale`).
`gompertz` has the rate constant gamma (`gam`) parameter in reciprocal
time, left on the natural scale and unconstrained in sign (a positive
value is a hazard rising with elapsed exposure time and a negative one a
hazard falling). The amplitude of the Gompertz hazard is not included,
because it is the same quantity as `loga` in M(T), and `mem` omits a
parameter for the same reason.

The scale of the `weibull` baseline is not estimated. It is fixed at the
reference time `tau_ref` set in
[build_hazard_data](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
so that h0(tau) = (k / tau_ref) (tau / tau_ref)^(k - 1) and the
cumulative baseline hazard at `tau_ref` is 1. A Weibull scale s and the
amplitude a of M(T) enter the hazard only through a / s^k, so the data
cannot tell them apart and only one of them can be estimated. With the
scale fixed, a is the cumulative hazard reached over the reference time
at the temperature where M(T) = a (the reference temperature `Tr`, for
the functions anchored there). The scale of `llogis` also enters where
its hazard turns from rising to falling, so it is estimated; the data
inform it only weakly when few failures come after that turn, and its
posterior then leans on its prior.

## References

Kingsolver, J. G. and Umbanhowar, J. (2018). The analysis and
interpretation of critical temperatures. Journal of Experimental Biology
221, jeb167858.
[doi:10.1242/jeb.167858](https://doi.org/10.1242/jeb.167858)
