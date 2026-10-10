# Failure rate functions M(T)

This catalogue contains every function the failure rate may take: the
nine failure rate functions M(T) of Kingsolver and Umbanhowar (2018,
Table 1), each catalogue key named for its model there. Without a
threshold: `const` (a constant rate), `exp` (exponential), `zeroexp`
(exponential minus its value at the reference temperature `Tr`, so the
rate is zero there), `power` (a power of temperature above a baseline
rate) and `zeropower` (the same power with no baseline, zero at `Tr`).
The threshold variants `const_threshold`, `exp_threshold`,
`zeroexp_threshold` and `zeropower_threshold` hold the rate at zero
below a fitted critical temperature Tc; above Tc it is constant for
`const_threshold` and accelerates for the other three. The exponential
functions give a straight thermal death time curve, whereas the power
functions give a curved one.

## Usage

``` r
M_catalogue
```

## Format

A named list, one entry per failure rate function M(T).

## Details

`Tr` (`T0` in Kingsolver and Umbanhowar's notation) is a single
reference temperature set for the whole model, not the per-trial
starting temperature. The exponential is anchored at `Tr` and the zero
forms vanish there, but because the failure rate must depend only on the
current temperature, `Tr` cannot be each trial's own starting
temperature. Trials that begin at different temperatures (acclimation
groups, field conditions) enter through their temperature trajectories,
not through `Tr`.

To keep fitting, simulation and reporting from diverging, all three
evaluate the same expression strings, which are held in
[brms_log_M](https://eogmartins.github.io/rTFT/reference/brms_log_M.md).

Each entry gives the integer the Stan code dispatches on (`stan_id`),
the parameters the function needs (`dpars`), and the steepness of the
threshold switch (`k_sw`). The rate parameters a, b and c are
represented on the log scale (`loga`, `logb`, `logc`), so the linear
predictor is unconstrained and their strictly positive values are
recovered by exponentiation; the critical temperature Tc is left on the
natural scale, since it is not sign-constrained. `k_sw` is roughly the
inverse of the width of the transition in temperature, so a larger value
is a sharper switch over a narrower band. The switch is smoothed rather
than hard because Stan fits the model with Hamiltonian Monte Carlo,
which explores the posterior by following the gradient of the
log-density with respect to each parameter; the gradient of a hard step
with respect to the critical temperature is zero on either side of it
and undefined at the step itself, so the sampler gets no slope to move
that temperature along, whereas a smooth switch supplies one.

## The units of M(T)

The hazard is h0(tau) M(T), and which of the two factors carries the
unit of reciprocal time depends on the baseline (see
[h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md)).
With the memoryless and Gompertz baselines h0(tau) has no unit, so M(T)
is a rate per unit time, and a is the failure rate at the temperature
where M(T) = a. With the Weibull and log-logistic baselines h0(tau) is
itself per unit time, so M(T) is a multiplier with no unit; with the
Weibull, a is the cumulative hazard reached over the reference time at
that temperature. The value of a therefore cannot be compared between
models with different baselines, nor between datasets timed in different
units when the baseline is memoryless or Gompertz.

Whatever the baseline, survival at a constant temperature is exp(-M(T)
H0(tau)), with H0 the cumulative baseline hazard.
[`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
derives the critical thermal maximum, the thermal death time curves and
the thermal tolerance landscape from that expression, and those can be
compared between models.

## References

Kingsolver, J. G. and Umbanhowar, J. (2018). The analysis and
interpretation of critical temperatures. Journal of Experimental Biology
221, jeb167858.
[doi:10.1242/jeb.167858](https://doi.org/10.1242/jeb.167858)
