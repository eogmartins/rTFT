# rTFT: Hazard Models for Thermal Failure Time Data

Fits discrete-time survival models to the times at which individual
organisms fail under thermal exposure (e.g., loss of equilibrium,
knockdown or death), using Bayesian inference in 'Stan' through 'brms'.
The hazard is the product of a failure rate that depends on temperature,
following Kingsolver and Umbanhowar (2018)
[doi:10.1242/jeb.167858](https://doi.org/10.1242/jeb.167858) , and a
baseline that may depend on the time under exposure. Failure times may
be exact, interval-censored or both, and may come from static assays,
ramping assays or both. The failure rate may depend on covariates,
group-level effects and other stressors such as dissolved oxygen. Static
and dynamic critical thermal maxima, thermal death time curves and
thermal tolerance landscapes are derived from the fitted model with
their posterior uncertainty. Models can be compared by cross-validation
that leaves out whole individuals or trials, averaged, and used to
predict survival along a record of temperature and of any co-stressors
recorded with it.

## See also

Useful links:

- <https://eogmartins.github.io/rTFT/>

- <https://github.com/eogmartins/rTFT>

- Report bugs at <https://github.com/eogmartins/rTFT/issues>

## Author

**Maintainer**: Eduardo G. Martins <eduardo.martins@unbc.ca>
([ORCID](https://orcid.org/0000-0003-3500-9331))

Authors:

- Eduardo G. Martins <eduardo.martins@unbc.ca>
  ([ORCID](https://orcid.org/0000-0003-3500-9331))
