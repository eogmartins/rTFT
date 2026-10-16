# rTFT (development version)

## Bug fixes

* A failure-rate function or a baseline given as a factor was read by its
  position among the factor's levels, not by its text, so another function of
  the catalogue could be fitted without an error. `fit_hazard()` and the other
  functions that take `m_form` or `baseline` now refuse anything but a single
  character string.

# rTFT 0.1.0

First release.

## Fitting

* `build_hazard_data()` builds the table of intervals at risk from failure times
  that are known-time, interval-censored or a mixture of both, under any
  temperature trajectory. `ramp_path()` builds the trajectory of a linear ramp,
  and `assess_ramp_linearity()` compares a logged ramp with a straight line.
* `fit_hazard()` fits the model using Bayesian inference in Stan through brms,
  with a choice of functions for the failure rate (`M_catalogue`) and for the
  baseline (`h0_catalogue`). The failure rate may depend on covariates and on
  group-level effects.
* `costressor()` declares a co-stressor (e.g., dissolved oxygen), with its main
  effect and, when asked for, its interaction with temperature. A model may have
  more than one co-stressor.
* `build_hazard_prior()` shows and sets the priors.

## Checking a fitted model

* `check_survival()` and `plot_check()` compare the proportion of individuals
  observed to survive with the proportions simulated from the model.
* `hazard_residuals()` and `plot_residuals()` give one randomised residual per
  individual.

## Derived quantities

* `derive_tdt_ttl()` derives the static critical thermal maximum, the thermal
  death time curves and the thermal tolerance landscape, and `derive_dctmax()`
  the dynamic critical thermal maximum under a ramp, each with its posterior
  uncertainty.
* `plot_tdt()`, `plot_ttl()` and `plot_ttl_3d()` draw them.
* `M_draws()` returns the posterior draws of the failure rate across
  temperatures.

## Comparing and averaging models

* `loo_hazard()` and `kfold_hazard()` compare models by cross-validation that
  leaves out whole individuals or trials.
* `hazard_weights()` gives the weights of the models, and `average_tdt_ttl()`,
  `average_dctmax()` and `average_survival()` average the derived quantities and
  the predictions.

## Predicting survival

* `predict_survival()` predicts survival along a record of temperature (and of
  any co-stressors recorded with it), and reports how much of the record and of
  the predicted risk lies outside the temperatures and exposure times of the
  assay. `plot_survival()` and `plot_risk()` draw the prediction.

## Data

* `static_assay`, `ramp_assay` and `stream_record` are simulated datasets used
  in the help pages and the vignettes.
