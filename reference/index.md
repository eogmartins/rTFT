# Package index

## Fitting

Building the table of intervals at risk from the failure times and the
temperature of each trial, and fitting the model.

- [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  : Build the row-per-interval table a failure-rate model is fitted to
- [`ramp_path()`](https://eogmartins.github.io/rTFT/reference/ramp_path.md)
  : Build a linear ramp temperature trajectory
- [`assess_ramp_linearity()`](https://eogmartins.github.io/rTFT/reference/assess_ramp_linearity.md)
  : Assess how linear a logged temperature trace is
- [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  : Fit a failure-rate hazard model with brms
- [`costressor()`](https://eogmartins.github.io/rTFT/reference/costressor.md)
  : Declare a co-stressor for a failure-rate hazard model
- [`build_hazard_prior()`](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)
  : Weakly-informative default priors for a failure-rate hazard model
- [`M_catalogue`](https://eogmartins.github.io/rTFT/reference/M_catalogue.md)
  : Failure rate functions M(T)
- [`h0_catalogue`](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md)
  : Baseline hazard functions h0(tau)

## Checking a fitted model

- [`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md)
  : Check a fitted model against the survival observed in each group
- [`plot_check()`](https://eogmartins.github.io/rTFT/reference/plot_check.md)
  : Plot the check of a fitted model against observed survival
- [`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
  : Residuals of a fitted model, one per individual
- [`plot_residuals()`](https://eogmartins.github.io/rTFT/reference/plot_residuals.md)
  : Plot the residuals of a fitted model

## Derived quantities

Critical thermal maxima, thermal death time curves and thermal tolerance
landscapes, each with its posterior uncertainty.

- [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
  : Thermal tolerance derived from a fitted failure-rate model
- [`derive_dctmax()`](https://eogmartins.github.io/rTFT/reference/derive_dctmax.md)
  : Dynamic critical thermal maximum under a ramp, from a fitted model
- [`plot_tdt()`](https://eogmartins.github.io/rTFT/reference/plot_tdt.md)
  : Plot the thermal death time curves
- [`plot_ttl()`](https://eogmartins.github.io/rTFT/reference/plot_ttl.md)
  : Plot the thermal tolerance landscape
- [`plot_ttl_3d()`](https://eogmartins.github.io/rTFT/reference/plot_ttl_3d.md)
  : Plot the thermal tolerance landscape as a surface to rotate
- [`M_draws()`](https://eogmartins.github.io/rTFT/reference/M_draws.md)
  : The fitted failure rate over temperature, for each posterior draw

## Comparing and averaging models

- [`build_hazard_models()`](https://eogmartins.github.io/rTFT/reference/build_hazard_models.md)
  : A table of failure-rate models to fit and compare
- [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  : Fit, check and cross-validate the models of a table
- [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  : Leave-one-out cross-validation that leaves out whole individuals or
  groups
- [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
  : Cross-validation by refitting, leaving out whole individuals or
  groups
- [`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
  : Weights for averaging failure-rate models
- [`average_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/average_tdt_ttl.md)
  : Thermal tolerance averaged over several fitted models
- [`average_dctmax()`](https://eogmartins.github.io/rTFT/reference/average_dctmax.md)
  : Dynamic critical thermal maximum averaged over several fitted models
- [`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md)
  : Survival along a record of temperature, averaged over several fitted
  models

## Running on a cluster

Fitting a table of models on a cluster that schedules its jobs with
Slurm.

- [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
  : Write a table of models into a folder to fit on a cluster
- [`check_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/check_hpc_bundle.md)
  : Check a bundle and the R that will fit it, before anything is
  submitted
- [`run_hpc_task()`](https://eogmartins.github.io/rTFT/reference/run_hpc_task.md)
  : Run one task of a bundle on the cluster
- [`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md)
  : What each task of a bundle has done, and what to do next
- [`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md)
  : Assemble the results of a bundle run on a cluster
- [`next_hpc_round()`](https://eogmartins.github.io/rTFT/reference/next_hpc_round.md)
  : Write the next round of a bundle, for some of its models

## Predicting survival

- [`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
  : Survival along a record of temperature, from a fitted failure-rate
  model
- [`plot_survival()`](https://eogmartins.github.io/rTFT/reference/plot_survival.md)
  : Plot survival along a record of temperature
- [`plot_risk()`](https://eogmartins.github.io/rTFT/reference/plot_risk.md)
  : Plot each period's share of the risk along a record

## Data

Simulated datasets used in the help pages and the vignettes.

- [`static_assay`](https://eogmartins.github.io/rTFT/reference/static_assay.md)
  : A simulated static assay of thermal tolerance
- [`ramp_assay`](https://eogmartins.github.io/rTFT/reference/ramp_assay.md)
  : A simulated ramping assay of thermal tolerance, with two levels of
  oxygen
- [`stream_record`](https://eogmartins.github.io/rTFT/reference/stream_record.md)
  : A simulated record of stream temperature

## Building blocks

The pieces that
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
and the derived quantities are assembled from. Most analyses do not call
them directly.

- [`build_hazard_spec()`](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md)
  : Assemble the brms formula and family for a failure-rate hazard model
- [`brms_log_M`](https://eogmartins.github.io/rTFT/reference/brms_log_M.md)
  : log M(T) for each failure rate function, as a brms expression string
- [`brms_log_h0`](https://eogmartins.github.io/rTFT/reference/brms_log_h0.md)
  : log h0(tau) for each baseline function, as a brms expression string
- [`H0_of()`](https://eogmartins.github.io/rTFT/reference/H0_of.md) :
  The cumulative baseline hazard of a given baseline function
- [`H0inv_of()`](https://eogmartins.github.io/rTFT/reference/H0inv_of.md)
  : The inverse of the cumulative baseline hazard
- [`cloglog_bern_family()`](https://eogmartins.github.io/rTFT/reference/cloglog_bern_family.md)
  : Complementary log-log Bernoulli response, on the identity link
- [`hazard_mixed_family()`](https://eogmartins.github.io/rTFT/reference/hazard_mixed_family.md)
  : Response distribution for a table mixing interval-censored and
  known-time failures

## The package

- [`rTFT`](https://eogmartins.github.io/rTFT/reference/rTFT-package.md)
  [`rTFT-package`](https://eogmartins.github.io/rTFT/reference/rTFT-package.md)
  : rTFT: Hazard Models for Thermal Failure Time Data
