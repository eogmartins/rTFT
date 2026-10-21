# rTFT (development version)

## Comparing and averaging models

* `build_hazard_models()` builds a table of models to compare, one row for each
  failure-rate function with a baseline and the formulas of its parameters. A
  formula is given to the models that have its parameter, and a name that is
  misspelt is reported before anything is fitted.
* `fit_hazard_models()` fits every model of the table to the same data, checks
  each fit for convergence by one rule, and cross-validates each with the same
  unit left out, by importance sampling and, when asked, by refitting with the
  same folds for every model. A model that fails does not stop the others. It
  returns the fits and the cross-validation results as `hazard_weights()` and
  the `average_*()` functions read them, with a summary of one row per model.
* `hazard_weights()` names a model that gives some unit left out no chance at
  all (a log predictive density of -Inf), as a model with a threshold can for a
  unit whose failures it places below the threshold. Pseudo-BMA gives such a
  model a weight of zero; stacking can still weight it, even the most, when
  other models predict that unit, and the message and the help say why.

## Running on a cluster

* `write_hpc_bundle()` writes a table of models, its data and the settings of
  the run into one folder, to be copied to a cluster that schedules its jobs
  with Slurm, with scripts written for the clusters of the Digital Research
  Alliance of Canada. The run is checked, the unit left out in cross-validation
  chosen and the folds assigned before anything is written, and the account,
  time, memory and modules are written once, in the folder's `settings.sh`.
* `check_hpc_bundle()`, run on the cluster by the folder's `check.sh`, checks by
  their checksums that the files arrived unchanged, that the installed rTFT is
  the version that wrote them (and, between releases, from the same commit),
  that cmdstanr finds CmdStan, and that every model can be assembled, before
  anything is submitted.
* `run_hpc_task()` runs one task of a bundle on the cluster: it fits one model,
  as `fit_hazard_models()` fits it, or refits one model without the units of one
  fold, starting from the model's fit. A task whose result is already there is
  not run again, and each task records how it ended, with the reason for a
  failure.
* The bundle's `submit.sh` submits the fits as one job array and the refits for
  cross-validation as another, queued to start when the fits have ended, and
  `task.sh` runs each task. A task that failed or ran out of time can be
  submitted again alone.
* `hpc_status()`, run on the cluster by the bundle's `status.sh`, says what each
  task has done and why a task is not done, including a task that Slurm stopped
  (out of time, out of memory, or for a reason of its own), and gives the
  command that submits again only the tasks not done, leaving out those still
  waiting in the queue or running.
* `collect_hpc()` assembles the results, once they have been copied back, into
  what `fit_hazard_models()` returns. Each result carries the checksums of the
  bundle it was computed from, so it is never assembled with another; models
  and folds without a result are named.
* `next_hpc_round()` writes the next round of a bundle for some of its models,
  typically those that did not converge, with new settings of the sampler and
  of the jobs. The data, the unit, the folds and the convergence rule are
  copied from the bundle, never chosen again, and a change that would make
  another model, such as a prior, is refused. A round can name any model of the
  first bundle, so a model left out of one round can be fitted again in a later
  one. `collect_hpc()`, given the folders of a bundle and of its rounds, takes
  each model from its latest round, even when that round failed or has no
  result yet, and says which round in a new column of the summary, `round`.

## Checking a fitted model

* `hazard_residuals()` draws the residuals several times with
  `n_randomisations`, so that a pattern can be told apart from the luck of one
  draw. Under the same seed the first randomisation is the one a single draw
  gives. The output gains a column, `randomisation`. The help says how far each
  kind of residual can move between draws: not at all for a failure at a known
  time, within its interval for a failure known to an interval, and up to 1 for
  a survivor.
* `plot_residuals()` draws several randomisations faintly, against the uniform
  distribution as a line each.
* `check_survival()` warns when a group holds fewer than five individuals, as
  every group does when each individual has a temperature trajectory of its
  own. In a group that small the band of simulated proportions is wide, often
  from 0 to 1, and a model wrong for every individual can pass; the residuals,
  which pool every individual, are then the check to judge by.

## Bug fixes

* A failure-rate function or a baseline given as a factor was read by its
  position among the factor's levels, not by its text, so another function of
  the catalogue could be fitted without an error. `fit_hazard()` and the other
  functions that take `m_form` or `baseline` now refuse anything but a single
  character string.
* `loo_hazard()`, and `fit_hazard_models()` through it, could hang for good
  where the option `mc.cores` was set above 1, as is common with brms: the loo
  package then copied the R process to work on several cores, and a copy made
  once a sampler's threads were running could wait forever. It now asks the
  loo package for one core; a matrix of one column per unit is quick on one.
* A temperature path with two readings at one time, as when a sudden change is
  written as two readings at the moment of the change, was accepted:
  `build_hazard_data()` replaced the two by their mean and followed a path that
  climbed steadily instead of jumping, with only a warning from base R. A path
  must now have one reading at each time, and a sudden change is written as two
  readings a moment apart.
* Where a temperature path bends within an interval at risk (a logged record, or
  a ramp held and then resumed, checked less often than it was logged),
  `build_hazard_data()` cut the interval into equal parts sized from its whole
  span, so the parts over the steep stretch could change by more than
  `grid_dT`, by over a degree in an example with a short steep climb. A part
  that still changes by more than `grid_dT` is now cut again. Paths that are
  straight between checks, such as those of the vignettes, give the same table
  as before.

## Other changes

* The help of `build_hazard_data()` says what the default time grid costs at a
  high thermal sensitivity (about 1% in b at b = 4.7, from simulations), and
  advises a step of 0.05 where b is above about 3. Its warning for a failure
  interval that spans more than half a degree now says how large the bias can
  be on a ramp checked at long intervals (b 20% low at b = 4.7 with intervals
  of up to 0.9 degrees), and the help advises checking more often.

# rTFT 0.1.1

## Bug fixes

* From version 2.10.0, `loo::loo_compare()` gives the compared models in a
  column rather than as row names. A test and the vignette "Comparing and
  averaging models" read the row names; both now read the models from
  either place. The functions of the package are unchanged.

## Other changes

* The package has a DOI, given by `citation("rTFT")` and the README.
* The package is checked on Linux, macOS and Windows with each push.
* The package can be installed from R-universe.
* The website has a page title and description for search engines.

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
