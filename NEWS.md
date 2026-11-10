# rTFT 0.2.2

The default prior on log a changes for two kinds of model, a model whose
formula of log a has terms that mark groups and a model with a threshold
function, for which it is wider, and it is built with the prior of the user's
own on log b (a normal) or on `Tc` (a normal or a Student-t) where one is given,
rather than with the defaults. A refit of such a model gives different
estimates. Every other model is given the same priors as in 0.2.1.

## Default priors

* The default prior on log a of `build_hazard_prior()` is scaled to the design.
  Where the formula of log a has terms that mark groups (a factor, a character
  or logical column, or a numeric column with two values, alone or in an
  interaction), the standard deviation of its intercept grows with the largest
  distance between the reference temperature `Tr` and the median failure
  temperature of any group, and the coefficients of each group term are given a
  prior whose standard deviation grows with the range of the groups' median
  failure temperatures (a term whose coefficients are the levels' own values,
  as in a formula without an intercept, takes the intercept's prior instead).
  With groups that failed far from `Tr` or from each other, the prior of 0.2.1
  pulled their differences towards zero and biased the shape of the Weibull
  baseline low (by about 0.7 in simulations of three species whose failures lay
  some degrees apart). A formula without such terms keeps the prior of 0.2.1.
* The default prior on log a is carried from the failures to `Tr` with the
  priors on log b and `Tc` that apply to the fit: a normal prior of the user's
  own on log b, and a normal or Student-t prior on `Tc`, are read, where 0.2.1
  used the defaults whatever was given. For a threshold function, the prior on
  the intercept of log a is carried across the range of temperatures the prior
  on `Tc` spans, so it is wider than in 0.2.1 (a standard deviation of 6.6
  rather than 2 for the zero-power threshold model of the vignette *Comparing
  and averaging models*).
* `slope_prior` can be given per parameter, as a named list. A prior named for
  a parameter whose formula has nothing for it to reach warns (for log a, a
  group term whose coefficients are differences between groups counts,
  intercept or not), and anything that is not a prior string, or a named list
  of them, is refused.
* The coefficients of a group term of log a are named as brms names them
  whatever order the formula is written in. An interaction written before one
  of its main effects (`~ 1 + species:do + do`) was named in R's order, which
  brms does not use, so `fit_hazard()` stopped in brms with a prior on a
  coefficient that did not exist.
* The data are refused, each case with its own condition, when they lack a
  column of the formula of log a, have missing values in one, or have a
  categorical column with a single value; these had given an infinite standard
  deviation or an error from R that did not say what was wrong.

## Checking the priors

* `check_prior()` measures how much the priors pull the estimates of a fitted
  model, without refitting it, by the power-scaling sensitivity analysis of
  Kallioinen et al. (2024, <https://doi.org/10.1007/s11222-023-10366-5>) in the
  package priorsense, now suggested at version 1.4.0 or later. It scales each
  block of priors in turn (one per parameter, one per co-stressor, the
  group-level priors, and all of them together) and reports the sensitivity of
  every population-level coefficient and, with `newdata`, of each quantity that
  `derive_tdt_ttl()` reports, sCTmax included when `t_ref` is given. Each
  sensitivity comes with the Pareto k of the power-scaling that produced it
  (`k_prior`, `k_likelihood`), with a warning when a scaling is unreliable; a
  block whose priors are all flat is left out, with a message. A second
  table sets the prior of each coefficient beside its posterior. The help of
  `build_hazard_prior()` and the vignette *Getting started with rTFT* point to
  it.

## Fitting

* The switch of the threshold functions is written in a form whose exponential
  cannot overflow. In the plain form, a temperature more than about 14 degrees
  below `Tc` (at the default `k_sw`) gave the rate an infinite logarithm and
  Stan an undefined gradient, so a chain whose `Tc` started or wandered that far
  above the coolest temperature of the data could not go on. The two forms give
  the same rate wherever the plain one was finite.
* `build_hazard_data()` refuses a `grid_dT` that would cut one individual's
  exposure into more than 100,000 intervals, before the table is built; nothing
  had bounded the rows a very fine grid makes.

## Documentation

* The models of the vignette *Comparing and averaging models* are fitted again
  under 0.2.2. With its wider prior on log a, the zero-power threshold model
  predicts the trials left out better than before (its elpd is 7.5 below that
  of the best model rather than 16.3) and receives a little more weight; the
  averaged sCTmax and dCTmax are almost unchanged.

# rTFT 0.2.1

## Bug fixes

* `loo_hazard()` stopped with "argument is not interpretable as logical" when
  loo could not estimate the Pareto k of a unit, as loo 2.10.0 cannot for a fit
  whose chains did not mix. On a cluster the fit was then lost with its task.
  Such a unit is now counted among those for which PSIS cannot be trusted, in
  the warning of `loo_hazard()` and in `high_k` of the summary that
  `fit_hazard_models()` and `collect_hpc()` return.

# rTFT 0.2.0

## Comparing and averaging models

* `build_hazard_models()` builds a table of models to compare, one row for each
  failure-rate function with a baseline and the formulas of its parameters. A
  formula is given to the models that have its parameter, and a name that is
  misspelt is reported before anything is fitted.
* `fit_hazard_models()` fits every model of the table to the same data, checks
  the sampling of each fit by one rule (Rhat, the effective sample sizes and the
  divergent transitions, met or not in the column `rule_met` of its summary), and
  cross-validates each with the same
  unit left out, by importance sampling and, when asked, by refitting with the
  same folds for every model. A model that fails does not stop the others. It
  returns the fits and the cross-validation results as `hazard_weights()` and
  the `average_*()` functions read them, with a summary of one row per model.
* `hazard_weights()` names a model that gives some unit left out no chance at
  all (a log predictive density of -Inf), as a model with a threshold can for a
  unit whose failures it places below the threshold. Pseudo-BMA gives such a
  model a weight of zero; stacking can still weight it, even the most, when
  other models predict that unit, and the message and the help say why. A unit
  that no model gives any chance stops the weighting with an error that names
  it, where the stacking would otherwise fail without saying why.

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
  failure. When the bundle was written with a `seed`, each refit seeds R from it
  and from the task's number, so that the effects it draws for the units left
  out, and so its estimate, are the same when it is run again.
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
  of the jobs. The data, the unit, the folds and the sampling rule are copied
  from the bundle, never chosen again, as are the place on the cluster the bundle
  is copied to and where rTFT is installed from, and a change that would make
  another model, such as a prior, is refused. A round can name any model of the
  first bundle, so a model left out of one round can be fitted again in a later
  one. `collect_hpc()`, given the folders of a bundle and of its rounds, takes
  each model from its latest round, even when that round failed or has no
  result yet, and says which round in a new column of the summary, `round`.
* The bundle's folder `inputs/` holds the data, the models, the settings and the
  tasks with their checksums; `write_hpc_bundle()` refuses a table in which the
  same model appears under two names, which would be weighted twice; each model
  is compiled into a folder of its own under `stan/`, so that two models with the
  same Stan code are never compiled into one file at the same time; and
  `install.sh` keeps its temporary files out of the login node's `/tmp`, which
  is small and shared.
* A fifth vignette, *Fitting models on a cluster*, takes a user who has never
  used a cluster from an account on the clusters of the Digital Research
  Alliance of Canada to the assembled results of a comparison, with the output
  of real runs on the Fir cluster.

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

* The estimate by refitting is no longer called exact. The effect of each unit
  left out is drawn at random, so two runs without a seed give different
  estimates (they differed by about 1 in elpd in a test on a cluster), and a seed
  gives the same estimate again without making it more precise. The help of
  `kfold_hazard()` and `fit_hazard_models()` says so, and the help of
  `fit_hazard_models()` says where the limits of its sampling rule come from.
* `tools` is imported (for the checksums of a bundle), `cmdstanr` and `withr`
  are suggested, and `Additional_repositories` names the repository of
  `cmdstanr`, which is not on CRAN.
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
