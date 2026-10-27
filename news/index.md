# Changelog

## rTFT 0.2.1

### Bug fixes

- [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  stopped with “argument is not interpretable as logical” when loo could
  not estimate the Pareto k of a unit, as loo 2.10.0 cannot for a fit
  whose chains did not mix. On a cluster the fit was then lost with its
  task. Such a unit is now counted among those for which PSIS cannot be
  trusted, in the warning of
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  and in `high_k` of the summary that
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  and
  [`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md)
  return.

## rTFT 0.2.0

### Comparing and averaging models

- [`build_hazard_models()`](https://eogmartins.github.io/rTFT/reference/build_hazard_models.md)
  builds a table of models to compare, one row for each failure-rate
  function with a baseline and the formulas of its parameters. A formula
  is given to the models that have its parameter, and a name that is
  misspelt is reported before anything is fitted.
- [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  fits every model of the table to the same data, checks the sampling of
  each fit by one rule (Rhat, the effective sample sizes and the
  divergent transitions, met or not in the column `rule_met` of its
  summary), and cross-validates each with the same unit left out, by
  importance sampling and, when asked, by refitting with the same folds
  for every model. A model that fails does not stop the others. It
  returns the fits and the cross-validation results as
  [`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
  and the `average_*()` functions read them, with a summary of one row
  per model.
- [`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
  names a model that gives some unit left out no chance at all (a log
  predictive density of -Inf), as a model with a threshold can for a
  unit whose failures it places below the threshold. Pseudo-BMA gives
  such a model a weight of zero; stacking can still weight it, even the
  most, when other models predict that unit, and the message and the
  help say why. A unit that no model gives any chance stops the
  weighting with an error that names it, where the stacking would
  otherwise fail without saying why.

### Running on a cluster

- [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
  writes a table of models, its data and the settings of the run into
  one folder, to be copied to a cluster that schedules its jobs with
  Slurm, with scripts written for the clusters of the Digital Research
  Alliance of Canada. The run is checked, the unit left out in
  cross-validation chosen and the folds assigned before anything is
  written, and the account, time, memory and modules are written once,
  in the folder’s `settings.sh`.
- [`check_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/check_hpc_bundle.md),
  run on the cluster by the folder’s `check.sh`, checks by their
  checksums that the files arrived unchanged, that the installed rTFT is
  the version that wrote them (and, between releases, from the same
  commit), that cmdstanr finds CmdStan, and that every model can be
  assembled, before anything is submitted.
- [`run_hpc_task()`](https://eogmartins.github.io/rTFT/reference/run_hpc_task.md)
  runs one task of a bundle on the cluster: it fits one model, as
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  fits it, or refits one model without the units of one fold, starting
  from the model’s fit. A task whose result is already there is not run
  again, and each task records how it ended, with the reason for a
  failure. When the bundle was written with a `seed`, each refit seeds R
  from it and from the task’s number, so that the effects it draws for
  the units left out, and so its estimate, are the same when it is run
  again.
- The bundle’s `submit.sh` submits the fits as one job array and the
  refits for cross-validation as another, queued to start when the fits
  have ended, and `task.sh` runs each task. A task that failed or ran
  out of time can be submitted again alone.
- [`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md),
  run on the cluster by the bundle’s `status.sh`, says what each task
  has done and why a task is not done, including a task that Slurm
  stopped (out of time, out of memory, or for a reason of its own), and
  gives the command that submits again only the tasks not done, leaving
  out those still waiting in the queue or running.
- [`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md)
  assembles the results, once they have been copied back, into what
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  returns. Each result carries the checksums of the bundle it was
  computed from, so it is never assembled with another; models and folds
  without a result are named.
- [`next_hpc_round()`](https://eogmartins.github.io/rTFT/reference/next_hpc_round.md)
  writes the next round of a bundle for some of its models, typically
  those that did not converge, with new settings of the sampler and of
  the jobs. The data, the unit, the folds and the sampling rule are
  copied from the bundle, never chosen again, as are the place on the
  cluster the bundle is copied to and where rTFT is installed from, and
  a change that would make another model, such as a prior, is refused. A
  round can name any model of the first bundle, so a model left out of
  one round can be fitted again in a later one.
  [`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md),
  given the folders of a bundle and of its rounds, takes each model from
  its latest round, even when that round failed or has no result yet,
  and says which round in a new column of the summary, `round`.
- The bundle’s folder `inputs/` holds the data, the models, the settings
  and the tasks with their checksums;
  [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
  refuses a table in which the same model appears under two names, which
  would be weighted twice; each model is compiled into a folder of its
  own under `stan/`, so that two models with the same Stan code are
  never compiled into one file at the same time; and `install.sh` keeps
  its temporary files out of the login node’s `/tmp`, which is small and
  shared.
- A fifth vignette, *Fitting models on a cluster*, takes a user who has
  never used a cluster from an account on the clusters of the Digital
  Research Alliance of Canada to the assembled results of a comparison,
  with the output of real runs on the Fir cluster.

### Checking a fitted model

- [`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
  draws the residuals several times with `n_randomisations`, so that a
  pattern can be told apart from the luck of one draw. Under the same
  seed the first randomisation is the one a single draw gives. The
  output gains a column, `randomisation`. The help says how far each
  kind of residual can move between draws: not at all for a failure at a
  known time, within its interval for a failure known to an interval,
  and up to 1 for a survivor.
- [`plot_residuals()`](https://eogmartins.github.io/rTFT/reference/plot_residuals.md)
  draws several randomisations faintly, against the uniform distribution
  as a line each.
- [`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md)
  warns when a group holds fewer than five individuals, as every group
  does when each individual has a temperature trajectory of its own. In
  a group that small the band of simulated proportions is wide, often
  from 0 to 1, and a model wrong for every individual can pass; the
  residuals, which pool every individual, are then the check to judge
  by.

### Bug fixes

- A failure-rate function or a baseline given as a factor was read by
  its position among the factor’s levels, not by its text, so another
  function of the catalogue could be fitted without an error.
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  and the other functions that take `m_form` or `baseline` now refuse
  anything but a single character string.
- [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md),
  and
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  through it, could hang for good where the option `mc.cores` was set
  above 1, as is common with brms: the loo package then copied the R
  process to work on several cores, and a copy made once a sampler’s
  threads were running could wait forever. It now asks the loo package
  for one core; a matrix of one column per unit is quick on one.
- A temperature path with two readings at one time, as when a sudden
  change is written as two readings at the moment of the change, was
  accepted:
  [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  replaced the two by their mean and followed a path that climbed
  steadily instead of jumping, with only a warning from base R. A path
  must now have one reading at each time, and a sudden change is written
  as two readings a moment apart.
- Where a temperature path bends within an interval at risk (a logged
  record, or a ramp held and then resumed, checked less often than it
  was logged),
  [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  cut the interval into equal parts sized from its whole span, so the
  parts over the steep stretch could change by more than `grid_dT`, by
  over a degree in an example with a short steep climb. A part that
  still changes by more than `grid_dT` is now cut again. Paths that are
  straight between checks, such as those of the vignettes, give the same
  table as before.

### Other changes

- The estimate by refitting is no longer called exact. The effect of
  each unit left out is drawn at random, so two runs without a seed give
  different estimates (they differed by about 1 in elpd in a test on a
  cluster), and a seed gives the same estimate again without making it
  more precise. The help of
  [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
  and
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  says so, and the help of
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
  says where the limits of its sampling rule come from.
- `tools` is imported (for the checksums of a bundle), `cmdstanr` and
  `withr` are suggested, and `Additional_repositories` names the
  repository of `cmdstanr`, which is not on CRAN.
- The help of
  [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  says what the default time grid costs at a high thermal sensitivity
  (about 1% in b at b = 4.7, from simulations), and advises a step of
  0.05 where b is above about 3. Its warning for a failure interval that
  spans more than half a degree now says how large the bias can be on a
  ramp checked at long intervals (b 20% low at b = 4.7 with intervals of
  up to 0.9 degrees), and the help advises checking more often.

## rTFT 0.1.1

### Bug fixes

- From version 2.10.0,
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
  gives the compared models in a column rather than as row names. A test
  and the vignette “Comparing and averaging models” read the row names;
  both now read the models from either place. The functions of the
  package are unchanged.

### Other changes

- The package has a DOI, given by `citation("rTFT")` and the README.
- The package is checked on Linux, macOS and Windows with each push.
- The package can be installed from R-universe.
- The website has a page title and description for search engines.

## rTFT 0.1.0

First release.

### Fitting

- [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md)
  builds the table of intervals at risk from failure times that are
  known-time, interval-censored or a mixture of both, under any
  temperature trajectory.
  [`ramp_path()`](https://eogmartins.github.io/rTFT/reference/ramp_path.md)
  builds the trajectory of a linear ramp, and
  [`assess_ramp_linearity()`](https://eogmartins.github.io/rTFT/reference/assess_ramp_linearity.md)
  compares a logged ramp with a straight line.
- [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  fits the model using Bayesian inference in Stan through brms, with a
  choice of functions for the failure rate (`M_catalogue`) and for the
  baseline (`h0_catalogue`). The failure rate may depend on covariates
  and on group-level effects.
- [`costressor()`](https://eogmartins.github.io/rTFT/reference/costressor.md)
  declares a co-stressor (e.g., dissolved oxygen), with its main effect
  and, when asked for, its interaction with temperature. A model may
  have more than one co-stressor.
- [`build_hazard_prior()`](https://eogmartins.github.io/rTFT/reference/build_hazard_prior.md)
  shows and sets the priors.

### Checking a fitted model

- [`check_survival()`](https://eogmartins.github.io/rTFT/reference/check_survival.md)
  and
  [`plot_check()`](https://eogmartins.github.io/rTFT/reference/plot_check.md)
  compare the proportion of individuals observed to survive with the
  proportions simulated from the model.
- [`hazard_residuals()`](https://eogmartins.github.io/rTFT/reference/hazard_residuals.md)
  and
  [`plot_residuals()`](https://eogmartins.github.io/rTFT/reference/plot_residuals.md)
  give one randomised residual per individual.

### Derived quantities

- [`derive_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/derive_tdt_ttl.md)
  derives the static critical thermal maximum, the thermal death time
  curves and the thermal tolerance landscape, and
  [`derive_dctmax()`](https://eogmartins.github.io/rTFT/reference/derive_dctmax.md)
  the dynamic critical thermal maximum under a ramp, each with its
  posterior uncertainty.
- [`plot_tdt()`](https://eogmartins.github.io/rTFT/reference/plot_tdt.md),
  [`plot_ttl()`](https://eogmartins.github.io/rTFT/reference/plot_ttl.md)
  and
  [`plot_ttl_3d()`](https://eogmartins.github.io/rTFT/reference/plot_ttl_3d.md)
  draw them.
- [`M_draws()`](https://eogmartins.github.io/rTFT/reference/M_draws.md)
  returns the posterior draws of the failure rate across temperatures.

### Comparing and averaging models

- [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  and
  [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
  compare models by cross-validation that leaves out whole individuals
  or trials.
- [`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
  gives the weights of the models, and
  [`average_tdt_ttl()`](https://eogmartins.github.io/rTFT/reference/average_tdt_ttl.md),
  [`average_dctmax()`](https://eogmartins.github.io/rTFT/reference/average_dctmax.md)
  and
  [`average_survival()`](https://eogmartins.github.io/rTFT/reference/average_survival.md)
  average the derived quantities and the predictions.

### Predicting survival

- [`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
  predicts survival along a record of temperature (and of any
  co-stressors recorded with it), and reports how much of the record and
  of the predicted risk lies outside the temperatures and exposure times
  of the assay.
  [`plot_survival()`](https://eogmartins.github.io/rTFT/reference/plot_survival.md)
  and
  [`plot_risk()`](https://eogmartins.github.io/rTFT/reference/plot_risk.md)
  draw the prediction.

### Data

- `static_assay`, `ramp_assay` and `stream_record` are simulated
  datasets used in the help pages and the vignettes.
