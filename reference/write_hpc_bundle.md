# Write a table of models into a folder to fit on a cluster

Writes what is needed to fit the models of a table built with
[`build_hazard_models()`](https://eogmartins.github.io/rTFT/reference/build_hazard_models.md)
on a cluster that schedules its jobs with Slurm, such as those of the
Digital Research Alliance of Canada, into one folder, the bundle. The
bundle is then copied to the cluster, where the scripts it contains
install what it needs (`install.sh`) and check it (`check.sh`, which
calls
[`check_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/check_hpc_bundle.md)).

## Usage

``` r
write_hpc_bundle(
  models,
  long,
  dir,
  account,
  observation = c("interval", "continuous", "mixed"),
  unit = NULL,
  kfold = FALSE,
  K = NULL,
  rhat = 1.01,
  ess = 400,
  divergent = 0.001,
  time = "12:00:00",
  mem_per_cpu = "4000M",
  email = NULL,
  modules = c("StdEnv/2023", "r/4.5.0"),
  r_library = NULL,
  remote = NULL,
  rtft_source = NULL,
  template = "alliance",
  ...
)
```

## Arguments

- models:

  A data frame, the table of models from
  [`build_hazard_models()`](https://eogmartins.github.io/rTFT/reference/build_hazard_models.md).

- long:

  A data frame, the table of intervals at risk from
  [`build_hazard_data()`](https://eogmartins.github.io/rTFT/reference/build_hazard_data.md),
  to which every model is fitted.

- dir:

  A character string, the folder to write the bundle into. It is created
  if it does not exist, and must be empty if it does.

- account:

  A character string, the allocation the jobs are charged to, as Slurm's
  `--account` takes it (for example `"def-someone"`).

- observation:

  A character string, how the failure times were observed, as in
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md).

- unit:

  An optional character string, the unit to leave out in
  cross-validation, as in
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md).
  `NULL` (the default) chooses it as
  [`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md)
  does.

- kfold:

  A single logical value. `TRUE` also cross-validates each model by
  refitting, one task per model and fold. Defaults to `FALSE`.

- K:

  An optional single whole number, the number of folds when `kfold` is
  `TRUE`, as in
  [`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md).
  `NULL` (the default) gives one fold per unit. With fewer folds than
  units the units are assigned to folds at random, once; call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) first to get the
  same folds again.

- rhat, ess, divergent:

  The limits of the sampling rule, as in
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md).

- time:

  A character string, the time each task may run, as Slurm's `--time`
  takes it: hours, minutes and seconds (`"12:00:00"`, the default), or
  days first (`"1-12:00:00"`). A task that runs past it is stopped and
  its work is lost. The allocation is charged for the time a task ran,
  times its cores, and not for the time asked; a longer time can mean a
  longer wait in the queue.

- mem_per_cpu:

  A character string, the memory for each core, with its unit, as
  Slurm's `--mem-per-cpu` takes it. Defaults to `"4000M"`, about the
  memory per core of the Alliance's standard nodes, so asking it costs
  no more than the cores do; asking more is charged as more cores. Slurm
  reads `"4G"` as 4096 M, which the clusters of the Alliance point out
  at each submission.

- email:

  An optional character string, the address Slurm writes to when a task
  fails. `NULL` (the default) sends no mail.

- modules:

  A character vector, the modules that give R on the cluster. The
  default, `c("StdEnv/2023", "r/4.5.0")`, is for the clusters of the
  Alliance; `module avail r` on the cluster lists the versions it
  offers.

- r_library:

  An optional character string, the R library on the cluster that holds
  rTFT, brms and cmdstanr. `NULL` (the default) uses R's own personal
  library.

- remote:

  An optional character string, where on the cluster to copy the bundle,
  as `"user@host:folder"` (for example
  `"someone@nibi.alliancecan.ca:scratch/my-run"`), without spaces. A
  folder written without a leading `/` is in your home folder on the
  cluster. It is used only to fill in the commands of `README.txt` and
  of the message, and nothing is sent to the cluster; `NULL` (the
  default) leaves placeholders there. The copy command creates the
  folder, with any folders above it that do not exist yet.

- rtft_source:

  An optional character string, where the cluster installs rTFT from, as
  `remotes::install_github()` takes it. `NULL` (the default) uses the
  commit the installed rTFT came from, when it was installed from
  GitHub, or the tag of its release; a development version installed
  from a copy on your machine must be given, such as
  `"eogmartins/rTFT@dev-0.2.0"`.

- template:

  A character string, `"alliance"` (the default) for the scripts written
  for the clusters of the Alliance, or the path of a folder of your own
  holding `settings.sh`, `install.sh`, `install.R`, `check.sh`,
  `submit.sh`, `task.sh`, `status.sh` and `README.txt`, with the same
  placeholders.
  [`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md)
  finds the log of a task as `logs/<stage>-<array job>_<task>.out`; a
  template that names its logs otherwise loses only the reason given for
  a task stopped by Slurm.

- ...:

  Further arguments for every fit, as in
  [`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md),
  such as `chains`, `iter`, `warmup`, `seed` or `control`, or the priors
  that
  [`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md)
  takes.

## Value

The path of the bundle, invisibly. A message says what was written and
what to do next.

## Details

Everything about the run is settled on the user's own machine, before
anything is copied: the arguments are checked as
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
checks them, the unit left out in cross-validation is chosen and, with
`kfold = TRUE`, the units are assigned to folds once for all the models.
A mistake therefore is quickly found here and not on the cluster hours
later.

The bundle holds:

- `settings.sh`, the account, time, memory, cores and modules, written
  once and read by every script of the bundle. It is the one file to
  edit when a task needs more time or memory.

- `install.sh` and `install.R`, which install what the bundle needs into
  the R of the cluster: rTFT at the version that wrote the bundle, from
  `rtft_source`, cmdstanr, and CmdStan. It is run once on each cluster,
  on the login node or, with `bash install.sh --job`, as a Slurm job.

- `check.sh`, which checks the bundle and the R of the cluster on the
  login node.

- `submit.sh`, which submits the tasks to Slurm as job arrays, one for
  the fits and one for the refits, queued to start when the fits have
  ended; `task.sh`, which Slurm runs for each task and which calls
  [`run_hpc_task()`](https://eogmartins.github.io/rTFT/reference/run_hpc_task.md);
  and `status.sh`, which calls
  [`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md)
  to say what each task has done.

- `README.txt`, the steps to take, with the commands filled in.

- `inputs/`, the data (`data.rds`), the models (`models.rds`), the
  settings of the run (`settings.rds`), the tasks as tables to read
  (`tasks-fit.csv`, and `tasks-kfold.csv` with `kfold = TRUE`), and
  `checksums.txt`, the checksum of each of these files, against which
  the cluster checks what it received.

- `stan/`, `results/` and `logs/`, empty, for the compiled models, the
  results and what each task prints.

The fits on the cluster use the cmdstanr backend, and each task runs its
chains on as many cores as there are chains, so `backend` and `cores`
cannot be given in `...`.

## See also

[`check_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/check_hpc_bundle.md),
which checks the bundle on the cluster,
[`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md),
which says what its tasks have done, and
[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md),
which assembles their results.

## Examples

``` r
# a small ramping assay in two trials (see fit_hazard())
path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
ind    <- data.frame(id = 1:20, path_id = "A", event = 1,
                     trial = rep(c("a", "b"), each = 10),
                     t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
long   <- build_hazard_data(ind, path, "interval", covariates = "trial",
                            checks = checks)
#> ℹ The intervals at risk were divided further, into 2570 rows from 870, so that
#>   temperature changes by at most 0.1 degrees within each.
#> ℹ The hazard is taken as constant within an interval, which biases the fit when
#>   temperature changes much within it. Failures at a known time are also given a
#>   short closing interval, so their hazard is that of the moment of failure.
#> ℹ Set `grid_dT` to change the step, or to `Inf` to keep the intervals as given.
models <- build_hazard_models("exp", c("mem", "weibull"),
                              dpar_formulas = list(loga = ~ 1 + (1 | trial)))

# the bundle, written into a temporary folder here; rtft_source says the cluster
# installs rTFT from its repository on GitHub (by default it is read from how rTFT
# was installed on this machine)
bundle <- file.path(tempdir(), "my-run")
write_hpc_bundle(models, long, bundle, account = "def-someone", kfold = TRUE,
                 remote = "someone@nibi.alliancecan.ca:scratch/my-run",
                 rtft_source = "eogmartins/rTFT", chains = 4, iter = 2000, seed = 1)
#> ℹ Leaving out one trial at a time: the individuals of each trial share its
#>   effect in the model, so leaving out one individual would leave that effect in
#>   the fit.
#> ℹ Use `unit = "id"` to leave out individuals instead.
#> ✔ Wrote a bundle of 2 models to /tmp/RtmpY6cfoi/my-run: 2 fits and 4 refits for
#>   cross-validation, leaving out one trial at a time.
#> ℹ Copy it to the cluster with the command in step 1 of README.txt, in that
#>   folder.
#> ℹ Then, on the cluster, in that folder, run these two (the first once on each
#>   cluster):
#>   bash install.sh
#>   bash check.sh
list.files(bundle, recursive = TRUE)
#>  [1] "README.txt"             "check.sh"               "inputs/checksums.txt"  
#>  [4] "inputs/data.rds"        "inputs/models.rds"      "inputs/settings.rds"   
#>  [7] "inputs/tasks-fit.csv"   "inputs/tasks-kfold.csv" "install.R"             
#> [10] "install.sh"             "settings.sh"            "status.sh"             
#> [13] "submit.sh"              "task.sh"               
```
