# Run one task of a bundle on the cluster

Runs one task of a bundle written by
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md):
the fit of one model, or one refit of a model for cross-validation,
without the units of one fold. It is called on a compute node of the
cluster by the bundle's script `task.sh`, which Slurm runs once for each
task of a job array, and is seldom called by hand.

## Usage

``` r
run_hpc_task(
  dir = ".",
  stage = "fit",
  task = Sys.getenv("SLURM_ARRAY_TASK_ID")
)
```

## Arguments

- dir:

  A character string, the folder of the bundle. Defaults to the working
  folder, which is where `task.sh` runs it.

- stage:

  A character string, `"fit"` (the default) for a fit, or `"kfold"` for
  a refit for cross-validation.

- task:

  A single whole number, or a character string holding one, the number
  of the task among those of its stage, as the tables
  `inputs/tasks-fit.csv` and `inputs/tasks-kfold.csv` number them.
  Defaults to the number Slurm gives the array task, in the environment
  variable `SLURM_ARRAY_TASK_ID`.

## Value

The path of the task's result, invisibly. A message says how the task
ended. A task that fails is an error, raised after its reason has been
written into `results/`.

## Details

A fit is made as
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
makes it: with
[`fit_hazard()`](https://eogmartins.github.io/rTFT/reference/fit_hazard.md),
refused if it lost a chain, and cross-validated by importance sampling
with
[`loo_hazard()`](https://eogmartins.github.io/rTFT/reference/loo_hazard.md),
leaving out the unit settled when the bundle was written. A refit starts
from the model's fit, read from the result of its task, so it runs only
once that fit has ended well, and it predicts the units left out as
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
does. Each model is compiled once, by its fit, into a folder of its own
under `stan/` in the bundle, named after its fit task, where its refits
find the compiled program; two models that share their Stan code, as two
whose formulas differ only in the name of a covariate do, are therefore
never compiled into one file at the same time. A refit draws the effects
of the units it leaves out at random, with R's own generator; when the
bundle was written with a `seed`, the refit sets R's seed from it and
from the task's number, so that a refit run again gives the same
estimate.

Before anything is fitted, the files of the bundle are checked against
their checksums and the installed rTFT against the version that wrote
the bundle; a task stopped by the second check is recorded as failed,
with the reason. A task whose result is already in `results/` is not run
again, so the tasks of a stage can be submitted again to run only those
that failed or did not finish.

Each task writes into `results/`: a record of when and where it started
(`fit-001.started`, for the first fit), and when it ends, a record of
how it ended (`fit-001.csv`) and, if it ended well, its result
(`fit-001.rds`). The record of how an earlier attempt ended is moved
aside when a task starts again (`fit-001.previous.csv`), so that a task
running again is seen as started. A result is written under another name
and renamed when complete, so a task stopped at the end of its wall time
leaves no result that looks complete. Every result is stamped with the
checksums of the bundle, so that
[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md)
never assembles it with another bundle.

## See also

[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md),
which writes the bundle, and
[`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md),
which says what each task has done.

## Examples

``` r
if (FALSE) { # \dontrun{
# on the cluster, from the folder of the bundle; task.sh runs the same, with the
# number of the task given by Slurm
run_hpc_task(stage = "fit")

# the second fit of the bundle in the folder my-run, run by hand
run_hpc_task("my-run", stage = "fit", task = 2)
} # }
```
