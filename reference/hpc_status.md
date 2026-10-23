# What each task of a bundle has done, and what to do next

Reads the records that the tasks of a bundle written by
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
have left in its folder `results/`, and says, for each stage of the run,
how many tasks are done, which are not and why, and what to do next. It
is run on the cluster by the bundle's script `status.sh`, and works as
well on a copy of the bundle brought back.

## Usage

``` r
hpc_status(dir = ".")
```

## Arguments

- dir:

  A character string, the folder of the bundle. Defaults to the working
  folder, which is where `status.sh` runs it.

## Value

A data frame with one row per task, invisibly: its `stage` (`"fit"` or
`"kfold"`), `task` number, `model`, `fold` (`NA` for a fit), `status`,
the `reason` a task is not done, and the `seconds` it took. A message
summarises it.

## Details

A task is `"done"` when it ended well, `"failed"` when it ended with an
error, whose message is its reason, `"started"` when it began and has no
record of how it ended, `"queued"` when `submit.sh` submitted it after
its files were last written, so that it waits in the queue, and
`"not started"` otherwise. A task that started and has no record is
either still running or was stopped by Slurm before it could write one.
Slurm says which in the log it wrote for the task, in the folder
`logs/`, and the reason given is read from it: out of time, out of
memory, stopped by Slurm for a reason it gives (a node that failed,
say), or cancelled. Slurm's command `squeue -u $USER` (`sq` on the
clusters of the Alliance) lists the tasks still waiting or running.

When tasks remain to be run, the message gives the command that submits
them again, and only them: it leaves out a task that may still be
running or that waits in the queue, and a refit whose fit is itself
still queued or running, since the refit would start without its fit and
fail. A task already done is never run twice.

## See also

[`run_hpc_task()`](https://eogmartins.github.io/rTFT/reference/run_hpc_task.md),
which runs a task, and
[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md),
which assembles the results.

## Examples

``` r
if (FALSE) { # \dontrun{
# on the cluster, from the folder of the bundle; status.sh runs the same
hpc_status()
} # }
```
