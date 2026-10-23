# Fitting models on a cluster

Every model in the previous vignettes was fitted on one computer, and
most of them in minutes. Cross-validation by refitting is slower,
because each model is refitted once for each group of units left out,
and a comparison of several models multiplies the refits again: four
models compared by leaving out one of twelve trials at a time take 48
refits besides the four fits. A cluster (i.e., a set of computers shared
by many users, each with many cores) can be used to run these fits and
refits side by side, so that a comparison that would take a whole day on
a personal computer can be done in about the time of its slowest fit and
slowest refit (assuming the cluster has enough cores free).

In this vignette, we fit two models of a ramping assay on a cluster of
the [Digital Research Alliance of Canada](https://www.alliancecan.ca/en)
(hereafter, the Alliance), with cross-validation by refitting. We
describe: (i) the account and its setup, which are done once; (ii) a run
from start to end, from writing the models into a folder on our own
computer to assembling the results there; (iii) a second round for
models that did not meet the convergence rule; (iv) how long the tasks
took and how much memory they used; and (v) what differs on another
cluster. The models and the cross-validation are those of the vignettes
[*Ramping
assays*](https://eogmartins.github.io/rTFT/articles/ramping-assays.md)
and [*Comparing and averaging
models*](https://eogmartins.github.io/rTFT/articles/comparing-and-averaging-models.md).

None of the work on the cluster can run when this vignette is built.
Every block of terminal output below was captured from real runs on the
cluster Fir and is shown as it was printed, apart from the changes
listed at the end of the vignette.

``` r

library(rTFT)
```

## Clusters and bundles

A cluster has two kinds of computers. The login nodes are those we
connect to: they are shared by everyone connected, and are meant for
copying files, installing software and short checks. The compute nodes
do the work, and are reached only through the scheduler, a program
(Slurm, on the clusters of the Alliance) that holds the work submitted
by every user in a queue and starts each piece of work when the cores
and memory it asked for become free. A piece of work submitted to the
scheduler is a job. A job array is a set of similar jobs submitted at
once under one number, and each of its members, an array task, runs as a
job of its own, numbered by the array and its place in it (e.g., task 2
of the array 66774007 is 66774007_2); in `rTFT`, each fit and each refit
is one array task. Each job states how long it may run by the clock (its
wall time), how many cores it needs and how much memory per core; the
scheduler stops a job that exceeds its wall time or its memory.

`rTFT` includes functions that settle everything about a run on our own
computer, where the data and the models are, and write it into one
folder, the **bundle**: the data, the table of models, the settings of
the sampler and of the jobs, and the scripts that run it on the cluster.
The bundle is copied to the cluster and run there with its scripts, and
its results are copied back and assembled into what
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
returns, as if the models had been fitted on our own computer. Each fit
and each refit is one array task, so the work is spread over as many
cores as the cluster can give at the time.

## The account and its setup

The steps of this section are done once, and each is described on a page
of the Alliance’s documentation, which should be followed when it
differs from what is written here (the pages were checked on
2026-10-22). The clusters of the Alliance are open to researchers at
Canadian institutions and to the collaborators they sponsor. Elsewhere,
the same steps apply to most clusters that run Slurm, with the
differences listed in [Other clusters](#other-clusters); a cluster with
another scheduler needs scripts of its own.

### Account and allocation

An account is requested in the Alliance’s database of users, CCDB
([ccdb.alliancecan.ca](https://ccdb.alliancecan.ca)), with a sponsor:
the principal investigator of the research group, whose own account must
exist first ([Apply for a CCDB
account](https://docs.alliancecan.ca/wiki/Apply_for_a_CCDB_account)).
The time used on the clusters is charged to an allocation of the
sponsor’s group, whose name is given to every job. Every group has a
default allocation named after its principal investigator, such as
`def-someone`; on Fir it is given to jobs as `def-someone_cpu`, the form
the command `sshare -U` lists after connecting.

### Second factor of authentication

Every connection to a cluster of the Alliance asks for a second factor
of authentication (e.g., a notification on a telephone, approved in an
app) besides the password or key, which is set up in CCDB before the
first connection ([Multifactor
authentication](https://docs.alliancecan.ca/wiki/Multifactor_authentication)).
Each copy of files to or from a cluster is a connection of its own, and
asks for it too.

### SSH keys and connecting

A connection is made with `ssh`, from a terminal on our own computer, as
`ssh someone@fir.alliancecan.ca`, and asks for the password of the
account and then for the second factor. An SSH key (i.e., a pair of
files made on our own computer, of which the public half is added to
CCDB) can be used in place of the password, but not of the second factor
([SSH Keys](https://docs.alliancecan.ca/wiki/SSH_Keys)).

### File systems

Each user has three places to keep files on a cluster ([Storage and file
management](https://docs.alliancecan.ca/wiki/Storage_and_file_management)):
the home folder, small and backed up, for software and settings; the
project folder, shared by the group; and scratch (`~/scratch`), large
and not backed up, for the files of running work. Files in scratch that
have not been used for some time are deleted ([Scratch purging
policy](https://docs.alliancecan.ca/wiki/Scratch_purging_policy)), so a
bundle is run in scratch and its results are copied back to our own
computer as soon as they are done.

### Modules and the R library

Software on the clusters of the Alliance is loaded as modules:
[R](https://docs.alliancecan.ca/wiki/R), for instance, with
`module load StdEnv/2023 r/4.5.0`, where the first module sets the
generation of the software environment and the second the version of R.
Packages that R does not include are installed by each user into a
personal library in the home folder, and the samplers of `brms` need
CmdStan, which is compiled once into the home folder. Both are done by
the script `install.sh` of the bundle, the folder that `rTFT` writes for
a run on a cluster (see [The bundle](#the-bundle)).

## A run from start to end

The example is the ramping assay of the vignette [*Ramping
assays*](https://eogmartins.github.io/rTFT/articles/ramping-assays.md),
restricted to the ramps at 0.1 and 0.3°C per min in normoxia: 8 trials
of fish whose failures were each observed at a known time. We compare
two models, both with the exponential failure rate and a random effect
of trial on the failure rate at the reference temperature, one with the
memoryless baseline and one with the Weibull baseline, by refitting each
model with one trial left out at a time, the trials grouped into 4
folds.

### The bundle

The data and the models are built as in the previous vignettes.

``` r

rates <- c(0.1, 0.3)
ramp  <- ramp_assay[ramp_assay$rate %in% rates & ramp_assay$do == 10, ]
path  <- ramp_path(path_id = paste0("r", rates), T_start = 14, rate = rates,
                   t_max = 15.5 / rates)
long  <- build_hazard_data(ramp, path, observation = "continuous",
                           covariates = c("mass", "trial"))
models <- build_hazard_models("exp", c("mem", "weibull"),
                              dpar_formulas = list(loga = ~ 1 + (1 | trial)))
```

The function
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
writes the bundle. It takes what
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
takes, and the allocation, the wall time and memory of each task, and
where on the cluster the bundle will be copied. Each task runs the four
chains of one fit on four cores, one chain per core. We gave each task a
wall time of one minute, far less than a fit needs, to show what happens
when a task runs out of time; a first run should be given ample time
instead (see [Time and memory](#time-and-memory)). The settings of the
sampler (1,000 iterations per chain, at the default step) are kept low
on purpose too, so that a second round is needed and can be shown below;
a real run starts from the settings of the other vignettes (2,000
iterations and `adapt_delta = 0.95`).

``` r

write_hpc_bundle(models, long, "ramp-run", account = "def-someone_cpu",
                 observation = "continuous", kfold = TRUE, K = 4,
                 time = "00:01:00",
                 remote = "someone@fir.alliancecan.ca:scratch/ramp-run",
                 chains = 4, iter = 1000, seed = 1, backend = "cmdstanr")
```

    #> ℹ Leaving out one trial at a time: the individuals of each trial share its
    #>   effect in the model, so leaving out one individual would leave that effect in
    #>   the fit.
    #> ℹ Use `unit = "id"` to leave out individuals instead.
    #> ✔ Wrote a bundle of 2 models to 'ramp-run': 2 fits and 8 refits for
    #>   cross-validation, leaving out one trial at a time.
    #> ℹ Copy it to the cluster with the command in step 1 of 'README.txt', in that
    #>   folder.
    #> ℹ Then, on the cluster, in that folder, run these two (the first once on each
    #>   cluster):
    #>   bash install.sh
    #>   bash check.sh

The bundle holds the scripts, the file `README.txt` and four folders.
The folder `inputs/` holds what the run is made of (the data, the table
of models and the settings, with their checksums), and nothing in it is
to be edited; the other three are filled by the tasks: `stan/` with the
compiled models, `results/` with what each task finds, and `logs/` with
what each task prints.

``` r

list.files("ramp-run")
#>  [1] "check.sh"    "inputs"      "install.R"   "install.sh"  "logs"       
#>  [6] "README.txt"  "results"     "settings.sh" "stan"        "status.sh"  
#> [11] "submit.sh"   "task.sh"
```

The file `settings.sh` keeps every setting of the jobs in one place, and
is the one file of the bundle that may be edited. The file `README.txt`
says what to do, step by step, with the commands to copy the bundle
filled in.

``` r

cat(readLines("ramp-run/README.txt"), sep = "\n")
#> A bundle of models to fit on a cluster
#> ======================================
#> 
#> Written by rTFT::write_hpc_bundle() on 2026-10-23 20:54 UTC, with rTFT 0.2.0.
#> It holds 2 models to fit, as 2 tasks, and 8 refits for cross-validation by
#> refitting.
#> 
#> What is in it
#> 
#>   settings.sh   the account, time, memory, cores and modules, in one place
#>   install.sh    installs rTFT, cmdstanr and CmdStan into the cluster's R
#>   install.R     the R part of install.sh
#>   check.sh      checks the bundle and the cluster's R, on the login node
#>   submit.sh     submits the tasks to Slurm, the cluster's scheduler
#>   task.sh       what Slurm runs for each task; not run by hand
#>   status.sh     says what each task has done, and what to do next
#>   inputs/       the data, the models and the settings of the run, with their
#>                 checksums; nothing in it is to be edited
#>   stan/         where the compiled Stan programs will be kept
#>   results/      where the tasks will write what they find
#>   logs/         where Slurm will write what each task printed
#>   tmp/          where install.sh keeps its temporary files on a login node
#> 
#> What to do
#> 
#>   1. Copy the folder to the cluster, from your own computer:
#> 
#>        rsync -av --rsync-path="mkdir -p scratch/ramp-run && rsync" \
#>          "/tmp/RtmpoYgL0o/ramp-run/" \
#>          someone@fir.alliancecan.ca:scratch/ramp-run/
#> 
#>      It is one command over three lines; copy all three. The --rsync-path part
#>      makes the folder on the cluster first, with any folders above it that do
#>      not exist yet; rsync alone makes only the last.
#> 
#>      Many clusters ask for a second factor of authentication besides your
#>      password, at each copy and each connection; those of the Digital Research
#>      Alliance of Canada do. If yours does, it must be set up first, as the
#>      cluster's documentation explains.
#> 
#>   2. Connect to the cluster and go to the folder:
#> 
#>        ssh someone@fir.alliancecan.ca
#>        cd scratch/ramp-run
#> 
#>   3. Install what the bundle needs, once on each cluster (and again for a
#>      bundle written by another version of rTFT):
#> 
#>        bash install.sh
#> 
#>      The first time, it compiles brms and CmdStan, which can take an hour or
#>      more. It prints a great deal; its last lines say whether it worked. If
#>      the login node stops it, run it as a job instead, and wait for the job
#>      to finish:
#> 
#>        bash install.sh --job
#> 
#>      The job downloads what it installs. On a cluster whose compute nodes
#>      cannot reach the internet it stops at once and says so; ask the
#>      cluster's support how R packages are installed there.
#> 
#>   4. Check the bundle and the cluster's R:
#> 
#>        bash check.sh
#> 
#>      check.sh reports anything missing and what to do about it. Run it again
#>      until it reports nothing missing.
#> 
#>   5. Submit the tasks:
#> 
#>        bash submit.sh
#> 
#>      Each fit is one task of a job array, a set of jobs that Slurm runs side
#>      by side as cores become free, and each refit for cross-validation is one
#>      task of a second array, which waits for the fits to end, since each refit
#>      starts from its model's fit. submit.sh prints the number Slurm gave each
#>      array. A task can wait in the queue for minutes or for hours, depending
#>      on how busy the cluster is. The tasks run whether or not you stay
#>      connected: you can log out, and run status.sh when you come back.
#> 
#>   6. Follow the tasks:
#> 
#>        sq                  the tasks waiting (PD) and running (R)
#>        bash status.sh      what each task has done, and what to do next
#> 
#>      What each task prints is in logs/, one file per task.
#> 
#>   7. When a task failed, ran out of time or memory, or was stopped by Slurm,
#>      status.sh says so and gives the command that submits again only the
#>      tasks not done. If a task ran out of time or memory, raise TIME or
#>      MEM_PER_CPU in settings.sh first. The command leaves out tasks still
#>      waiting or running; do not submit those again by hand, or they will run
#>      twice.
#> 
#>   8. When status.sh says every task is done, copy the results back, in a
#>      terminal on your own computer, not the one connected to the cluster:
#> 
#>        rsync -av \
#>          someone@fir.alliancecan.ca:scratch/ramp-run/results/ \
#>          "/tmp/RtmpoYgL0o/ramp-run/results/"
#> 
#>      results/ holds every fit, and can take some gigabytes. If you have moved
#>      the bundle's folder on your computer since writing it, change the path
#>      in the command, and in step 9.
#> 
#>   9. Assemble them in R, on your own computer:
#> 
#>        out <- rTFT::collect_hpc("/tmp/RtmpoYgL0o/ramp-run")
#>        out$summary
#> 
#>      You can collect before every task is done: models without a result yet
#>      are named, and step 8 can be run again to bring back the rest.
```

### Copying the bundle

Step 1 of `README.txt` gives the command that copies the bundle to the
cluster, run from a terminal on our own computer, in the form that suits
the computer that wrote the bundle. On macOS and Linux it is one `rsync`
command over three lines, whose second line names the folder where the
bundle was written (here a temporary folder, since this vignette writes
it there). Windows has no `rsync`, and there the command is two: `ssh`,
which makes the folders above the bundle’s on the cluster, and `scp`,
which copies the bundle; graphical tools that the Alliance documents
(e.g., Globus, WinSCP) do the same.

### Installing rTFT, brms and CmdStan

After connecting and going to the bundle’s folder, `install.sh` installs
what the bundle needs into R on the cluster: `rTFT` at the version that
wrote the bundle, with `brms` and the other packages it needs,
`cmdstanr`, and CmdStan. A first install compiles a great deal, and a
login node may stop a long compilation, so it can be run as a job
instead.

``` bash
ssh someone@fir.alliancecan.ca
cd scratch/ramp-run
bash install.sh --job
```

``` text
Submitted batch job 66770460
Submitted. When the job has finished, read logs/install-<job number>.out,
then run:  bash check.sh
```

The job’s log ends by saying whether the install worked. A first install
into an empty library installed 77 packages in 1 h 48 min (the 72 that
the log lists in one line below, with `remotes`, `rTFT`, `cmdstanr` and
the 2 packages it needs), so a first install should be expected to take
an hour or more; selected lines of its log are below.

``` text
Installing into /home/someone/scratch/rtft-lib-empty
rTFT: installing version 0.1.1.9000 from eogmartins/rTFT@3016af22758c0faf57d62cae07e19e80072c9971 (installed: none)
Installing 72 packages: rlang, lifecycle, glue, cli, vctrs, utf8, pkgconfig, magrittr, numDeriv, generics, backports, matrixStats, distributional, pillar, tensorA, tibble, checkmate, abind, viridisLite, RColorBrewer, R6, labeling, farver, cpp11, withr, parallelly, listenv, digest, globals, future, stringi, scales, coda, stringr, Brobdingnag, mvtnorm, purrr, plyr, S7, isoband, Rcpp, quadprog, gtable, tidyselect, tidyr, reshape2, posterior, ggridges, ggplot2, ggdist, dplyr, QuickJSR, RcppParallel, desc, ps, otel, processx, callr, BH, RcppEigen, pkgbuild, loo, gridExtra, inline, StanHeaders, nleqslv, future.apply, bridgesampling, bayesplot, rstantools, rstan, brms
* DONE (brms)
* DONE (rTFT)
cmdstanr: installing
* DONE (cmdstanr)
CmdStan is installed at /home/someone/.cmdstan/cmdstan-2.40.0
Done. Now run:  bash check.sh
```

CmdStan had been installed by an earlier install, so this one found it.
When it is missing, `install.sh` compiles it too, which ends with these
lines (from an install on the login node of Fir on 2026-10-21):

``` text
--- CmdStan v2.40.0 built ---
* Finished installing CmdStan to /home/someone/.cmdstan/cmdstan-2.40.0
CmdStan path set to: /home/someone/.cmdstan/cmdstan-2.40.0
```

The install is needed once on each cluster, and again for a bundle
written by another version of `rTFT`, when it installs only `rTFT`.

### Checking the bundle

`check.sh` checks that the files of the bundle are those that were
written, that the installed `rTFT` is the version that wrote the bundle,
that `cmdstanr` finds CmdStan, and that the models can be assembled. Its
first lines, printed by the cluster each time R is loaded, say which
modules were loaded in place of those loaded at connection; they are not
shown in the blocks that follow.

``` bash
bash check.sh
```

``` text
The following have been reloaded with a version change:
  1) StdEnv/2020 => StdEnv/2023
  2) flexiblas/3.0.4 => flexiblas/3.3.1
  3) gcc/9.3.0 => gcc/12.3
  4) gcccore/.9.3.0 => gcccore/.12.3
  5) gentoo/2020 => gentoo/2023
  6) libfabric/1.10.1 => libfabric/1.18.0
  7) openmpi/4.0.3 => openmpi/4.1.5
  8) ucx/1.8.0 => ucx/1.14.1

✔ Every file of the bundle is as it was written (5 files).
✔ rTFT 0.1.1.9000, the version that wrote the bundle; brms 2.23.0; CmdStan
  found.
✔ 2 models assembled, with 8 refits for cross-validation.
ℹ To submit the tasks:
  bash submit.sh
```

### Submitting the tasks

`submit.sh` submits the fits as one job array and the refits as a
second, which waits for the fits to end, since each refit starts from
its model’s fit.

``` bash
bash submit.sh
```

``` text
Submitted the fits 1-2 of 2 as job 66774007.
Submitted the refits 1-8 of 8 as job 66774008,
to start when the fits have ended.

Follow them with:  sq               the tasks waiting (PD) and running (R)
                   bash status.sh   what each task has done, and what to do next
```

### Following the tasks

The command `sq` lists our jobs waiting (PD, pending) and running (R),
with the reason a pending job waits: `(None)` before the scheduler has
considered it, `(Priority)` behind the jobs of other users, and
`(Dependency)` for the refits waiting for the fits. We listed only the
jobs of this run, by their names, since `sq` lists every job of the
user.

``` bash
sq -n rtft-fit,rtft-kfold
```

``` text
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
 66774007_[1-2] someone def-someone       rtft-fit  PD       1:00     1    4        N/A   4000M  (None)
 66774008_[1-8] someone def-someone     rtft-kfold  PD       1:00     1    4        N/A   4000M  (None)
```

`status.sh` says what each task has done, from what the tasks wrote into
the bundle, and what to do next. It can be run at any time, as it only
reports the status of the tasks and changes nothing.

### A task that ran out of time

With a wall time of one minute, both fits were stopped by the scheduler
after 1 min 28 s. `status.sh` read the reason from the log of each task,
and gave the command to submit the fits again. The refits were still
waiting.

``` bash
bash status.sh
```

``` text
• Fits: 0 of 2 done; 2 started.
• Refits for cross-validation: 0 of 8 done; 8 queued.
! fit 1 (exp_mem): started, it ran out of time (raise TIME in settings.sh); see logs/fit-66774007_1.out
! fit 2 (exp_weibull): started, it ran out of time (raise TIME in settings.sh); see logs/fit-66774007_2.out
• refit 1 (exp_mem, fold 1): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 2 (exp_mem, fold 2): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 3 (exp_mem, fold 3): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 4 (exp_mem, fold 4): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 5 (exp_weibull, fold 1): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 6 (exp_weibull, fold 2): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 7 (exp_weibull, fold 3): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
• refit 8 (exp_weibull, fold 4): queued, submitted at 2026-10-22 07:07 as job 66774008, waiting to start
ℹ To submit again the tasks not done, and only those:
  bash submit.sh --fit 1-2
ℹ A task submitted and waiting to start is not submitted again. One that the command sq no longer lists was cancelled before it started; to submit those:
  bash submit.sh --kfold 1-8
```

`sq` then listed the refits alone, waiting on the fits that had ended:

``` text
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
 66774008_[1-8] someone def-someone     rtft-kfold  PD       1:00     1    4        N/A   4000M  (Dependency)
```

The refits then started, since they wait for the fits to end whether or
not they ended well, and each stopped within 3 s, having no fit to start
from.

``` text
• Fits: 0 of 2 done; 2 started.
• Refits for cross-validation: 0 of 8 done; 8 failed.
! fit 1 (exp_mem): started, it ran out of time (raise TIME in settings.sh); see logs/fit-66774007_1.out
! fit 2 (exp_weibull): started, it ran out of time (raise TIME in settings.sh); see logs/fit-66774007_2.out
✖ refit 1 (exp_mem, fold 1): failed. The refits of "exp_mem" start from its fit, task 1 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 1 --kfold 1 submits the two, the refit to start when the fit has ended.
✖ refit 2 (exp_mem, fold 2): failed. The refits of "exp_mem" start from its fit, task 1 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 1 --kfold 2 submits the two, the refit to start when the fit has ended.
✖ refit 3 (exp_mem, fold 3): failed. The refits of "exp_mem" start from its fit, task 1 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 1 --kfold 3 submits the two, the refit to start when the fit has ended.
✖ refit 4 (exp_mem, fold 4): failed. The refits of "exp_mem" start from its fit, task 1 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 1 --kfold 4 submits the two, the refit to start when the fit has ended.
✖ refit 5 (exp_weibull, fold 1): failed. The refits of "exp_weibull" start from its fit, task 2 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 2 --kfold 5 submits the two, the refit to start when the fit has ended.
✖ refit 6 (exp_weibull, fold 2): failed. The refits of "exp_weibull" start from its fit, task 2 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 2 --kfold 6 submits the two, the refit to start when the fit has ended.
✖ refit 7 (exp_weibull, fold 3): failed. The refits of "exp_weibull" start from its fit, task 2 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 2 --kfold 7 submits the two, the refit to start when the fit has ended.
✖ refit 8 (exp_weibull, fold 4): failed. The refits of "exp_weibull" start from its fit, task 2 of the fits, which has not ended well. Run that task first; bash submit.sh --fit 2 --kfold 8 submits the two, the refit to start when the fit has ended.
ℹ To submit again the tasks not done, and only those:
  bash submit.sh --fit 1-2 --kfold 1-8
```

The remedy is to raise the wall time in `settings.sh` (here, with a text
editor directly in the terminal, such as `nano settings.sh`, from
`TIME="00:01:00"` to `TIME="01:00:00"`) and to submit again the tasks
that `status.sh` names. A task that has ended well is never run again,
so the same command can be used after any number of failures.

``` bash
bash submit.sh --fit 1-2 --kfold 1-8
```

``` text
Submitted the fits 1-2 of 2 as job 66774407.
Submitted the refits 1-8 of 8 as job 66774408,
to start when the fits have ended.

Follow them with:  sq               the tasks waiting (PD) and running (R)
                   bash status.sh   what each task has done, and what to do next
```

While the fits ran, `sq` and `status.sh` said so.

``` text
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
     66774407_1 someone def-someone       rtft-fit   R      59:29     1    4        N/A   4000M fc30565 (None)
     66774407_2 someone def-someone       rtft-fit   R      59:29     1    4        N/A   4000M fc30622 (None)
 66774408_[1-8] someone def-someone     rtft-kfold  PD    1:00:00     1    4        N/A   4000M  (Dependency)
```

``` text
• Fits: 0 of 2 done; 2 started.
• Refits for cross-validation: 0 of 8 done; 8 queued.
! fit 1 (exp_mem): started, no result yet: it may still be running, or it stopped without one; see logs/fit-66774407_1.out
! fit 2 (exp_weibull): started, no result yet: it may still be running, or it stopped without one; see logs/fit-66774407_2.out
• refit 1 (exp_mem, fold 1): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 2 (exp_mem, fold 2): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 3 (exp_mem, fold 3): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 4 (exp_mem, fold 4): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 5 (exp_weibull, fold 1): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 6 (exp_weibull, fold 2): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 7 (exp_weibull, fold 3): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
• refit 8 (exp_weibull, fold 4): queued, submitted at 2026-10-22 07:10 as job 66774408, waiting to start
ℹ A task started without a record of its end may still be running; the command sq lists the tasks waiting or running.
ℹ A task submitted and waiting to start is not submitted again. One that the command sq no longer lists was cancelled before it started; to submit those:
  bash submit.sh --kfold 1-8
```

When every task had ended well:

``` text
• Fits: 2 of 2 done.
• Refits for cross-validation: 8 of 8 done.
✔ Every task is done. Copy 'results/' back to your own computer as step 8 of 'README.txt' says, and assemble the results there with `collect_hpc()`.
```

### Copying the results back

The command of step 8 of `README.txt`, run in a terminal on our own
computer, copies the folder `results/` back into the bundle there. It
holds each fit and, for each refit, the log predictive density of each
trial of its fold, which take some megabytes for these models and can
take some gigabytes for larger ones. The results should be copied back
as soon as the tasks are done, since scratch is purged.

### Assembling the results

[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md)
reads the results into what
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
returns: the fitted models, their summary, and the cross-validation by
refitting, which the functions of the vignette [*Comparing and averaging
models*](https://eogmartins.github.io/rTFT/articles/comparing-and-averaging-models.md)
take as they are.

``` r

fits <- collect_hpc("ramp-run")
fits$summary[, c("model", "rule_met", "rhat", "ess_bulk", "ess_tail", "divergent",
                 "elpd_kfold", "se_kfold")]
```

``` text
✔ Collected 2 of the 2 models, 2 cross-validated by refitting.
Warning messages:
1: ! 2 models did not meet the sampling rule: "exp_mem" and "exp_weibull".
ℹ The rule is Rhat below 1.01, bulk and tail effective sample sizes above 400,
  and at most 0.1% of the draws from divergent transitions.
ℹ They are kept in the fits returned; see rhat, ess_bulk, ess_tail and
  divergent in the summary.
2: ! PSIS cannot be trusted for some of the units left out, in 2 models: "exp_mem"
  and "exp_weibull".
ℹ The number of such units is in high_k of the summary.
ℹ The estimates by refitting, in kfold, do not rely on it.
# A tibble: 2 × 8
  model       rule_met  rhat ess_bulk ess_tail divergent elpd_kfold se_kfold
  <chr>       <lgl>    <dbl>    <dbl>    <dbl>     <int>      <dbl>    <dbl>
1 exp_mem     FALSE     1.01     283.     372.         3      -590.     18.4
2 exp_weibull FALSE     1.01     398.     363.         0      -588.     18.7
```

Neither model met the sampling rule (`rule_met` in the summary: Rhat
below 1.01, bulk and tail effective sample sizes above 400 and at most
0.1% of the draws from divergent transitions; see
[`?fit_hazard_models`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)).
Their Rhat, printed as 1.01, was 1.015 and 1.014, above the limit of
1.01, and their bulk effective sample sizes fell short of 400.

## Rounds for models that did not converge

A model that did not meet the sampling rule is fitted again with other
settings of the sampler: more draws when the effective sample sizes are
low, and a smaller step (a higher `adapt_delta`) when divergent
transitions remain. The function
[`next_hpc_round()`](https://eogmartins.github.io/rTFT/reference/next_hpc_round.md)
writes a new bundle, a round, for the models named, beside the bundle it
follows, with the same data, unit and folds, so that the results of
every round can be assembled together. Only the settings of the sampler
and of the jobs can change in a round, since anything else would make
another model.

Both models of the run above missed the rule. We first fitted the
Weibull model again with three times as many iterations.

``` r

next_hpc_round("ramp-run", "exp_weibull", iter = 3000)
```

``` text
✔ Wrote round 2 of the bundle
  '~/my-runs/ramp-run' to
  '~/my-runs/ramp-run-round-2': 1 fit and
  4 refits for cross-validation, of "exp_weibull".
ℹ Copy it to the cluster with the command in step 1 of its 'README.txt', and
  run it as you ran the first.
ℹ Then assemble the results of every round with `collect_hpc(c(...))`, giving
  the folders of the rounds.
```

The round is written beside the bundle it follows, into a folder named
after the first bundle and the number of the round: `ramp-run-round-2`
here, and `ramp-run-round-3` for a round written from it. It has its own
`README.txt`, whose commands copy the round to the cluster and its
results back, and it is run there with the same scripts as the first
bundle (`check.sh`, `submit.sh` and `status.sh`), without a new install,
since it was written by the same version of `rTFT`. The round is copied
to a folder of its own on the cluster, beside the first, and its scripts
are run from that folder, which step 2 of its `README.txt` goes to
(here, `cd scratch/ramp-run-round-2`).

``` text
✔ Every file of the bundle is as it was written (5 files).
✔ rTFT 0.1.1.9000, the version that wrote the bundle; brms 2.23.0; CmdStan
  found.
✔ 1 model assembled, with 4 refits for cross-validation.
ℹ To submit the tasks:
  bash submit.sh
```

When its tasks were done, its results were copied back with the command
of step 8 of its `README.txt`. Its bulk effective sample size rose from
398 to 1,298, but 15 of its 6,000 draws came from divergent transitions
(0.25%, against the limit of 0.1%). A third round, written from the
second, fitted both models with a smaller step; the number of iterations
was kept from the round before.

``` r

next_hpc_round("ramp-run-round-2", c("exp_mem", "exp_weibull"),
               control = list(adapt_delta = 0.99))
```

``` text
✔ Wrote round 3 of the bundle
  '~/my-runs/ramp-run-round-2' to
  '~/my-runs/ramp-run-round-3': 2 fits
  and 8 refits for cross-validation, of "exp_mem" and "exp_weibull".
ℹ Copy it to the cluster with the command in step 1 of its 'README.txt', and
  run it as you ran the first.
ℹ Then assemble the results of every round with `collect_hpc(c(...))`, giving
  the folders of the rounds.
```

The third round was copied to the cluster, checked, submitted and copied
back in the same way, with the commands of its own `README.txt`.
[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md),
given the folders of every round, keeps for each model the result of its
latest round.

``` r

fits <- collect_hpc(c("ramp-run", "ramp-run-round-2", "ramp-run-round-3"))
fits$summary[, c("model", "round", "rule_met", "rhat", "ess_bulk", "ess_tail",
                 "divergent", "draws", "elpd_kfold")]
```

``` text
✔ Collected 2 of the 2 models, 2 cross-validated by refitting, from 3 rounds.
Warning message:
! PSIS cannot be trusted for some of the units left out, in 2 models: "exp_mem"
  and "exp_weibull".
ℹ The number of such units is in high_k of the summary.
ℹ The estimates by refitting, in kfold, do not rely on it.
# A tibble: 2 × 9
  model       round rule_met  rhat ess_bulk ess_tail divergent draws elpd_kfold
  <chr>       <int> <lgl>    <dbl>    <dbl>    <dbl>     <int> <dbl>      <dbl>
1 exp_mem         3 TRUE      1.01    1059.    1815.         1  6000      -589.
2 exp_weibull     3 TRUE      1.00    1306.    1902.         0  6000      -588.
```

Both models met the rule in the third round. The Rhat of the memoryless
model, printed as 1.01, was 1.006, below the limit, and 1 of its 6,000
draws came from a divergent transition. The assembled models are
compared and averaged as in the vignette [*Comparing and averaging
models*](https://eogmartins.github.io/rTFT/articles/comparing-and-averaging-models.md):

``` r

hazard_weights(fits$kfold)
```

``` text
ℹ Weights by stacking, from 8 units left out (one trial each).
# A tibble: 2 × 6
  model        elpd weight_all cumulative kept     weight
  <chr>       <dbl>      <dbl>      <dbl> <lgl>     <dbl>
1 exp_weibull -588.  1.000          1.000 TRUE  1.000
2 exp_mem     -589.  0.0000314      1     TRUE  0.0000314
```

The weights should be read beside the standard errors. The two models
predict the trials left out equally well (their elpd by refitting
differs by 1.0, with standard errors of about 19), and the estimate by
refitting itself moves by about 1 between runs, because the effect of a
trial left out is drawn at random (see
[`?kfold_hazard`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)).
With eight trials, stacking puts nearly all of its weight on one of two
near-equivalent models, where pseudo-BMA (`method = "pseudobma"`) would
split it; the vignette [*Comparing and averaging
models*](https://eogmartins.github.io/rTFT/articles/comparing-and-averaging-models.md)
discusses both.

## Time and memory

In the run above, after the wall time was raised, the two fits took 3
min 48 s and 4 min 17 s and the eight refits 18 s to 36 s each, as the
tasks measured them, from the start of R to the end of the
cross-validation by importance sampling; the time of a fit includes
compiling its model, which its refits reuse. The scheduler’s clock, by
which the allocation is charged, adds the loading of R and the modules
to each task, as `sacct` listed for the jobs of this run (the job, its
name, how it ended and how long it ran):

``` text
JobID           JobName      State    Elapsed
------------ ---------- ---------- ----------
66774007_1     rtft-fit    TIMEOUT   00:01:28
66774007_2     rtft-fit    TIMEOUT   00:01:28
66774008_1   rtft-kfold     FAILED   00:00:03
66774008_2   rtft-kfold     FAILED   00:00:03
66774008_3   rtft-kfold     FAILED   00:00:02
66774008_4   rtft-kfold     FAILED   00:00:02
66774008_5   rtft-kfold     FAILED   00:00:03
66774008_6   rtft-kfold     FAILED   00:00:03
66774008_7   rtft-kfold     FAILED   00:00:03
66774008_8   rtft-kfold     FAILED   00:00:03
66774407_1     rtft-fit  COMPLETED   00:03:50
66774407_2     rtft-fit  COMPLETED   00:04:20
66774408_1   rtft-kfold  COMPLETED   00:00:25
66774408_2   rtft-kfold  COMPLETED   00:00:24
66774408_3   rtft-kfold  COMPLETED   00:00:25
66774408_4   rtft-kfold  COMPLETED   00:00:24
66774408_5   rtft-kfold  COMPLETED   00:00:38
66774408_6   rtft-kfold  COMPLETED   00:00:42
66774408_7   rtft-kfold  COMPLETED   00:00:39
66774408_8   rtft-kfold  COMPLETED   00:00:38
```

The tasks that completed cost 0.83 core-hours with 4 cores per task, and
the fits that ran out of time and the refits that failed another 0.22. A
larger run can be judged from a small one: the time of a fit grows
roughly in proportion to the number of rows of the table of intervals at
risk and to the number of iterations, and a model that converges poorly
(with many steps per iteration) can take many times longer than its
neighbour. A wall time of two to three times the longest fit of a small
run leaves a margin; a task that asks for more waits longer in the
queue, but a task that runs out of time is lost, and is run again from
its start.

The fits used 3.0 and 3.5 GB and the refits about 1.0 GB, against the
3.9 GB per core asked for (`MEM_PER_CPU` in `settings.sh` is 4000M, in
megabytes), about 16 GB per task. The command `seff` with the number of
a finished job reports its time and memory, as the scheduler measured
them. A task that runs out of memory is reported by `status.sh`, which
gives the command to submit it again after `MEM_PER_CPU` is raised in
`settings.sh`.

## Other clusters

The bundle’s scripts are written for the clusters of the Alliance, and
most of what they assume holds on any cluster that runs Slurm. The table
lists what we checked on Fir and Nibi, two clusters of the Alliance, and
what to look for elsewhere.

|  | Fir | Nibi | Another cluster |
|----|----|----|----|
| Account given to jobs | `def-someone_cpu` | `def-someone_cpu` | the cluster’s documentation |
| Second factor at each connection | yes | yes | the cluster’s documentation |
| Modules for R | `StdEnv/2023 r/4.5.0` | `StdEnv/2023 r/4.5.0` | `module avail r` |
| Compute nodes reach the internet (`install.sh --job`) | yes | yes | the job stops at once and says so when they do not |
| Command listing one’s jobs | `sq` | `sq` | `squeue -u $USER` |
| Largest job array | 10,000 | 10,000 | `scontrol show config` (`MaxArraySize`) |
| Longest wall time | 7 days | 7 days | the cluster’s documentation |
| Scratch | `~/scratch` | `~/scratch` | the cluster’s documentation |
| Checked | 2026-10-22 | 2026-10-19 |  |

The scripts themselves can be replaced:
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
takes, in `template`, a folder of scripts of one’s own with the same
names, written for another scheduler or another cluster, and fills them
in as it fills its own (see its help page).

> **Captured output**
>
> The output of this vignette was captured on the cluster Fir and on the
> computer that wrote the bundles: the run from start to end on
> 2026-10-22 (`rTFT` at commit 0c30a81 of the development version that
> became 0.2.0), the first install on 2026-10-22 (at 3016af2), the
> rounds on 2026-10-21 (at 96187ae), from an earlier run of the same
> bundle on the same cluster, which gave the same fits (the sampler is
> seeded; a seed does not reproduce draws across machines) and estimates
> by refitting that differ by under 1 in elpd, since the effects of the
> trials left out were drawn with an unseeded generator in that run; and
> the build of CmdStan on 2026-10-21, from a second install on the login
> node; all with `brms` 2.23.0. Six changes were made to it: (i) the
> user name is replaced by `someone`, (ii) the folder of the runs on
> that computer by `~/my-runs`, and (iii) the name of every bundle by
> `ramp-run`; (iv) the notice of the modules reloaded is shown once; (v)
> the logs of the installs are shown by selected lines; and (vi) the
> column `converged` of the summaries and the warning’s “convergence
> rule” have the names given to them afterwards, `rule_met` and
> “sampling rule”. `sq` cuts the account to 12 characters, so
> `def-someone` in its output stands for `def-someone_cpu`. The bundle
> and its `README.txt` are written when the vignette is built.
