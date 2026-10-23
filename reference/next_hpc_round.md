# Write the next round of a bundle, for some of its models

Writes a new bundle that fits again some of the models of a bundle
written by
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md)
(or of an earlier round), with new settings of the sampler and of the
jobs: typically the models that did not meet the sampling rule, given
more draws or a smaller step. The round is a bundle of its own, in a new
folder, copied to the cluster and run with the same scripts as the
first.

## Usage

``` r
next_hpc_round(
  dir,
  models,
  new_dir = NULL,
  ...,
  time = NULL,
  mem_per_cpu = NULL,
  email = NULL,
  remote = NULL,
  rtft_source = NULL,
  template = "alliance"
)
```

## Arguments

- dir:

  A character string, the folder of the bundle the round follows: the
  first bundle, or the latest round.

- models:

  A character vector, the names of the models to fit again, as in the
  column `model` of the first bundle's table of models.

- new_dir:

  An optional character string, the folder to write the round into. It
  is created if it does not exist, and must be empty if it does. `NULL`
  (the default) writes it beside `dir`, named after the first bundle and
  the number of the round: round 2 of `"my-run"` into
  `"my-run-round-2"`, and round 3, written from that, into
  `"my-run-round-3"`.

- ...:

  The settings of the sampler that change, among `chains`, `iter`,
  `warmup`, `thin`, `seed`, `control` and `init`.

- time, mem_per_cpu, email:

  Optional settings of the jobs that change, as in
  [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md).
  `NULL` (the default) keeps those of the bundle before.

- remote:

  An optional character string, where on the cluster to copy the round,
  as in
  [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md).
  `NULL` (the default) takes the bundle's, with its last folder renamed
  as the round's folder is (`scratch/my-run` becomes
  `scratch/my-run-round-2`).

- rtft_source:

  An optional character string, where the cluster installs rTFT from, as
  in
  [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md).
  `NULL` (the default) takes the bundle's.

- template:

  A character string, the scripts of the round, as in
  [`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md).

## Value

The path of the round, invisibly. A message says what was written and
what to do next.

## Details

Everything that makes the results of the rounds comparable is copied
from the bundle before, never chosen again: the data (the file itself),
the unit left out in cross-validation, the folds, the observation mode
and the sampling rule, as well as the account, the modules and the R
library of the cluster, where the bundle is copied to (with the round's
own folder name) and where rTFT is installed from. With cross-validation
by refitting, each model of the round is refitted on the same folds,
starting from its new fit.
[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md),
given the folders of every round, keeps for each model the result of its
latest round.

Every round holds the whole table of models of the first bundle, and
fits those named in `models`, so a model left out of one round can be
fitted again in a later one: a round can name any model of the first
bundle, whichever round it is written from.

Only the settings of the sampler and of the jobs can change, since
anything else would make another model: `chains`, `iter`, `warmup`,
`thin`, `seed`, `control` and `init` in `...`, and `time`, `mem_per_cpu`
and `email`. A setting not given is kept from the bundle before;
`control` is merged into the bundle's own, so that
`control = list(adapt_delta = 0.99)` keeps a `max_treedepth` given
before. To change a prior, write a new bundle with
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md).
The settings are one set for every model the round names, so a model
that needs more draws and one that needs a smaller step are two rounds.

The round must be written by the version of rTFT that wrote the first
bundle (and, between releases, from the same commit), so that every
round fits the models with the same code.

## See also

[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md),
which writes the first bundle, and
[`collect_hpc()`](https://eogmartins.github.io/rTFT/reference/collect_hpc.md),
which assembles the results of every round.

## Examples

``` r
if (FALSE) { # \dontrun{
# the model "exp_weibull" of the bundle "my-run" did not meet the convergence
# rule: fit it again with more draws and a smaller step, and longer jobs
# (written into the folder "my-run-round-2", beside "my-run")
next_hpc_round("my-run", "exp_weibull", iter = 4000,
               control = list(adapt_delta = 0.99), time = "24:00:00")

# once the results of both are copied back
out <- collect_hpc(c("my-run", "my-run-round-2"))
} # }
```
