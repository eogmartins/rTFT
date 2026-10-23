# Check a bundle and the R that will fit it, before anything is submitted

Checks a bundle written by
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md),
where it is to be run: on the login node of the cluster, from the
bundle's script `check.sh`, which loads R as the bundle's settings say
and calls this function. It takes seconds, fits nothing and installs
nothing.

## Usage

``` r
check_hpc_bundle(dir = ".")
```

## Arguments

- dir:

  A character string, the folder of the bundle. Defaults to the working
  folder, which is where `check.sh` runs it.

## Value

`TRUE`, invisibly, when every check passes. A message lists what was
checked.

## Details

It checks, in order, that the files of the bundle are all there and
unchanged since they were written (by their checksums, so a copy that
was cut short or a file edited by hand is found); that the installed
rTFT is the version that wrote the bundle; that cmdstanr is installed
and finds CmdStan; and that every model of the bundle can be assembled
on its data. The first problem stops the check, with an error that says
what to do about it.

The version is compared by its number and, when both the rTFT that wrote
the bundle and the one installed came from GitHub, by the commit each
was installed from, since copies of a development version share a number
and can differ in their code. A copy installed from a folder on your own
machine has no commit to compare, so for it only the number is checked.
[`run_hpc_task()`](https://eogmartins.github.io/rTFT/reference/run_hpc_task.md)
makes the same check at the start of each task, since the library can
change while tasks wait in the queue.

## See also

[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md),
which writes the bundle.

## Examples

``` r
if (FALSE) { # \dontrun{
# on the cluster, from the folder of the bundle; check.sh runs the same
check_hpc_bundle()
} # }
```
