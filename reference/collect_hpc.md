# Assemble the results of a bundle run on a cluster

Assembles the results of the tasks of a bundle written by
[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md),
once its folder `results/` has been copied back from the cluster, into
what
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
returns: the fits, their cross-validation results and a summary of one
row per model, ready for
[`hazard_weights()`](https://eogmartins.github.io/rTFT/reference/hazard_weights.md)
and the functions that average models.

## Usage

``` r
collect_hpc(dir = ".")
```

## Arguments

- dir:

  A character vector, the folder of the bundle, on your own computer,
  with the folder `results/` copied back into it, or the folders of the
  bundle and of its rounds. Defaults to the working folder.

## Value

A list with the four elements that
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
returns: `fits`, `loo`, `kfold` and `summary`. In `summary`, `round` is
the round each model's result comes from (1 for the first bundle),
`status` is `"fitted"`, `"failed"` or `"missing"`, and `seconds` is the
time each fit and its cross-validation took on the cluster.

## Details

Given the folders of a bundle and of its rounds, written by
[`next_hpc_round()`](https://eogmartins.github.io/rTFT/reference/next_hpc_round.md),
it assembles each model from its latest round, the last round that
fitted it, whether the result of that round is there, failed or missing,
so that a model fitted again is never given its earlier fit without a
word. The folders must be the first bundle and every round from the
second to the last, all holding the same data, and can be given in any
order.

Nothing is assembled from a bundle whose files have changed since it was
written, or with a result computed from another bundle: the files are
checked against their checksums, and each result against the stamp of
the bundle that computed it.

The sampling rule of the bundle is applied to each fit and reported as
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
reports it. A model whose fit failed is reported with its reason, and
one whose fit has no result yet (its task did not run, or has not
finished) is marked `"missing"` in the summary and named in a warning,
so the results can be looked at before every task is done. A model is
cross-validated by refitting only when every one of its folds has a
result with a log predictive density defined for each of its units;
otherwise the reason is given in the summary. A density of zero, a log
predictive density of `-Inf`, is kept, as
[`kfold_hazard()`](https://eogmartins.github.io/rTFT/reference/kfold_hazard.md)
keeps it.

## See also

[`write_hpc_bundle()`](https://eogmartins.github.io/rTFT/reference/write_hpc_bundle.md),
which writes the bundle, and
[`hpc_status()`](https://eogmartins.github.io/rTFT/reference/hpc_status.md),
which says what each task has done.

## Examples

``` r
if (FALSE) { # \dontrun{
# on your own computer, once results/ has been copied back into the bundle
out <- collect_hpc("my-run")
out$summary
weights <- hazard_weights(out$kfold)

# with the round written by next_hpc_round() for a model that did not converge
out <- collect_hpc(c("my-run", "my-run-round-2"))
} # }
```
