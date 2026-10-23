# A table of failure-rate models to fit and compare

Builds a table with one row per model, a failure-rate function with a
baseline and the formulas of its parameters, for
[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md)
to fit. Each model is assembled once here, so a misspelt name or a
formula for a parameter no model has is reported at once, before
anything is fitted.

## Usage

``` r
build_hazard_models(
  m_form,
  baseline = "mem",
  dpar_formulas = list(),
  costressors = list(),
  cross = TRUE,
  model_names = NULL
)
```

## Arguments

- m_form:

  A character vector, the failure-rate functions, each a key of
  [M_catalogue](https://eogmartins.github.io/rTFT/reference/M_catalogue.md).

- baseline:

  A character vector, the baselines, each a key of
  [h0_catalogue](https://eogmartins.github.io/rTFT/reference/h0_catalogue.md).
  Defaults to `"mem"`, the memoryless baseline.

- dpar_formulas:

  A named list of one-sided formulas, one per parameter that depends on
  covariates or has a group-level effect, as in
  [`build_hazard_spec()`](https://eogmartins.github.io/rTFT/reference/build_hazard_spec.md).
  Each is given to the models that have the parameter.

- costressors:

  An optional list of
  [costressor](https://eogmartins.github.io/rTFT/reference/costressor.md)
  objects, given to every model. Defaults to none.

- cross:

  A single logical value. `TRUE` (the default) combines every
  failure-rate function with every baseline; `FALSE` pairs them in the
  order given, a single value being repeated.

- model_names:

  An optional character vector, a name for each model, in the order of
  the rows. By default a model is named by its failure-rate function and
  its baseline, as in `"exp_weibull"`.

## Value

A data frame with one row per model and the columns `model` (its name),
`m_form`, `baseline`, `dpar_formulas` and `costressors`, the last two
holding a list for each model.

## Details

With `cross = TRUE` every failure-rate function is combined with every
baseline. With `cross = FALSE` they are paired in the order given, for a
set of models that is not every combination. The table is an ordinary
data frame: two tables can be bound by rows, a row dropped, or the
formulas of one model replaced, as long as each model keeps a name of
its own.

The models differ in the parameters they have (a threshold function has
`Tc`, a Weibull baseline has `logshape`), so a formula in
`dpar_formulas` is given to each model that has the parameter and to no
other. A formula for `loga`, which every failure-rate function has, is
therefore shared by all the models, and one for `logshape` reaches only
those with a Weibull or a log-logistic baseline.

## See also

[`fit_hazard_models()`](https://eogmartins.github.io/rTFT/reference/fit_hazard_models.md),
which fits the models of the table.

## Examples

``` r
# Two failure-rate functions with each of two baselines, all with a random
# intercept per trial on log a. The Weibull shape also depends on body mass, in
# the two models that have a Weibull baseline.
build_hazard_models(
  m_form = c("exp", "zeropower_threshold"), baseline = c("mem", "weibull"),
  dpar_formulas = list(loga = ~ 1 + (1 | trial), logshape = ~ 1 + mass))
#> # A tibble: 4 × 5
#>   model                       m_form          baseline dpar_formulas costressors
#>   <chr>                       <chr>           <chr>    <list>        <list>     
#> 1 exp_mem                     exp             mem      <named list>  <list [0]> 
#> 2 exp_weibull                 exp             weibull  <named list>  <list [0]> 
#> 3 zeropower_threshold_mem     zeropower_thre… mem      <named list>  <list [0]> 
#> 4 zeropower_threshold_weibull zeropower_thre… weibull  <named list>  <list [0]> 

# A set that is not every combination, with names of its own.
build_hazard_models(
  m_form = c("exp", "exp", "zeropower_threshold"),
  baseline = c("mem", "llogis", "mem"), cross = FALSE,
  model_names = c("Exponential, memoryless", "Exponential, log-logistic",
            "Zero-power threshold, memoryless"))
#> # A tibble: 3 × 5
#>   model                            m_form     baseline dpar_formulas costressors
#>   <chr>                            <chr>      <chr>    <list>        <list>     
#> 1 Exponential, memoryless          exp        mem      <list [0]>    <list [0]> 
#> 2 Exponential, log-logistic        exp        llogis   <list [0]>    <list [0]> 
#> 3 Zero-power threshold, memoryless zeropower… mem      <list [0]>    <list [0]> 
```
