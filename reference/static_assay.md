# A simulated static assay of thermal tolerance

Times to loss of equilibrium of 120 juvenile salmonids held at three
constant temperatures, simulated from a known model so that what a
fitted model estimates can be set against the values the data were made
from. Twelve trials of ten fish, four at each temperature (23, 25 and 27
degrees C), were checked every hour for five days, so each failure is
known to have happened between two checks. Fish that had not failed by
the end of the fifth day are survivors. The coolest temperature is one a
summer heat wave can reach, so that the assay informs predictions for
the field (see
[stream_record](https://eogmartins.github.io/rTFT/reference/stream_record.md)).

## Usage

``` r
static_assay
```

## Format

A tibble with 120 rows, one per fish, and 7 columns:

- id:

  Integer, the fish.

- trial:

  Character, the trial the fish was in (`"S01"` to `"S12"`).

- path_id:

  Character, the temperature group (`"T23"`, `"T25"`, `"T27"`), the key
  to its temperature trajectory.

- temp:

  Numeric, the constant temperature of the trial, degrees C.

- mass:

  Numeric, body mass, g.

- event:

  Integer, 1 if the fish lost equilibrium, 0 if it survived to the end
  of the trial.

- t_upper:

  Numeric, the time of the check at which the fish was found to have
  failed, or of the last check for a survivor, minutes.

## Source

Simulated by `data-raw/make-example-data.R` in the package's repository.

## Details

The data are simulated, not measured. The values are plausible for a
juvenile rainbow trout acclimated to 14 degrees C, and stand for no
particular study. The model they come from is described under
[ramp_assay](https://eogmartins.github.io/rTFT/reference/ramp_assay.md),
which shares it, and is stored with the data as
`attr(static_assay, "truth")`.

## See also

[ramp_assay](https://eogmartins.github.io/rTFT/reference/ramp_assay.md),
[stream_record](https://eogmartins.github.io/rTFT/reference/stream_record.md)

## Examples

``` r
head(static_assay)
#> # A tibble: 6 × 7
#>      id trial path_id  temp  mass event t_upper
#>   <int> <chr> <chr>   <dbl> <dbl> <int>   <dbl>
#> 1     1 S01   T23        23  19.9     0    7200
#> 2     2 S01   T23        23  14.3     0    7200
#> 3     3 S01   T23        23   8       0    7200
#> 4     4 S01   T23        23  16.6     1    2040
#> 5     5 S01   T23        23  29.5     1    2400
#> 6     6 S01   T23        23   9       1    5820

# The temperature trajectory of each group, and the checks.
path <- data.frame(path_id = rep(c("T23", "T25", "T27"), each = 2),
                   time = rep(c(0, 7200), 3),
                   temp = rep(c(23, 25, 27), each = 2))
checks <- expand.grid(path_id = c("T23", "T25", "T27"),
                      time = seq(0, 7200, by = 60))
long <- build_hazard_data(static_assay, path, "interval",
                          covariates = c("mass", "trial"), checks = checks)
head(long)
#> # A tibble: 6 × 12
#>      id  mass trial path_id   tau exposure    Ta failed    Tr log_tau tau_ref
#>   <int> <dbl> <chr> <chr>   <dbl>    <dbl> <dbl>  <int> <dbl>   <dbl>   <dbl>
#> 1     1  19.9 S01   T23        30       60    23      0    25    3.40     360
#> 2     1  19.9 S01   T23        90       60    23      0    25    4.50     360
#> 3     1  19.9 S01   T23       150       60    23      0    25    5.01     360
#> 4     1  19.9 S01   T23       210       60    23      0    25    5.35     360
#> 5     1  19.9 S01   T23       270       60    23      0    25    5.60     360
#> 6     1  19.9 S01   T23       330       60    23      0    25    5.80     360
#> # ℹ 1 more variable: log_exposure <dbl>
```
