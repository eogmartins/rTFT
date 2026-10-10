# A simulated record of stream temperature

Thirty days of hourly water temperature, with a daily cycle and a heat
wave in the third week that takes the afternoon maximum to about 24
degrees C. It is simulated, to give
[`predict_survival()`](https://eogmartins.github.io/rTFT/reference/predict_survival.md)
a record that lies below the temperatures of the example assays for most
of its length and reaches the coolest of them, 23 degrees C, on the
afternoons of the heat wave.

## Usage

``` r
stream_record
```

## Format

A tibble with 721 rows, one per hourly reading, and 3 columns:

- time:

  Numeric, minutes from the first reading.

- day:

  Numeric, the day of the record, 1 to 30.

- temp:

  Numeric, water temperature, degrees C.

## Source

Simulated by `data-raw/make-example-data.R` in the package's repository.

## See also

[static_assay](https://eogmartins.github.io/rTFT/reference/static_assay.md),
[ramp_assay](https://eogmartins.github.io/rTFT/reference/ramp_assay.md)

## Examples

``` r
range(stream_record$temp)
#> [1] 14.62 24.15
tapply(stream_record$temp, stream_record$day, max)[15:21]
#>    15    16    17    18    19    20    21 
#> 21.92 23.01 24.01 24.15 23.72 22.85 21.74 
```
