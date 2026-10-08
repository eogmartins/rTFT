# data.R - documentation of the example datasets

#' A simulated static assay of thermal tolerance
#'
#' Times to loss of equilibrium of 120 juvenile salmonids held at three constant
#' temperatures, simulated from a known model so that what a fitted model estimates
#' can be set against the values the data were made from. Twelve trials of ten fish,
#' four at each temperature (23, 25 and 27 degrees C), were checked every hour for
#' five days, so each failure is known to have happened between two checks. Fish
#' that had not failed by the end of the fifth day are survivors. The coolest
#' temperature is one a summer heat wave can reach, so that the assay informs
#' predictions for the field (see [stream_record]).
#'
#' The data are simulated, not measured. The values are plausible for a juvenile
#' rainbow trout acclimated to 14 degrees C, and stand for no particular study. The
#' model they come from is described under [ramp_assay], which shares it, and is
#' stored with the data as `attr(static_assay, "truth")`.
#'
#' @format A tibble with 120 rows, one per fish, and 7 columns:
#' \describe{
#'   \item{id}{Integer, the fish.}
#'   \item{trial}{Character, the trial the fish was in (`"S01"` to `"S12"`).}
#'   \item{path_id}{Character, the temperature treatment (`"T23"`, `"T25"`,
#'     `"T27"`), the key to its temperature trajectory.}
#'   \item{temp}{Numeric, the constant temperature of the trial, degrees C.}
#'   \item{mass}{Numeric, body mass, g.}
#'   \item{event}{Integer, 1 if the fish lost equilibrium, 0 if it survived to the
#'     end of the trial.}
#'   \item{t_upper}{Numeric, the time of the check at which the fish was found to
#'     have failed, or of the last check for a survivor, minutes.}
#' }
#' @source Simulated by `data-raw/make-example-data.R` in the package's repository.
#' @seealso [ramp_assay], [stream_record]
#' @examples
#' head(static_assay)
#'
#' # The temperature trajectory of each treatment, and the checks.
#' path <- data.frame(path_id = rep(c("T23", "T25", "T27"), each = 2),
#'                    time = rep(c(0, 7200), 3),
#'                    temp = rep(c(23, 25, 27), each = 2))
#' checks <- expand.grid(path_id = c("T23", "T25", "T27"),
#'                       time = seq(0, 7200, by = 60))
#' long <- build_hazard_data(static_assay, path, "interval",
#'                           covariates = c("mass", "trial"), checks = checks)
#' head(long)
"static_assay"

#' A simulated ramping assay of thermal tolerance, with two levels of oxygen
#'
#' Times to loss of equilibrium of 280 juvenile salmonids warmed from 14 degrees C
#' at a constant rate until they failed or the water reached 29.5 degrees C,
#' simulated from a known model. The fish were run in 28 trials of ten. Three ramp
#' rates (0.1, 0.2 and 0.3 degrees C per minute) are crossed with two levels of
#' dissolved oxygen, four trials of each combination, and a fourth, much slower
#' rate (0.025 degrees C per minute) was run in normoxia only, in four trials. The
#' fish were watched throughout, so each time of failure is known exactly.
#'
#' The slow rate allows the same data to show what the choice of ramp rates lets a
#' model estimate: the three faster rates span 3-fold, and the rates of 0.025, 0.1
#' and 0.3 span 12-fold (see the vignette *Ramping assays*).
#'
#' The data are simulated, not measured. The values are plausible for a juvenile
#' rainbow trout acclimated to 14 degrees C, and stand for no particular study.
#'
#' @section The model the data were simulated from:
#' The hazard of losing equilibrium is h0(tau) M(T), with an exponential failure
#' rate M(T) = a exp(b (T - 26)) and a Weibull baseline
#' h0(tau) = (k / s) (tau / s)^(k - 1), with times in minutes. The values, stored as
#' `attr(ramp_assay, "truth")`, are b = 1.7 per degree C, k = 1.5, s = 60 minutes,
#' and log a = -2.13 at 26 degrees C for a 15 g fish in normoxia. log a varies among
#' trials with a standard deviation of 0.3, rises by 0.03 per gram of body mass and
#' by 0.3 per mg/L of dissolved oxygen below 10. At a constant temperature half the
#' fish then fail in about 1 hour at 27 degrees C, 10 hours at 25 and 4 days at 23.
#' [static_assay] was simulated from the same model.
#'
#' @format A tibble with 280 rows, one per fish, and 8 columns:
#' \describe{
#'   \item{id}{Integer, the fish.}
#'   \item{trial}{Character, the trial the fish was in (`"R01"` to `"R28"`).}
#'   \item{path_id}{Character, the ramp (`"r0.025"`, `"r0.1"`, `"r0.2"`, `"r0.3"`),
#'     the key to its temperature trajectory.}
#'   \item{rate}{Numeric, the ramp rate, degrees C per minute.}
#'   \item{do}{Numeric, dissolved oxygen in the trial, mg/L (10 or 5).}
#'   \item{mass}{Numeric, body mass, g.}
#'   \item{event}{Integer, 1 if the fish lost equilibrium, 0 if it survived to the
#'     end of the ramp.}
#'   \item{t_event}{Numeric, the time of failure, or of the end of the ramp for a
#'     survivor, minutes from its start.}
#' }
#' @source Simulated by `data-raw/make-example-data.R` in the package's repository.
#' @seealso [static_assay], [stream_record]
#' @examples
#' head(ramp_assay)
#'
#' # The four ramps, from 14 to 29.5 degrees C.
#' rates <- c(0.025, 0.1, 0.2, 0.3)
#' path  <- ramp_path(c("r0.025", "r0.1", "r0.2", "r0.3"), T_start = 14, rate = rates,
#'                    t_max = 15.5 / rates)
#' # No grid is given: the intervals follow the ramp in steps of 0.1 degrees C.
#' long  <- build_hazard_data(ramp_assay, path, "continuous",
#'                            covariates = c("mass", "do", "trial"))
#' head(long)
"ramp_assay"

#' A simulated record of stream temperature
#'
#' Thirty days of hourly water temperature, with a daily cycle and a heat wave in
#' the third week that takes the afternoon maximum to about 24 degrees C. It is
#' simulated, to give [predict_survival()] a record that lies below the
#' temperatures of the example assays for most of its length and reaches the coolest
#' of them, 23 degrees C, on the afternoons of the heat wave.
#'
#' @format A tibble with 721 rows, one per hourly reading, and 3 columns:
#' \describe{
#'   \item{time}{Numeric, minutes from the first reading.}
#'   \item{day}{Numeric, the day of the record, 1 to 30.}
#'   \item{temp}{Numeric, water temperature, degrees C.}
#' }
#' @source Simulated by `data-raw/make-example-data.R` in the package's repository.
#' @seealso [static_assay], [ramp_assay]
#' @examples
#' range(stream_record$temp)
#' tapply(stream_record$temp, stream_record$day, max)[15:21]
"stream_record"
