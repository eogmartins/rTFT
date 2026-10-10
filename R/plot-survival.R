# plot-survival.R

# Basic plots of what predict_survival() and average_survival() return.
# plot_survival() draws survival along the record, with the record itself above it
# when asked, and plot_risk() each period's share of the risk. Both return ggplot
# objects, so the user can refine them with further ggplot2 layers, scales, labels
# and themes. They share the panels, colours and theme of the plots of thermal death
# time curves, whose helpers they use.

# ===========================================================================
# 1. Shared helpers
# ===========================================================================

# .setting_colours() / .setting_labels() ---------------------------------------

# `n` colours for curves drawn together, one per setting of newdata. They are taken
# at `n` evenly spaced points along `.covariate_colours`, the four colours from teal
# to deep purple that the thermal death time plots use for a covariate (defined with
# those plots, in plot-tdt-ttl.R), so a plot coloured by setting and one coloured by
# a covariate's value look alike. Beyond `.max_settings` neighbouring colours are too
# close to tell apart in a legend.
.setting_colours <- function(n) grDevices::colorRampPalette(.covariate_colours)(n)
.max_settings <- 8L

# A label for each row of a table naming its setting of newdata: its covariate
# values, separated by commas, as a factor in the order the settings appear.
.setting_labels <- function(tbl, cols) {
  label <- do.call(paste, c(lapply(tbl[cols], format, trim = TRUE), sep = ", "))
  factor(label, levels = unique(label))
}


# ===========================================================================
# 2. Survival along a record
# ===========================================================================

# plot_survival() --------------------------------------------------------------

#' Plot survival along a record of temperature
#'
#' Draws the survival predicted by [predict_survival()] or [average_survival()]: the
#' probability of surviving to each time of the record, as the posterior median with
#' its credible interval. The result is a ggplot object, to be refined by the user
#' with further ggplot2 layers, scales, labels and themes.
#'
#' By default each setting of `newdata` has its own panel, following the same rule
#' as in [plot_tdt()]. There are two ways to draw the settings together in one
#' panel. With `facet = FALSE` each setting has its own colour, from teal to deep
#' purple in the order of `newdata`, and an entry in the legend, which suits a few
#' settings or a covariate that is not numeric; with more than eight settings the
#' colours are too close to tell apart, and a message says to use panels or
#' `colour_by`. With `colour_by`, the curves of every value of one
#' numeric covariate are coloured on a continuous scale from teal for its lowest
#' value to deep purple for its highest, and the intervals are left out by default,
#' since they overlap.
#'
#' With `temperature = TRUE` the record of temperature is drawn above the survival
#' curve on the same time axis, so that a drop in survival can be read against the
#' temperatures that caused it. A co-stressor the record follows through time is
#' drawn in a panel of its own between the two, named by its column, when it
#' changes along the record; one that stays constant is left out.
#'
#' @param x A list, the output of [predict_survival()] or [average_survival()].
#' @param temperature A single logical value, whether to draw the record of
#'   temperature, and of any co-stressor that changes along it, in panels above the
#'   survival curve. Defaults to `FALSE`.
#' @param facet A single logical value or the character string `"all"`. `TRUE` (the
#'   default) gives a panel per setting of `newdata`, cut down to quartile values
#'   when there are many; `"all"` a panel for every setting; `FALSE` no panels, with
#'   a colour per setting.
#' @param colour_by An optional character string naming a numeric covariate of
#'   `newdata`, to colour the curves by.
#' @param interval A single logical value, whether to draw the credible interval.
#'   Defaults to `TRUE`, and to `FALSE` with `colour_by`.
#' @param time_unit An optional character string giving the time unit for the axis
#'   label, e.g. `"min"`. The fit does not record the unit of its times, so none is
#'   shown unless given.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
#' # degrees over 60 minutes and checked every 10, in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
#'                             checks = checks)
#' fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#' # a two-day record read every hour, with a daily cycle
#' hours  <- 0:48
#' record <- data.frame(time = hours * 60, day = pmin(hours %/% 24, 1),
#'                      temp = 24 + 3 * sin(2 * pi * (hours - 9) / 24))
#'
#' pred <- predict_survival(fit, record, newdata = data.frame(mass = c(8, 12)))
#' plot_survival(pred, time_unit = "min")
#' plot_survival(pred, temperature = TRUE, time_unit = "min")
#'
#' # both masses in one panel, the record drawn once
#' plot_survival(pred, temperature = TRUE, facet = FALSE, time_unit = "min")
#' }
#' @export
plot_survival <- function(x, temperature = FALSE, facet = TRUE, colour_by = NULL,
                          interval = is.null(colour_by), time_unit = NULL) {
  tbl  <- .plot_table(x, "survival", from = "predict_survival")
  cols <- setdiff(.newdata_cols(tbl), "time")
  if (!(isTRUE(temperature) || isFALSE(temperature)))
    cli::cli_abort("{.arg temperature} must be {.code TRUE} or {.code FALSE}.",
                   class = "rTFT_bad_temperature")
  if (!is.null(colour_by) &&
      !(is.character(colour_by) && length(colour_by) == 1L && colour_by %in% cols &&
        is.numeric(tbl[[colour_by]])))
    cli::cli_abort(c(
      "{.arg colour_by} must name one numeric covariate of {.arg newdata}.",
      "x" = "You supplied {.val {colour_by}}.",
      "i" = "The covariates are {.val {cols}}."
    ), class = "rTFT_bad_colour_by")

  if (is.null(colour_by)) {
    tbl <- .choose_panels(tbl, cols, facet)
    panel_by <- if (!isFALSE(facet)) cols
  } else {
    .choose_panels(tbl, character(0), facet)             # checks `facet`
    panel_by <- if (!isFALSE(facet)) setdiff(cols, colour_by)
  }
  tbl$.group <- if (length(cols)) .curve_group(tbl, cols) else factor(1)

  # Drawn together without panels, each setting of newdata has its own colour, named
  # in the legend by its covariate values.
  n_settings <- nlevels(tbl$.group)
  by_setting <- is.null(colour_by) && isFALSE(facet) && n_settings > 1L
  if (by_setting && n_settings > .max_settings) {
    cli::cli_inform(c(
      "!" = "{n_settings} settings of {.arg newdata} are too many to tell apart by colour, so they share one.",
      "i" = "Use panels ({.code facet = TRUE}) or {.arg colour_by} for a numeric covariate."
    ), class = "rTFT_plot_colours")
    by_setting <- FALSE
  }
  if (by_setting) tbl$.setting <- .setting_labels(tbl, cols)

  # With the record above, it and survival share the time axis as rows of panels,
  # each named by what it shows, and the name stands in for the axis title:
  # temperature first, then each co-stressor that changes along the record, then
  # survival.
  g <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$time, group = .data$.group))
  if (temperature) {
    record  <- .plot_table(x, "record", from = "predict_survival")
    varying <- purrr::keep(setdiff(names(record), c("time", "temp")),
                           \(col) length(unique(record[[col]])) > 1L)
    rows <- c(temp = "Temperature (\u00b0C)", stats::setNames(varying, varying),
              "Survival")
    tbl$.row <- factor("Survival", levels = rows)
    g$data   <- tbl
    traces <- dplyr::bind_rows(purrr::map(c("temp", varying), \(col)
      tibble::tibble(time = record$time, value = record[[col]],
                     .row = factor(rows[[col]], levels = rows))))
    # The record is the same for every setting of newdata, so it is repeated in
    # each column of panels.
    if (length(panel_by)) traces <- dplyr::cross_join(unique(tbl[panel_by]), traces)
    limits <- tbl[rep(1L, 2L), , drop = FALSE]
    limits$median <- c(0, 1)
    g <- g +
      ggplot2::geom_line(data = traces, ggplot2::aes(y = .data$value, group = NULL),
                         colour = "grey30") +
      ggplot2::geom_blank(data = limits, ggplot2::aes(y = .data$median))
  }
  if (by_setting) {
    colours <- .setting_colours(n_settings)
    if (isTRUE(interval))
      g <- g + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$.setting),
        alpha = 0.2)
    g <- g +
      ggplot2::geom_line(ggplot2::aes(y = .data$median, colour = .data$.setting)) +
      ggplot2::scale_colour_manual(values = colours, name = paste(cols, collapse = ", ")) +
      ggplot2::scale_fill_manual(values = colours, name = paste(cols, collapse = ", "))
  } else if (is.null(colour_by)) {
    if (isTRUE(interval))
      g <- g + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper), fill = "#225EA8", alpha = 0.2)
    g <- g + ggplot2::geom_line(ggplot2::aes(y = .data$median), colour = "#225EA8")
  } else {
    if (isTRUE(interval))
      g <- g + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data[[colour_by]]),
        alpha = 0.2)
    g <- g +
      ggplot2::geom_line(ggplot2::aes(y = .data$median, colour = .data[[colour_by]])) +
      ggplot2::scale_colour_gradientn(colours = .covariate_colours, name = colour_by) +
      ggplot2::scale_fill_gradientn(colours = .covariate_colours, name = colour_by)
  }
  g <- g + ggplot2::labs(x = .time_label("Time", time_unit), y = "Survival") + .plot_theme()

  if (!temperature) {
    g <- g + ggplot2::scale_y_continuous(limits = c(0, 1))
    if (length(panel_by))
      g <- g + ggplot2::facet_wrap(panel_by, labeller = ggplot2::label_both)
    return(g)
  }
  layout <- stats::as.formula(paste(".row ~", if (length(panel_by))
    paste(panel_by, collapse = " + ") else "."))
  g + ggplot2::facet_grid(layout, scales = "free_y", switch = "y",
                          labeller = ggplot2::labeller(.rows = ggplot2::label_value,
                                                       .cols = ggplot2::label_both)) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme(strip.placement = "outside",
                   strip.background.y = ggplot2::element_blank(),
                   strip.text.y.left = ggplot2::element_text(colour = "grey20"))
}


# ===========================================================================
# 3. Each period's share of the risk
# ===========================================================================

# plot_risk() ------------------------------------------------------------------

#' Plot each period's share of the risk along a record
#'
#' Draws the `risk` table of [predict_survival()] or [average_survival()]: the share
#' of the hazard accrued over the record that falls in each period, as a bar at the
#' posterior median with its credible interval. It shows where along the record the
#' risk is concentrated. The result is a ggplot object, which can be refined by the
#' user with further ggplot2 layers, scales, labels and themes.
#'
#' The shares are of whatever risk there is, so they say where it falls and not how
#' large it is; read them with the `mortality` table. A period has no bar where no
#' posterior draw has any hazard along the record.
#'
#' By default each setting of `newdata` has its own panel, following the same rule
#' as in [plot_tdt()]. With `facet = FALSE` the settings are drawn side by side
#' within each period, each in its own colour, from teal to deep purple in the order
#' of `newdata`; this needs eight settings or fewer.
#'
#' @param x A list, the output of [predict_survival()] or [average_survival()]
#'   called with `by`.
#' @param facet A single logical value or the character string `"all"`. `TRUE` (the
#'   default) gives a panel per setting of `newdata`, cut down to quartile values
#'   when there are many; `"all"` a panel for every setting; `FALSE` no panels, with
#'   the settings side by side in each period.
#' @param interval A single logical value, whether to draw the credible interval.
#'   Defaults to `TRUE`.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
#' # degrees over 60 minutes and checked every 10, in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 10))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
#'                             checks = checks)
#' fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#' # a two-day record read every hour, with a daily cycle
#' hours  <- 0:48
#' record <- data.frame(time = hours * 60, day = pmin(hours %/% 24, 1),
#'                      temp = 24 + 3 * sin(2 * pi * (hours - 9) / 24))
#'
#' pred <- predict_survival(fit, record, newdata = data.frame(mass = c(8, 12)),
#'                          by = "day")
#' plot_risk(pred)
#' plot_risk(pred, facet = FALSE)
#' }
#' @export
plot_risk <- function(x, facet = TRUE, interval = TRUE) {
  if (is.list(x) && is.data.frame(x$survival) && !is.data.frame(x$risk))
    cli::cli_abort(c(
      "{.arg x} has no {.field risk} table to plot.",
      "i" = "Give {.arg by} to {.fn predict_survival}, naming the column of the record that holds the periods."
    ), class = "rTFT_no_risk")
  tbl  <- .plot_table(x, "risk", from = "predict_survival")
  cols <- setdiff(.newdata_cols(.plot_table(x, "survival", from = "predict_survival")),
                  "time")
  by   <- setdiff(.newdata_cols(tbl), cols)
  tbl  <- .choose_panels(tbl, cols, facet)
  # Periods that are not numbers keep the order they have along the record.
  if (!is.numeric(tbl[[by]])) tbl[[by]] <- factor(tbl[[by]], levels = unique(tbl[[by]]))
  tbl$.group <- if (length(cols)) .curve_group(tbl, cols) else factor(1)

  n_settings <- nlevels(tbl$.group)
  by_setting <- isFALSE(facet) && n_settings > 1L
  if (by_setting && n_settings > .max_settings)
    cli::cli_abort(c(
      "{n_settings} settings of {.arg newdata} are too many to draw side by side.",
      "i" = "Use panels ({.code facet = TRUE}), or predict at fewer settings."
    ), class = "rTFT_plot_colours")

  if (by_setting) tbl$.setting <- .setting_labels(tbl, cols)

  # A period with no risk in any draw has no share to draw. It keeps its place on
  # the axis, held by a layer that draws nothing, and gets no bar.
  shown <- tbl[!is.na(tbl$median), , drop = FALSE]
  g <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data[[by]], group = .data$.group)) +
    ggplot2::geom_blank(ggplot2::aes(y = 0))
  side_by_side <- ggplot2::position_dodge(width = 0.9)
  if (by_setting) {
    g <- g +
      ggplot2::geom_col(data = shown,
                        ggplot2::aes(y = .data$median, fill = .data$.setting),
                        position = side_by_side) +
      ggplot2::scale_fill_manual(values = .setting_colours(n_settings),
                                 name = paste(cols, collapse = ", "), drop = FALSE)
  } else {
    g <- g + ggplot2::geom_col(data = shown, ggplot2::aes(y = .data$median),
                               fill = "#225EA8")
  }
  if (isTRUE(interval))
    g <- g + ggplot2::geom_linerange(
      data = shown, ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      colour = "grey20", position = if (by_setting) side_by_side else "identity")
  g <- g +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = by, y = "Share of the risk") +
    .plot_theme()
  if (!isFALSE(facet) && length(cols))
    g <- g + ggplot2::facet_wrap(cols, labeller = ggplot2::label_both)
  g
}
