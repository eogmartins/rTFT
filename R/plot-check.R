# plot-check.R

# Plots of the two checks of a fitted model: plot_check() draws what
# check_survival() returns, and plot_residuals() what hazard_residuals() returns.
# Both return ggplot objects, so the user can refine them with further ggplot2 layers,
# scales, labels and themes, and share the theme of the package's other plots.

# ===========================================================================
# 1. The posterior predictive check of survival
# ===========================================================================

# plot_check() -----------------------------------------------------------------

#' Plot the check of a fitted model against observed survival
#'
#' Draws the output of [check_survival()]: for each group, the proportion of
#' individuals surviving in datasets simulated from the fitted model, as its median
#' and band, with the proportion observed to survive as a stepped line over it.
#' Where the model describes the data, the stepped line stays within the band. The
#' result is a ggplot object, which can be refined by the user with further ggplot2
#' layers, scales, labels and themes.
#'
#' @param x A tibble, the output of [check_survival()].
#' @param time_unit An optional character string giving the time unit for the axis
#'   label, e.g. `"min"`. The fit does not record the unit of its times, so none is
#'   shown unless given.
#' @return A ggplot object, with a panel per group.
#' @examples
#' \dontrun{
#' plot_check(check_survival(fit), time_unit = "min")
#' }
#' @export
plot_check <- function(x, time_unit = NULL) {
  needed <- c("path_id", "time", "observed", "median", "lower", "upper")
  if (!is.data.frame(x) || !all(needed %in% names(x)))
    cli::cli_abort("{.arg x} must be the output of {.fn check_survival}.",
                   class = "rTFT_bad_plot_input")
  panel_by <- setdiff(names(x), c("n", "Ta", "prob", needed[-1]))
  ggplot2::ggplot(x, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "#225EA8", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median), colour = "#225EA8") +
    ggplot2::geom_step(ggplot2::aes(y = .data$observed), colour = "grey20") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = .time_label("Time", time_unit), y = "Survival") +
    ggplot2::facet_wrap(panel_by, scales = "free_x", labeller = ggplot2::label_both) +
    .plot_theme()
}


# ===========================================================================
# 2. Residuals by individual
# ===========================================================================

# plot_residuals() -------------------------------------------------------------

#' Plot the residuals of a fitted model
#'
#' Draws the output of [hazard_residuals()]. With `against = NULL`, the residuals
#' are set against the values expected of a uniform distribution, which is what
#' they follow when the model describes the data: the points then lie along the
#' diagonal. With `against` naming a column, the residuals are drawn against it, to
#' look for a pattern the model has missed: points for a numeric variable, such as
#' body mass, and boxes for a grouping, such as the trial or the treatment. A
#' horizontal line marks 0.5, the middle of the residuals under the model; residuals
#' below it are individuals that failed earlier than the model predicts. The time
#' of failure is not a variable to draw them against (see [hazard_residuals()]).
#' The result is a ggplot object, which can be refined by the user with further
#' ggplot2 layers, scales, labels and themes.
#'
#' @param x A tibble, the output of [hazard_residuals()].
#' @param against An optional character string naming a column of `x` to draw the
#'   residuals against. `NULL` (the default) compares them with the uniform
#'   distribution.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' res <- hazard_residuals(fit)
#' plot_residuals(res)
#' plot_residuals(res, against = "mass")
#' plot_residuals(res, against = "trial")
#' }
#' @export
plot_residuals <- function(x, against = NULL) {
  if (!is.data.frame(x) || !"residual" %in% names(x))
    cli::cli_abort("{.arg x} must be the output of {.fn hazard_residuals}.",
                   class = "rTFT_bad_plot_input")
  if (!is.null(against) &&
      !(is.character(against) && length(against) == 1L && against %in% names(x) &&
        against != "residual"))
    cli::cli_abort(c(
      "{.arg against} must name one column of {.arg x} to draw the residuals against.",
      "x" = "You supplied {.val {against}}.",
      "i" = "The columns are {.val {setdiff(names(x), 'residual')}}."
    ), class = "rTFT_bad_against")

  if (is.null(against)) {
    n <- nrow(x)
    qq <- tibble::tibble(expected = (seq_len(n) - 0.5) / n, residual = sort(x$residual))
    return(
      ggplot2::ggplot(qq, ggplot2::aes(x = .data$expected, y = .data$residual)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey50",
                             linetype = "dashed") +
        ggplot2::geom_point(colour = "#225EA8") +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "Expected under the model", y = "Residual") +
        .plot_theme())
  }
  g <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[against]], y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0.5, colour = "grey50", linetype = "dashed")
  g <- if (is.numeric(x[[against]])) {
    g + ggplot2::geom_point(colour = "#225EA8")
  } else {
    g + ggplot2::geom_boxplot(fill = "#225EA8", alpha = 0.2, colour = "grey20",
                              outlier.shape = NA) +
      ggplot2::geom_point(colour = "#225EA8",
                          position = ggplot2::position_jitter(width = 0.15, height = 0))
  }
  g + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = against, y = "Residual") +
    .plot_theme()
}
