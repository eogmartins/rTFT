# plot-tdt-ttl.R

# Basic plots of what derive_tdt_ttl() returns: the thermal death time curves and
# the thermal tolerance landscape. plot_tdt() and plot_ttl() return ggplot objects,
# so the user refines them with further ggplot2 layers, scales, labels and themes.
# plot_ttl_3d() returns a plotly object, a surface to rotate.

# ===========================================================================
# 1. Shared helpers
# ===========================================================================

# .plot_table() / .newdata_cols() / .time_label() ------------------------------

# The columns derive_tdt_ttl() and predict_survival() write themselves. Any other
# column of one of their tables is a covariate carried over from `newdata`.
.derive_cols <- c("quantity", "p", "t_ref", "Ta", "tau", "mean", "median", "lower",
                  "upper", "prob", "n_draws", "n_unreached")

# Take one table out of the output of derive_tdt_ttl() (or of the function named in
# `from`), or stop when `x` is not such an output.
.plot_table <- function(x, which, from = "derive_tdt_ttl") {
  if (!is.list(x) || !is.data.frame(x[[which]]))
    cli::cli_abort(c(
      "{.arg x} must be the output of {.fn {from}}.",
      "x" = "It has no {.field {which}} table."
    ), class = "rTFT_bad_plot_input")
  x[[which]]
}

# The covariate columns of a table from derive_tdt_ttl().
.newdata_cols <- function(tbl) setdiff(names(tbl), .derive_cols)

# An axis label for a time, with its unit when the user gives one. The fit does not
# record the unit of its times, so none is assumed.
.time_label <- function(label, time_unit) {
  if (is.null(time_unit)) label else paste0(label, " (", time_unit, ")")
}

# .nearest_quartiles() / .choose_panels() --------------------------------------

# The values of `v` nearest its minimum, lower quartile, median, upper quartile and
# maximum: up to five values that are in `v`, so that each is a setting the
# quantities were derived at.
.nearest_quartiles <- function(v) {
  v <- sort(unique(v))
  q <- stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  unique(vapply(q, \(x) v[which.min(abs(v - x))], numeric(1)))
}

# Choose which settings of newdata get a panel. Each row of newdata has its own
# curve, an exact result at those covariate values, so rows cannot be pooled into
# groups; they can only be shown or left out. Up to six settings are all shown. With
# more, each numeric covariate with more than five values is cut down to the values
# nearest its quartiles, and a message says which were kept. `facet = "all"` shows
# every setting and `facet = FALSE` shows them all in one panel, for the user to
# arrange.
.choose_panels <- function(tbl, cols, facet) {
  if (!(isTRUE(facet) || isFALSE(facet) || identical(facet, "all")))
    cli::cli_abort(c(
      "{.arg facet} must be {.code TRUE}, {.code FALSE} or {.val all}.",
      "x" = "You supplied {.val {facet}}."
    ), class = "rTFT_bad_facet")
  if (!length(cols) || !isTRUE(facet) || nrow(unique(tbl[cols])) <= 6)
    return(tbl)
  for (col in cols) {
    values <- unique(tbl[[col]])
    if (!is.numeric(values) || length(values) <= 5) next
    keep <- .nearest_quartiles(values)
    tbl  <- tbl[tbl[[col]] %in% keep, , drop = FALSE]
    cli::cli_inform(c(
      "i" = "Showing {length(keep)} of {length(values)} values of {.field {col}} (its quartiles): {.val {keep}}.",
      "i" = "Use {.code facet = \"all\"} for a panel each."
    ), class = "rTFT_plot_panels")
  }
  tbl
}

# The ColorBrewer RdYlBu colours, from red through pale yellow to blue, for
# survival: red where most individuals have failed, blue where most survive, and the
# pale middle on a survival of one half. A diverging scale suits survival because
# one half is a meaningful midpoint, and as a fill its pale middle stays visible.
.rdylbu <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF",
             "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

# `n` colours from cool blue to warm red, for lines that have an order: the
# fractions failed of the thermal death time curves, fewest to most, and the
# temperatures of the survival curves, cool to hot. The ramp runs between the two
# ends of RdYlBu without passing through its pale middle, which would vanish as a
# line. So on every plot blue reads as safe and red as lethal: the curve for a
# fraction failed of 0.1 is the landscape's contour at a survival of 0.9, and both
# are blue.
.cool_warm <- function(n) grDevices::colorRampPalette(c("#2C7BB6", "#D7191C"))(n)

# A matching colour and fill scale of `n` cool-to-warm colours, for the lines and
# their interval ribbons.
.cool_warm_scales <- function(n, name) {
  list(ggplot2::scale_colour_manual(values = .cool_warm(n), name = name),
       ggplot2::scale_fill_manual(values = .cool_warm(n), name = name))
}

# The colours for curves coloured by a covariate, from its lowest value to its
# highest: teal, through blue, to deep purple. Both the hue and the lightness
# change, so neighbouring curves differ and the order survives in greyscale. The
# colours stay on the cool side of the other plots' palette and use no red, since a
# covariate such as body mass is neither safe nor lethal.
.covariate_colours <- c("#7FCDBB", "#41B6C4", "#225EA8", "#3F007D")

# Axis labels for a time on a log scale as plain numbers (0.1, 1, 10, 100), not in
# scientific notation.
.plain_number <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

# The look shared by the plots: a light theme, with the legend above the panels.
.plot_theme <- function() {
  list(ggplot2::theme_light(), ggplot2::theme(legend.position = "top"))
}

# One group per curve: a curve is one setting of newdata at one further value (a
# fraction failing, a temperature), so that lines are not joined across settings.
.curve_group <- function(tbl, cols) interaction(tbl[cols], drop = TRUE)

# .cell_edges() ----------------------------------------------------------------

# The lower and upper edge of the cell around each grid value, halfway to its
# neighbours, measured on a log scale when `log` is TRUE. The landscape's grids need
# not be evenly spaced (the default times are evenly spaced only on a log scale), so
# each cell is drawn from its own edges.
.cell_edges <- function(v, log = FALSE) {
  u <- sort(unique(v))
  w <- if (log) log(u) else u
  n <- length(w)
  if (n == 1L) {
    lo <- w - 0.5
    hi <- w + 0.5
  } else {
    mid <- (w[-1] + w[-n]) / 2
    lo  <- c(w[1] - (mid[1] - w[1]), mid)
    hi  <- c(mid, w[n] + (w[n] - mid[n - 1]))
  }
  if (log) {
    lo <- exp(lo)
    hi <- exp(hi)
  }
  at <- match(v, u)
  list(lo = lo[at], hi = hi[at])
}


# .contour_labels() ------------------------------------------------------------

# Where to label each contour line of the landscape map: for every setting of
# newdata and every survival level, one point on that level's line, with the level
# as its label. The line is looked for first along elapsed time held fixed: at a
# time near the middle of those at which the level lies between the lowest and
# highest survival across temperatures, the temperature is interpolated between the
# two grid temperatures either side of the level. On a coarse grid a line can pass
# between two grid times at every temperature, so it is then looked for along
# temperature held fixed, interpolating the time (on a log scale when the time axis
# is one). A level whose line does not cross the landscape gets no label.
.contour_labels <- function(tbl, cols, survival, log_time = FALSE) {
  # The grid values of `fixed` at which the level lies within the range of survival.
  crossed_at <- function(d, fixed, level) {
    purrr::keep(sort(unique(d[[fixed]])), \(value) {
      m <- d$median[d[[fixed]] == value]
      length(m) > 1 && min(m) <= level && level <= max(m)
    })
  }
  middle <- function(v) v[ceiling(length(v) / 2)]
  settings <- if (length(cols)) split(tbl, .curve_group(tbl, cols)) else list(tbl)

  out <- dplyr::bind_rows(purrr::map(settings, \(d) {
    dplyr::bind_rows(purrr::map(survival, \(level) {
      times <- crossed_at(d, "tau", level)
      temps <- crossed_at(d, "Ta", level)
      if (length(times)) {
        tau  <- middle(times)
        line <- d[d$tau == tau, , drop = FALSE]
        Ta   <- stats::approx(line$median, line$Ta, xout = level, ties = mean)$y
      } else if (length(temps)) {
        Ta   <- middle(temps)
        line <- d[d$Ta == Ta, , drop = FALSE]
        time <- if (log_time) log(line$tau) else line$tau
        tau  <- stats::approx(line$median, time, xout = level, ties = mean)$y
        if (log_time) tau <- exp(tau)
      } else {
        return(NULL)
      }
      dplyr::bind_cols(d[1, cols, drop = FALSE],
                       tibble::tibble(tau = tau, Ta = Ta, label = format(level)))
    }))
  }))
  if (!nrow(out)) out <- tibble::tibble(tau = numeric(0), Ta = numeric(0), label = character(0))
  out
}


# ===========================================================================
# 2. Thermal death time curves
# ===========================================================================

# plot_tdt() -------------------------------------------------------------------

#' Plot the thermal death time curves
#'
#' Draws the thermal death time curves from [derive_tdt_ttl()]: the time until a
#' fraction `p` of individuals has failed, against constant temperature, as the
#' posterior median with its credible interval. The result is a ggplot object, which
#' can be refined by the user with further ggplot2 layers, scales, labels and themes.
#'
#' By default each fraction `p` has its own colour, from blue for the fewest failed
#' to red for the most, and each setting of `newdata` its own panel. Each setting has
#' its own curve, so settings are shown or left out, never averaged: with more than 
#' six, a numeric covariate with more than five values is cut down to the values 
#' nearest its minimum, quartiles and maximum, and a message says which were kept.
#'
#' With `colour_by`, the curves of every value of one numeric covariate are drawn
#' together, coloured on a continuous scale from teal for its lowest value to deep
#' purple for its highest, with a panel per fraction `p`. This shows how the curve
#' shifts with the covariate. The intervals are then left out by default, since they
#' overlap. Where the covariate is itself a temperature, such as an acclimation
#' temperature, a cool-to-warm scale may read better, and can be added to the result:
#' `+ ggplot2::scale_colour_gradient(low = "#2C7BB6", high = "#D7191C")`.
#'
#' Time is on a log scale, on which an exponential failure rate gives straight
#' lines. Where a share of the posterior never reaches the fraction, the interval
#' runs to the top of the panel, and where the median itself is never reached the
#' line stops.
#'
#' @param x A list, the output of [derive_tdt_ttl()].
#' @param facet A single logical value or the character string `"all"`. `TRUE` (the
#'   default) gives a panel per setting of `newdata`, cut down to quartile values
#'   when there are many; `"all"` a panel for every setting; `FALSE` no panels.
#' @param colour_by An optional character string naming a numeric covariate of
#'   `newdata`, to colour the curves by instead of by `p`.
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
#'
#' tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
#' plot_tdt(tol, time_unit = "min")
#'
#' # how the curves shift with body mass
#' tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = seq(5, 15, by = 0.5)),
#'                       t_ref = 60)
#' plot_tdt(tol, colour_by = "mass", time_unit = "min")
#' }
#' @importFrom ggplot2 .data
#' @export
plot_tdt <- function(x, facet = TRUE, colour_by = NULL, interval = is.null(colour_by),
                     time_unit = NULL) {
  tbl  <- .plot_table(x, "tdt_curves")
  cols <- .newdata_cols(tbl)
  if (!is.null(colour_by) &&
      !(is.character(colour_by) && length(colour_by) == 1L && colour_by %in% cols &&
        is.numeric(tbl[[colour_by]])))
    cli::cli_abort(c(
      "{.arg colour_by} must name one numeric covariate of {.arg newdata}.",
      "x" = "You supplied {.val {colour_by}}.",
      "i" = "The covariates are {.val {cols}}."
    ), class = "rTFT_bad_colour_by")

  # A median that is never reached is left out of the line, which stops there.
  tbl$median[!is.finite(tbl$median)] <- NA_real_

  if (is.null(colour_by)) {
    tbl <- .choose_panels(tbl, cols, facet)
    tbl$p <- factor(tbl$p)
    tbl$.group <- .curve_group(tbl, c("p", cols))
    panel_by <- if (!isFALSE(facet)) cols
    colour <- "p"
    scales <- .cool_warm_scales(nlevels(tbl$p), "Fraction failed")
  } else {
    .choose_panels(tbl, character(0), facet)             # checks `facet`
    tbl$.group <- .curve_group(tbl, c("p", cols))
    panel_by <- if (!isFALSE(facet)) c("p", setdiff(cols, colour_by))
    colour <- colour_by
    scales <- list(
      ggplot2::scale_colour_gradientn(colours = .covariate_colours, name = colour_by),
      ggplot2::scale_fill_gradientn(colours = .covariate_colours, name = colour_by))
  }

  g <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$Ta, group = .data$.group))
  if (isTRUE(interval))
    g <- g + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data[[colour]]),
      alpha = 0.2)
  g <- g +
    ggplot2::geom_line(ggplot2::aes(y = .data$median, colour = .data[[colour]]),
                       na.rm = TRUE) +
    scales +
    ggplot2::scale_y_log10(labels = .plain_number) +
    ggplot2::labs(x = "Temperature (\u00b0C)",
                  y = .time_label("Time to failure", time_unit)) +
    .plot_theme()
  if (length(panel_by))
    g <- g + ggplot2::facet_wrap(panel_by, labeller = ggplot2::label_both)
  g
}


# ===========================================================================
# 3. Thermal tolerance landscape
# ===========================================================================

# plot_ttl() -------------------------------------------------------------------

#' Plot the thermal tolerance landscape
#'
#' Draws the thermal tolerance landscape from [derive_tdt_ttl()]: survival over
#' constant temperature and elapsed exposure time. The result is a ggplot object, to
#' be refined by the user with further ggplot2 layers, scales, labels and themes.
#'
#' With `type = "surface"` the whole landscape is drawn as a map, elapsed time across
#' and temperature up, coloured by the posterior median survival from red (most have
#' failed) through pale yellow (a survival of one half) to blue (most survive), with
#' contour lines at the survival levels in `survival`, each labelled with its level.
#' The credible interval is not shown on it. With `type = "curves"` survival is
#' drawn against elapsed time at a few temperatures, each as the posterior median
#' with its credible interval.
#'
#' Elapsed time is on a log scale, unless the times include zero. The panels follow
#' the same rule as in [plot_tdt()], except that the map needs a panel for each
#' setting of `newdata`, so `facet = FALSE` is refused for it when there are several.
#'
#' @param x A list, the output of [derive_tdt_ttl()].
#' @param type A character string, either `"surface"` (the default) for the
#'   landscape as a map, or `"curves"` for survival curves at a few temperatures.
#' @param facet A single logical value or the character string `"all"`. `TRUE` (the
#'   default) gives a panel per setting of `newdata`, cut down to quartile values
#'   when there are many; `"all"` a panel for every setting; `FALSE` no panels.
#' @param Ta An optional numeric vector of the temperatures to draw for
#'   `type = "curves"`, each replaced by the nearest temperature of the landscape's
#'   grid. `NULL` (the default) takes the grid temperatures nearest its minimum,
#'   quartiles and maximum.
#' @param survival A numeric vector of survival levels between 0 and 1, at which to
#'   draw the contour lines for `type = "surface"`.
#' @param labels A single logical value, whether to label each contour line of
#'   `type = "surface"` with its survival level.
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
#'
#' tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
#' plot_ttl(tol, time_unit = "min")
#' plot_ttl(tol, type = "curves", Ta = c(22, 25, 28), time_unit = "min")
#' }
#' @export
plot_ttl <- function(x, type = c("surface", "curves"), facet = TRUE, Ta = NULL,
                     survival = c(0.1, 0.5, 0.9), labels = TRUE, time_unit = NULL) {
  type <- match.arg(type)
  tbl  <- .plot_table(x, "ttl_surface")
  cols <- .newdata_cols(tbl)
  tbl  <- .choose_panels(tbl, cols, facet)
  panel_by <- if (!isFALSE(facet)) cols
  log_time <- all(tbl$tau > 0)
  x_lab    <- .time_label("Elapsed exposure time", time_unit)

  if (type == "surface") {
    # A map has one cell per temperature and time, so two settings in one panel
    # would be drawn on top of each other, with only the last one seen.
    if (isFALSE(facet) && length(cols) && nrow(unique(tbl[cols])) > 1)
      cli::cli_abort(c(
        "The landscape as a map shows one setting of {.arg newdata} per panel.",
        "x" = "{.code facet = FALSE} would draw {nrow(unique(tbl[cols]))} settings on top of each other.",
        "i" = "Use {.code facet = TRUE}, or {.code type = \"curves\"}, which can show several settings in one panel."
      ), class = "rTFT_surface_needs_panels")
    tau_edge <- .cell_edges(tbl$tau, log = log_time)
    Ta_edge  <- .cell_edges(tbl$Ta)
    tbl$tau_lo <- tau_edge$lo
    tbl$tau_hi <- tau_edge$hi
    tbl$Ta_lo  <- Ta_edge$lo
    tbl$Ta_hi  <- Ta_edge$hi
    g <- ggplot2::ggplot(tbl) +
      ggplot2::geom_rect(ggplot2::aes(xmin = .data$tau_lo, xmax = .data$tau_hi,
                                      ymin = .data$Ta_lo, ymax = .data$Ta_hi,
                                      fill = .data$median)) +
      ggplot2::scale_fill_gradientn(colours = .rdylbu, limits = c(0, 1),
                                    name = "Survival") +
      ggplot2::labs(x = x_lab, y = "Temperature (\u00b0C)")
    # Contour lines need at least two temperatures and two times to run between.
    if (length(unique(tbl$tau)) > 1 && length(unique(tbl$Ta)) > 1 && length(survival)) {
      g <- g + ggplot2::geom_contour(
        ggplot2::aes(x = .data$tau, y = .data$Ta, z = .data$median),
        breaks = survival, colour = "grey20", linewidth = 0.3)
      if (isTRUE(labels)) {
        at <- .contour_labels(tbl, cols, survival, log_time)
        if (nrow(at))
          # Small plain text, set just above and to the right of its point. The
          # contours of a tolerance landscape fall from left to right (a higher
          # temperature is survived for a shorter time), so that corner is clear of
          # the line.
          g <- g + ggplot2::geom_text(
            data = at,
            ggplot2::aes(x = .data$tau, y = .data$Ta, label = .data$label),
            size = 2.5, colour = "grey20", hjust = -0.15, vjust = -0.4)
      }
    }
  } else {
    grid <- sort(unique(tbl$Ta))
    show <- if (is.null(Ta)) .nearest_quartiles(grid)
            else unique(vapply(Ta, \(v) grid[which.min(abs(grid - v))], numeric(1)))
    tbl  <- tbl[tbl$Ta %in% show, , drop = FALSE]
    tbl$Ta <- factor(format(tbl$Ta, digits = 3), levels = format(sort(show), digits = 3))
    tbl$.group <- .curve_group(tbl, c("Ta", cols))
    g <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$tau, group = .data$.group)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                        fill = .data$Ta), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = .data$median, colour = .data$Ta)) +
      .cool_warm_scales(nlevels(tbl$Ta), "Temperature (\u00b0C)") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = x_lab, y = "Survival")
  }

  if (log_time) g <- g + ggplot2::scale_x_log10(labels = .plain_number)
  g <- g + .plot_theme()
  if (length(panel_by))
    g <- g + ggplot2::facet_wrap(panel_by, labeller = ggplot2::label_both)
  g
}


# ===========================================================================
# 4. Thermal tolerance landscape in three dimensions
# ===========================================================================

# .has_plotly() ----------------------------------------------------------------

# Whether the plotly package is installed. plot_ttl_3d() is the one function that
# needs it, so it is a suggested package rather than a required one.
.has_plotly <- function() requireNamespace("plotly", quietly = TRUE)

# .log_ticks() -----------------------------------------------------------------

# Round tick values for a log axis covering `range`: 1, 2 and 5 times each power of
# ten that fall within it (2, 5, 10, 20, 50, ...). plotly's own choice for a log
# axis spanning less than two powers of ten labels every multiple with a single
# digit, so that 20, 30 and 40 read as "2 3 4". Where that gives fewer than four
# ticks, 3 times each power of ten is added as well, and on a range narrower still,
# every whole multiple.
.log_ticks <- function(range) {
  powers <- 10^seq(floor(log10(range[1])), ceiling(log10(range[2])))
  inside <- function(steps) {
    v <- sort(as.vector(outer(steps, powers)))
    v[v >= range[1] & v <= range[2]]
  }
  ticks <- inside(c(1, 2, 5))
  if (length(ticks) < 4) ticks <- inside(c(1, 2, 3, 5))
  if (length(ticks) < 3) ticks <- inside(1:9)
  ticks
}

# plot_ttl_3d() ----------------------------------------------------------------

#' Plot the thermal tolerance landscape as a surface to rotate
#'
#' Draws the thermal tolerance landscape from [derive_tdt_ttl()] in three
#' dimensions: elapsed exposure time and constant temperature across the floor, and
#' the posterior median survival as the height. The figure is interactive, so it can
#' be rotated, zoomed and read by hovering, in a viewer, a browser or an HTML
#' document.
#'
#' It needs the plotly package, and returns a plotly object rather than a ggplot, so
#' it is refined with plotly's own functions (for example [plotly::layout()]), not
#' with ggplot2 layers. A three-dimensional surface hides part of itself in a static
#' figure; for print, [plot_ttl()] shows the same landscape as a map.
#'
#' One setting of `newdata` is drawn at a time. Elapsed time is on a log scale,
#' unless the times include zero.
#'
#' @param x A list, the output of [derive_tdt_ttl()].
#' @param row An optional single whole number, the setting of `newdata` to draw, as
#'   its row number among the settings. `NULL` (the default) draws the first, with a
#'   message when there are several.
#' @param interval A single logical value, whether to add the lower and upper bounds
#'   of the credible interval as two faint surfaces around the median.
#' @param grid A single logical value, or one or two whole numbers. `TRUE` draws a
#'   grid on the median surface, ten lines at constant temperature and ten at constant
#'   elapsed time. A number sets how many lines each way, and two numbers the lines at
#'   constant temperature and at constant time, e.g. `c(5, 15)`. `FALSE` (the default)
#'   draws none.
#' @param time_unit An optional character string giving the time unit for the axis
#'   label, e.g. `"min"`.
#' @return A plotly object.
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
#'
#' tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
#' plot_ttl_3d(tol, row = 2, time_unit = "min")
#' }
#' @export
plot_ttl_3d <- function(x, row = NULL, interval = FALSE, grid = FALSE,
                        time_unit = NULL) {
  if (!.has_plotly())
    cli::cli_abort(c(
      "{.fn plot_ttl_3d} needs the {.pkg plotly} package.",
      "i" = "Install it with {.code install.packages(\"plotly\")}, or use {.fn plot_ttl} for the landscape as a map."
    ), class = "rTFT_needs_plotly")
  if (!(isTRUE(grid) || isFALSE(grid) ||
        (is.numeric(grid) && length(grid) %in% 1:2 && all(is.finite(grid)) &&
         all(grid >= 1) && all(grid == round(grid)))))
    cli::cli_abort(c(
      "{.arg grid} must be {.code TRUE}, {.code FALSE}, or one or two whole numbers of lines.",
      "x" = "You supplied {.val {grid}}.",
      "i" = "For example {.code grid = 20}, or {.code grid = c(5, 15)} for 5 lines at constant temperature and 15 at constant time."
    ), class = "rTFT_bad_grid")
  tbl  <- .plot_table(x, "ttl_surface")
  cols <- .newdata_cols(tbl)

  # One setting of newdata per figure: a surface cannot be drawn in panels.
  settings <- unique(tbl[cols])
  n <- max(nrow(settings), 1L)
  chosen <- if (is.null(row)) 1L else row
  if (!is.numeric(chosen) || length(chosen) != 1L || !chosen %in% seq_len(n))
    cli::cli_abort(c(
      "{.arg row} must be the number of one setting of {.arg newdata}, from 1 to {n}.",
      "x" = "You supplied {.val {row}}."
    ), class = "rTFT_bad_row")
  if (length(cols)) {
    this <- settings[chosen, , drop = FALSE]
    tbl  <- tbl[.curve_group(tbl, cols) %in% .curve_group(this, cols), , drop = FALSE]
    if (is.null(row) && n > 1)
      cli::cli_inform(c(
        "i" = "Drawing the first of {n} settings of {.arg newdata}: {paste(cols, format(unlist(this)), sep = ' = ', collapse = ', ')}.",
        "i" = "Use {.arg row} to draw another."
      ), class = "rTFT_plot_setting")
  }

  # plotly takes the heights as a matrix with a row per temperature and a column
  # per elapsed time.
  tau <- sort(unique(tbl$tau))
  Ta  <- sort(unique(tbl$Ta))
  tbl <- tbl[order(tbl$Ta, tbl$tau), , drop = FALSE]
  heights <- function(col) matrix(tbl[[col]], nrow = length(Ta), ncol = length(tau),
                                  byrow = TRUE)

  # plotly takes a colour scale as pairs of a position from 0 to 1 and a colour.
  rdylbu <- purrr::map2(seq(0, 1, length.out = length(.rdylbu)), .rdylbu,
                        \(at, colour) list(at, colour))
  g <- plotly::plot_ly(x = tau, y = Ta, z = heights("median"), type = "surface",
                       colorscale = rdylbu, cmin = 0, cmax = 1,
                       colorbar = list(title = list(text = "Survival", side = "top"),
                                       orientation = "h", y = 1, yanchor = "bottom",
                                       len = 0.5))
  if (isTRUE(interval))
    for (col in c("lower", "upper"))
      g <- plotly::add_surface(g, x = tau, y = Ta, z = heights(col), opacity = 0.3,
                               colorscale = rdylbu, cmin = 0, cmax = 1,
                               showscale = FALSE)
  # A grid on the median surface: lines at constant temperature and at constant
  # elapsed time, each following the surface. `grid` gives how many each way: ten for
  # TRUE, the number given, or two numbers (temperature lines, then time lines).
  # There cannot be more lines than the landscape has grid values.
  if (!isFALSE(grid)) {
    n_lines <- if (isTRUE(grid)) c(10, 10) else rep_len(grid, 2)
    z <- heights("median")
    spread <- function(n, k) unique(round(seq(1, n, length.out = min(n, k))))
    grid_line <- function(g, x, y, z) {
      plotly::add_trace(g, x = x, y = y, z = z, type = "scatter3d", mode = "lines",
                        line = list(color = "#333333", width = 2), showlegend = FALSE,
                        hoverinfo = "skip", inherit = FALSE)
    }
    for (i in spread(length(Ta), n_lines[1]))
      g <- grid_line(g, tau, rep(Ta[i], length(tau)), z[i, ])
    for (j in spread(length(tau), n_lines[2]))
      g <- grid_line(g, rep(tau[j], length(Ta)), Ta, z[, j])
  }
  time_axis <- list(title = .time_label("Elapsed exposure time", time_unit))
  if (all(tau > 0)) {
    ticks <- .log_ticks(range(tau))
    time_axis <- c(time_axis, list(type = "log", tickmode = "array", tickvals = ticks,
                                   ticktext = vapply(ticks, .plain_number, character(1))))
  } else {
    time_axis$type <- "linear"
  }
  plotly::layout(g, scene = list(
    xaxis = time_axis,
    yaxis = list(title = "Temperature (\u00b0C)"),
    zaxis = list(title = "Survival", range = c(0, 1))))
}
