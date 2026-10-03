# ===========================================================================
# 1. Shared helpers
# ===========================================================================

# Small tables in the layout derive_tdt_ttl() returns, written by hand, so the plots
# are tested without a fitted model. `tdt_tbl()` has one curve per mass and p over
# three temperatures; `ttl_tbl()` has survival over three temperatures and the given
# elapsed times.
tdt_tbl <- function(mass = c(8, 12)) {
  g <- expand.grid(Ta = c(25, 30, 35), p = c(0.5, 0.9), mass = mass)
  tibble::tibble(mass = g$mass, p = g$p, Ta = g$Ta,
                 median = 100 / (g$Ta - 20) * g$p, lower = 50 / (g$Ta - 20) * g$p,
                 upper = 200 / (g$Ta - 20) * g$p, prob = 0.95, n_draws = 100L,
                 n_unreached = 0L)
}
ttl_tbl <- function(tau = c(1, 10, 100), Ta = c(25, 30, 35)) {
  g <- expand.grid(tau = tau, Ta = Ta, mass = 8)
  s <- exp(-g$tau * (g$Ta - 20) / 500)
  tibble::tibble(mass = g$mass, Ta = g$Ta, tau = g$tau, median = s, lower = s * 0.9,
                 upper = pmin(1, s * 1.1), prob = 0.95, n_draws = 100L, n_unreached = 0L)
}
geoms    <- function(g) unname(vapply(g$layers, \(l) class(l$geom)[1], character(1)))
y_trans  <- function(g) g$scales$get_scales("y")$get_transformation()$name
has_x_log <- function(g) {
  sc <- g$scales$get_scales("x")
  !is.null(sc) && sc$get_transformation()$name == "log-10"
}

test_that("the covariate columns are those derive_tdt_ttl() did not write", {
  expect_equal(.newdata_cols(tdt_tbl()), "mass")
  expect_equal(.newdata_cols(tdt_tbl()[, -1]), character(0))
})


test_that("a time axis label carries its unit only when one is given", {
  expect_equal(.time_label("Time to failure", NULL), "Time to failure")
  expect_equal(.time_label("Time to failure", "min"), "Time to failure (min)")
})


test_that("times on a log axis are labelled as plain numbers", {
  expect_equal(.plain_number(c(0.001, 0.1, 1, 10, 1000, 1e5)),
               c("0.001", "0.1", "1", "10", "1000", "100000"))
})


test_that("the palettes keep hard-to-see colours away from lines", {
  # lines run from cool blue to warm red, without a pale colour between
  heat <- .cool_warm(5)
  expect_equal(toupper(heat[c(1, 5)]), c("#2C7BB6", "#D7191C"))
  expect_true(all(colSums(grDevices::col2rgb(heat)) < 600))   # none near white (765)
  expect_equal(toupper(.cool_warm(1)), "#2C7BB6")             # a single curve is blue
  # survival runs from red, through pale yellow at one half, to blue
  expect_equal(.rdylbu[c(1, 6, 11)], c("#A50026", "#FFFFBF", "#313695"))
})


test_that(".nearest_quartiles picks values that are in the data", {
  expect_equal(.nearest_quartiles(seq(5, 15, by = 0.5)), c(5, 7.5, 10, 12.5, 15))
  # an uneven set: each quartile is replaced by the nearest value present
  v <- c(1, 2, 3, 10, 50, 51, 100)
  expect_true(all(.nearest_quartiles(v) %in% v))
  expect_equal(range(.nearest_quartiles(v)), c(1, 100))
  expect_equal(.nearest_quartiles(c(3, 3, 3)), 3)
})


test_that("up to six settings each get a panel, and more are cut to quartile values", {
  few <- tdt_tbl(mass = c(8, 10, 12))
  expect_silent(kept <- .choose_panels(few, "mass", TRUE))
  expect_equal(kept, few)

  many <- tdt_tbl(mass = seq(5, 15, by = 0.5))
  expect_message(kept <- .choose_panels(many, "mass", TRUE), class = "rTFT_plot_panels")
  expect_equal(sort(unique(kept$mass)), c(5, 7.5, 10, 12.5, 15))
  # asking for all, or for no panels, keeps every setting and says nothing
  expect_silent(all_kept <- .choose_panels(many, "mass", "all"))
  expect_equal(all_kept, many)
  expect_equal(.choose_panels(many, "mass", FALSE), many)
})


test_that("the panel limit is six settings: six all shown, seven cut", {
  expect_silent(kept <- .choose_panels(tdt_tbl(mass = 1:6), "mass", TRUE))
  expect_equal(sort(unique(kept$mass)), 1:6)
  expect_message(kept <- .choose_panels(tdt_tbl(mass = 1:7), "mass", TRUE),
                 class = "rTFT_plot_panels")
  expect_equal(sort(unique(kept$mass)), c(1, 2, 4, 5, 7))   # nearest its quartiles
})


test_that("with two covariates, a numeric one is cut and a text one kept whole", {
  g <- expand.grid(Ta = c(25, 30), p = 0.5, mass = seq(5, 15, by = 0.5),
                   sex = c("F", "M"), stringsAsFactors = FALSE)
  tbl <- tibble::tibble(mass = g$mass, sex = g$sex, p = g$p, Ta = g$Ta, median = 1,
                        lower = 1, upper = 1, prob = 0.95, n_draws = 10L,
                        n_unreached = 0L)
  expect_message(kept <- .choose_panels(tbl, c("mass", "sex"), TRUE),
                 class = "rTFT_plot_panels")
  expect_equal(sort(unique(kept$mass)), c(5, 7.5, 10, 12.5, 15))
  expect_equal(sort(unique(kept$sex)), c("F", "M"))
  expect_equal(nrow(unique(kept[c("mass", "sex")])), 5 * 2)  # every combination kept
})


test_that("a covariate that is not numeric keeps all its levels", {
  g <- expand.grid(Ta = c(25, 30), p = 0.5, species = letters[1:8], stringsAsFactors = FALSE)
  tbl <- tibble::tibble(species = g$species, p = g$p, Ta = g$Ta, median = 1, lower = 1,
                        upper = 1, prob = 0.95, n_draws = 10L, n_unreached = 0L)
  expect_silent(kept <- .choose_panels(tbl, "species", TRUE))
  expect_equal(kept, tbl)
})


test_that("facet must be TRUE, FALSE or \"all\"", {
  expect_error(.choose_panels(tdt_tbl(), "mass", "some"), class = "rTFT_bad_facet")
  expect_error(.choose_panels(tdt_tbl(), "mass", NA), class = "rTFT_bad_facet")
})


test_that(".cell_edges puts each edge halfway to the neighbouring grid value", {
  e <- .cell_edges(c(1, 2, 4))
  expect_equal(e$lo, c(0.5, 1.5, 3))
  expect_equal(e$hi, c(1.5, 3, 5))
  # on a log scale halfway is the geometric mean
  e <- .cell_edges(c(1, 10, 100), log = TRUE)
  expect_equal(e$lo, 10^c(-0.5, 0.5, 1.5))
  expect_equal(e$hi, 10^c(0.5, 1.5, 2.5))
  # repeated and unordered values each get the edges of their own grid value
  expect_equal(.cell_edges(c(4, 1, 4, 2))$lo, c(3, 0.5, 3, 1.5))
  expect_equal(.cell_edges(7), list(lo = 6.5, hi = 7.5))
})


# ===========================================================================
# 2. Thermal death time curves (plot_tdt)
# ===========================================================================

test_that("plot_tdt draws a ribbon and a line, with time on a log scale", {
  g <- plot_tdt(list(tdt_curves = tdt_tbl()), time_unit = "min")
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomRibbon", "GeomLine"))
  expect_equal(y_trans(g), "log-10")
  labs <- ggplot2::get_labs(g)
  expect_equal(labs$x, "Temperature (\u00b0C)")
  expect_equal(labs$y, "Time to failure (min)")
  expect_null(labs$title)                                # no title
  expect_equal(g$theme$legend.position, "top")           # the legend above the panels
  # the time axis is labelled 1, 10, 100, not 1e+00, 1e+01, 1e+02
  expect_equal(g$scales$get_scales("y")$labels(c(1, 10, 100)), c("1", "10", "100"))
  expect_no_error(ggplot2::ggplot_build(g))
})


test_that("plot_tdt colours by p and gives each setting of newdata a panel", {
  g <- plot_tdt(list(tdt_curves = tdt_tbl()))
  expect_s3_class(g$facet, "FacetWrap")
  expect_equal(names(g$facet$params$facets), "mass")
  expect_s3_class(g$data$p, "factor")                    # one colour per fraction
  # the fewest failed in blue, the most in red, as on the landscape
  line <- ggplot2::ggplot_build(g)$data[[2]]
  expect_equal(toupper(unique(line$colour[order(line$group)])), c("#2C7BB6", "#D7191C"))
  expect_equal(nlevels(g$data$.group), 4)                # 2 masses x 2 fractions
  # without panels the curves are still separate
  g2 <- plot_tdt(list(tdt_curves = tdt_tbl()), facet = FALSE)
  expect_s3_class(g2$facet, "FacetNull")
  expect_equal(nlevels(g2$data$.group), 4)
})


test_that("plot_tdt can leave out the interval", {
  g <- plot_tdt(list(tdt_curves = tdt_tbl()), interval = FALSE)
  expect_equal(geoms(g), "GeomLine")
})


test_that("a median that is never reached is left out of the line", {
  tbl <- tdt_tbl()
  tbl$median[1] <- Inf
  tbl$upper[1:2] <- Inf
  g <- plot_tdt(list(tdt_curves = tbl))
  expect_true(is.na(g$data$median[1]))
  expect_equal(g$data$upper[1:2], c(Inf, Inf))           # the interval stays unbounded
  expect_no_error(ggplot2::ggplot_build(g))
})


test_that("colour_by draws every value of a covariate together, a panel per p", {
  many <- tdt_tbl(mass = seq(5, 15, by = 0.5))
  expect_silent(g <- plot_tdt(list(tdt_curves = many), colour_by = "mass"))
  expect_equal(geoms(g), "GeomLine")                     # no interval by default
  expect_equal(length(unique(g$data$mass)), 21)          # every value kept
  expect_equal(names(g$facet$params$facets), "p")
  expect_equal(nlevels(g$data$.group), 21 * 2)
  # the lowest value of the covariate in teal, the highest in deep purple
  line <- ggplot2::ggplot_build(g)$data[[1]]
  expect_true(all(c("#7FCDBB", "#3F007D") %in% toupper(line$colour)))
  expect_equal(length(unique(line$colour)), 21)          # a colour per value
  expect_no_error(ggplot2::ggplot_build(g))
  # the interval can be asked for
  expect_equal(geoms(plot_tdt(list(tdt_curves = many), colour_by = "mass",
                              interval = TRUE)), c("GeomRibbon", "GeomLine"))
})


test_that("plot_tdt refuses a bad colour_by or an input that is not derive_tdt_ttl's", {
  x <- list(tdt_curves = tdt_tbl())
  expect_error(plot_tdt(x, colour_by = "size"), class = "rTFT_bad_colour_by")
  x$tdt_curves$sex <- "F"
  expect_error(plot_tdt(x, colour_by = "sex"), class = "rTFT_bad_colour_by")
  expect_error(plot_tdt(list(params = 1)), class = "rTFT_bad_plot_input")
  expect_error(plot_tdt(tdt_tbl()), class = "rTFT_bad_plot_input")
})


# ===========================================================================
# 3. Thermal tolerance landscape (plot_ttl)
# ===========================================================================

test_that("plot_ttl draws the landscape as cells with contour lines", {
  g <- plot_ttl(list(ttl_surface = ttl_tbl()), time_unit = "min")
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomRect", "GeomContour", "GeomText"))
  expect_equal(g$layers[[3]]$aes_params$size, 2.5)       # small labels, no box
  expect_true(has_x_log(g))
  labs <- ggplot2::get_labs(g)
  expect_equal(labs$x, "Elapsed exposure time (min)")
  expect_equal(labs$y, "Temperature (\u00b0C)")
  expect_null(labs$title)
  # the cells meet halfway between grid values, on the log scale for time
  expect_equal(sort(unique(g$data$tau_lo)), 10^c(-0.5, 0.5, 1.5))
  expect_equal(sort(unique(g$data$Ta_lo)), c(22.5, 27.5, 32.5))
  expect_equal(g$theme$legend.position, "top")
  expect_no_error(ggplot2::ggplot_build(g))
})


test_that("each contour line is labelled with its survival level, on the line", {
  tbl <- ttl_tbl(tau = c(1, 10, 100), Ta = c(25, 30, 35))
  at  <- .contour_labels(tbl, "mass", c(0.1, 0.5, 0.9), log_time = TRUE)
  expect_equal(at$label, c("0.1", "0.5", "0.9"))
  expect_equal(at$mass, rep(8, 3))
  # 0.1 and 0.9 are crossed at a grid time, where the temperature is interpolated:
  # interpolating that time's column back at the label's temperature gives the level
  for (i in c(1, 3)) {
    column <- tbl[tbl$tau == at$tau[i], ]
    expect_equal(stats::approx(column$Ta, column$median, xout = at$Ta[i])$y,
                 c(0.1, 0.5, 0.9)[i])
  }
  # 0.5 passes between the grid times 10 and 100 at every temperature, so it is
  # found along a grid temperature, the time interpolated on the log scale
  expect_true(at$Ta[2] %in% c(25, 30, 35))
  row <- tbl[tbl$Ta == at$Ta[2], ]
  expect_equal(stats::approx(log(row$tau), row$median, xout = log(at$tau[2]))$y, 0.5)
  # a level the landscape never reaches gets no label
  expect_equal(.contour_labels(tbl, "mass", 0.9999)$label, character(0))
  # and the labels can be left off
  expect_equal(geoms(plot_ttl(list(ttl_surface = tbl), labels = FALSE)),
               c("GeomRect", "GeomContour"))
})


test_that("the time axis is not on a log scale when the times include zero", {
  g <- plot_ttl(list(ttl_surface = ttl_tbl(tau = c(0, 10, 20))))
  expect_false(has_x_log(g))
  expect_no_error(ggplot2::ggplot_build(g))
})


test_that("a landscape with a single temperature has no contour lines", {
  g <- plot_ttl(list(ttl_surface = ttl_tbl(Ta = 30)))
  expect_equal(geoms(g), "GeomRect")
})


test_that("the map refuses several settings in one panel, and the curves allow it", {
  two <- dplyr::bind_rows(ttl_tbl(), dplyr::mutate(ttl_tbl(), mass = 12))
  x <- list(ttl_surface = two)
  expect_error(plot_ttl(x, facet = FALSE), class = "rTFT_surface_needs_panels")
  expect_no_error(plot_ttl(x))                           # a panel each is fine
  expect_no_error(plot_ttl(list(ttl_surface = ttl_tbl()), facet = FALSE))  # one setting
  g <- plot_ttl(x, type = "curves", facet = FALSE)
  expect_s3_class(g$facet, "FacetNull")
  expect_no_error(ggplot2::ggplot_build(g))
})


test_that("plot_ttl draws survival curves at chosen temperatures", {
  x <- list(ttl_surface = ttl_tbl(Ta = seq(20, 40, by = 1)))
  g <- plot_ttl(x, type = "curves")
  expect_equal(geoms(g), c("GeomRibbon", "GeomLine"))
  expect_equal(levels(g$data$Ta), c("20", "25", "30", "35", "40"))   # the quartiles
  expect_equal(ggplot2::get_labs(g)$y, "Survival")
  # asked-for temperatures are replaced by the nearest on the grid
  g2 <- plot_ttl(x, type = "curves", Ta = c(22.4, 33.6))
  expect_equal(levels(g2$data$Ta), c("22", "34"))
  expect_no_error(ggplot2::ggplot_build(g2))
})


test_that("plot_ttl refuses an input that is not derive_tdt_ttl's, or a bad type", {
  expect_error(plot_ttl(list(tdt_curves = tdt_tbl())), class = "rTFT_bad_plot_input")
  expect_error(plot_ttl(list(ttl_surface = ttl_tbl()), type = "lines"))
})


# ===========================================================================
# 4. From a fitted model
# ===========================================================================

test_that("the plots build from what derive_tdt_ttl() returns for a fitted model", {
  fit <- readRDS(test_path("fixtures", "hazard-fit-exp.rds"))
  tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60,
                        n_draws = 50)
  expect_no_error(ggplot2::ggplot_build(plot_tdt(tol)))
  expect_no_error(ggplot2::ggplot_build(plot_ttl(tol)))
  expect_no_error(ggplot2::ggplot_build(plot_ttl(tol, type = "curves")))
})


# ===========================================================================
# 5. Thermal tolerance landscape in three dimensions (plot_ttl_3d)
# ===========================================================================

# The built plotly figure holds its surfaces as traces: `x` the elapsed times, `y`
# the temperatures, and `z` the heights, a row per temperature.
surfaces <- function(g) plotly::plotly_build(g)$x$data

test_that("plot_ttl_3d draws the median survival as a surface over time and temperature", {
  tbl <- ttl_tbl(tau = c(1, 10, 100), Ta = c(25, 30))
  g <- plot_ttl_3d(list(ttl_surface = tbl), time_unit = "min")
  expect_s3_class(g, "plotly")
  s <- surfaces(g)
  expect_length(s, 1)
  expect_equal(s[[1]]$type, "surface")
  expect_equal(as.vector(s[[1]]$x), c(1, 10, 100))
  expect_equal(as.vector(s[[1]]$y), c(25, 30))
  # the height at [temperature i, time j] is that cell's median survival
  expect_equal(dim(s[[1]]$z), c(2L, 3L))
  expect_equal(s[[1]]$z[2, 3], tbl$median[tbl$Ta == 30 & tbl$tau == 100])
  expect_equal(s[[1]]$z[1, 2], tbl$median[tbl$Ta == 25 & tbl$tau == 10])
  scene <- plotly::plotly_build(g)$x$layout$scene
  expect_equal(scene$xaxis$type, "log")
  # round ticks written as plain numbers, not plotly's single digits
  expect_equal(as.vector(scene$xaxis$tickvals), c(1, 2, 5, 10, 20, 50, 100))
  expect_equal(as.vector(scene$xaxis$ticktext), c("1", "2", "5", "10", "20", "50", "100"))
  expect_equal(scene$xaxis$title, "Elapsed exposure time (min)")
  expect_equal(scene$zaxis$range, c(0, 1))
  # coloured from red at a survival of 0 to blue at 1
  expect_equal(s[[1]]$colorscale[[1]][[2]], "#A50026")
  expect_equal(s[[1]]$colorscale[[11]][[2]], "#313695")
})


test_that(".log_ticks gives round values within the range", {
  expect_equal(.log_ticks(c(0.3, 8)), c(0.5, 1, 2, 5))
  # fewer than four on the 1-2-5 pattern: 3 is added
  expect_equal(.log_ticks(c(2.5, 42.5)), c(3, 5, 10, 20, 30))
  # narrower still: every whole multiple
  expect_equal(.log_ticks(c(2.5, 8)), 3:8)
})


test_that("plot_ttl_3d can add the interval as two more surfaces", {
  tbl <- ttl_tbl()
  s <- surfaces(plot_ttl_3d(list(ttl_surface = tbl), interval = TRUE))
  expect_length(s, 3)
  expect_equal(s[[2]]$z[1, 1], tbl$lower[tbl$Ta == 25 & tbl$tau == 1])
  expect_equal(s[[3]]$z[1, 1], tbl$upper[tbl$Ta == 25 & tbl$tau == 1])
})


test_that("plot_ttl_3d can draw a grid of lines that follow the surface", {
  tbl <- ttl_tbl(tau = c(1, 10, 100), Ta = c(25, 30))
  s <- surfaces(plot_ttl_3d(list(ttl_surface = tbl), grid = TRUE))
  types <- vapply(s, \(trace) trace$type, character(1))
  expect_equal(types, c("surface", rep("scatter3d", 2 + 3)))   # 2 temperatures, 3 times
  # the first line is at the lowest temperature and follows the surface along time
  expect_equal(as.vector(s[[2]]$y), rep(25, 3))
  expect_equal(as.vector(s[[2]]$z), as.vector(s[[1]]$z[1, ]))
  # the last is at the longest time and follows the surface across temperature
  expect_equal(as.vector(s[[6]]$x), rep(100, 2))
  expect_equal(as.vector(s[[6]]$z), as.vector(s[[1]]$z[, 3]))
  # a large landscape gets ten lines each way, not one per grid value
  big <- list(ttl_surface = ttl_tbl(tau = 1:40, Ta = 21:60))
  expect_length(surfaces(plot_ttl_3d(big, grid = TRUE)), 1 + 10 + 10)
  # the number of lines can be chosen, the same each way or one for each
  expect_length(surfaces(plot_ttl_3d(big, grid = 20)), 1 + 20 + 20)
  s2 <- surfaces(plot_ttl_3d(big, grid = c(5, 15)))
  expect_length(s2, 1 + 5 + 15)
  expect_equal(length(unique(as.vector(s2[[2]]$y))), 1)       # a temperature line first
  # never more lines than grid values
  expect_length(surfaces(plot_ttl_3d(list(ttl_surface = tbl), grid = 50)), 1 + 2 + 3)
  for (bad in list(0, 2.5, "ten", c(1, 2, 3), NA))
    expect_error(plot_ttl_3d(big, grid = bad), class = "rTFT_bad_grid")
})


test_that("the 3D time axis is not on a log scale when the times include zero", {
  g <- plot_ttl_3d(list(ttl_surface = ttl_tbl(tau = c(0, 10, 20))))
  expect_equal(plotly::plotly_build(g)$x$layout$scene$xaxis$type, "linear")
})


test_that("plot_ttl_3d draws one setting of newdata, the first unless told", {
  two <- dplyr::bind_rows(ttl_tbl(), dplyr::mutate(ttl_tbl(), mass = 12, median = median / 2))
  x <- list(ttl_surface = two)
  expect_message(g1 <- plot_ttl_3d(x), class = "rTFT_plot_setting")
  expect_equal(surfaces(g1)[[1]]$z[1, 1], two$median[two$mass == 8][1])
  expect_silent(g2 <- plot_ttl_3d(x, row = 2))
  expect_equal(surfaces(g2)[[1]]$z[1, 1], two$median[two$mass == 12][1])
  expect_error(plot_ttl_3d(x, row = 3), class = "rTFT_bad_row")
  expect_error(plot_ttl_3d(x, row = "a"), class = "rTFT_bad_row")
})


test_that("plot_ttl_3d says so when plotly is not installed", {
  local_mocked_bindings(.has_plotly = function() FALSE)
  expect_error(plot_ttl_3d(list(ttl_surface = ttl_tbl())), class = "rTFT_needs_plotly")
})


test_that("the 3D landscape builds from a fitted model's output", {
  fit <- readRDS(test_path("fixtures", "hazard-fit-exp.rds"))
  tol <- derive_tdt_ttl(fit, newdata = data.frame(mass = 10), t_ref = 60, n_draws = 50)
  s <- surfaces(plot_ttl_3d(tol))
  expect_equal(dim(s[[1]]$z), c(50L, 50L))               # the default grids
  expect_true(all(s[[1]]$z >= 0 & s[[1]]$z <= 1))
})
