# ===========================================================================
# 1. Helpers
# ===========================================================================

geoms <- function(g) unname(vapply(g$layers, \(l) class(l$geom)[1], character(1)))


# ===========================================================================
# 2. Survival along a record
# ===========================================================================

# A small table in the layout predict_survival() returns: survival falling over four
# times, one curve per mass.
surv_tbl <- function(mass = c(8, 12)) {
  g <- expand.grid(time = c(0, 10, 20, 30), mass = mass)
  s <- exp(-g$time * g$mass / 400)
  tibble::tibble(mass = g$mass, time = g$time, median = s, lower = s * 0.8,
                 upper = pmin(1, s * 1.2), prob = 0.95, n_draws = 100L, n_unreached = 0L)
}

test_that("plot_survival draws the median and its interval against time", {
  g <- plot_survival(list(survival = surv_tbl()), time_unit = "min")
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomRibbon", "GeomLine"))
  expect_equal(g$labels$x, "Time (min)")
  expect_equal(g$labels$y, "Survival")
  expect_equal(plot_survival(list(survival = surv_tbl()))$labels$x, "Time")
  built <- ggplot2::ggplot_build(g)
  expect_equal(nrow(unique(built$data[[2]]["PANEL"])), 2L)     # a panel per mass
  expect_equal(built$layout$panel_params[[1]]$y.range[1] <= 0, TRUE)
  expect_equal(geoms(plot_survival(list(survival = surv_tbl()), interval = FALSE)),
               "GeomLine")
  # no panels: both curves in one, still two lines
  one <- ggplot2::ggplot_build(plot_survival(list(survival = surv_tbl()), facet = FALSE))
  expect_equal(nrow(unique(one$data[[2]]["PANEL"])), 1L)
  expect_equal(length(unique(one$data[[2]]$group)), 2L)
})


test_that("plot_survival draws a table with no covariates, and colours by one when asked", {
  bare <- surv_tbl(mass = 8)[, -1]
  expect_no_error(ggplot2::ggplot_build(plot_survival(list(survival = bare))))
  many <- surv_tbl(mass = 5:15)
  expect_silent(g <- plot_survival(list(survival = many), colour_by = "mass"))
  expect_equal(geoms(g), "GeomLine")                           # no intervals by default
  built <- ggplot2::ggplot_build(g)
  expect_equal(length(unique(built$data[[1]]$group)), 11L)
  expect_equal(nrow(unique(built$data[[1]]["PANEL"])), 1L)
  # without colour_by, many settings are cut down to the quartile values
  expect_message(plot_survival(list(survival = many)), class = "rTFT_plot_panels")
})


test_that("plot_survival draws the temperature record above the survival curve", {
  rec <- tibble::tibble(time = c(0, 10, 20, 30), temp = c(22, 26, 25, 23))
  x <- list(survival = surv_tbl(), record = rec)
  g <- plot_survival(x, temperature = TRUE, time_unit = "min")
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomLine", "GeomBlank", "GeomRibbon", "GeomLine"))
  built <- ggplot2::ggplot_build(g)
  # two rows (temperature above survival) by two masses
  expect_equal(nrow(built$layout$layout), 4L)
  expect_equal(levels(built$layout$layout$.row), c("Temperature (\u00b0C)", "Survival"))
  expect_equal(as.integer(built$layout$layout$ROW[built$layout$layout$.row == "Survival"]),
               c(2L, 2L))
  # the record is drawn once in each column, at its own temperatures
  expect_equal(nrow(built$data[[1]]), 2L * nrow(rec))
  expect_equal(range(built$data[[1]]$y), c(22, 26))
  # the survival row runs from 0 to 1 while the temperature row keeps its own range
  survival_panel <- built$layout$layout$PANEL[built$layout$layout$.row == "Survival"][1]
  y_range <- built$layout$panel_params[[as.integer(survival_panel)]]$y.range
  expect_true(y_range[1] <= 0 && y_range[2] >= 1)
  expect_null(g$labels$y)
  expect_equal(g$labels$x, "Time (min)")
  # with no covariates, and with the curves coloured by one
  bare <- list(survival = surv_tbl(mass = 8)[, -1], record = rec)
  expect_equal(nrow(ggplot2::ggplot_build(plot_survival(bare, temperature = TRUE))$layout$layout), 2L)
  many <- list(survival = surv_tbl(mass = 5:15), record = rec)
  coloured <- ggplot2::ggplot_build(plot_survival(many, temperature = TRUE, colour_by = "mass"))
  expect_equal(nrow(coloured$layout$layout), 2L)
  expect_equal(nrow(coloured$data[[1]]), nrow(rec))            # one record line
  # a co-stressor the record follows gets a row of its own when it changes, and
  # none when it stays constant
  with_do <- list(survival = surv_tbl(),
                  record = cbind(rec, do = c(8, 6, 4, 6), ph = 7))
  three <- ggplot2::ggplot_build(plot_survival(with_do, temperature = TRUE))
  expect_equal(levels(three$layout$layout$.row),
               c("Temperature (\u00b0C)", "do", "Survival"))
  expect_equal(nrow(three$layout$layout), 6L)                  # three rows by two masses
  do_panel <- three$layout$layout$PANEL[three$layout$layout$.row == "do"][1]
  expect_equal(range(three$data[[1]]$y[three$data[[1]]$PANEL == do_panel]), c(4, 8))
  # asked for without a record to draw
  expect_error(plot_survival(list(survival = surv_tbl()), temperature = TRUE),
               class = "rTFT_bad_plot_input")
  expect_error(plot_survival(x, temperature = "yes"), class = "rTFT_bad_temperature")
})


test_that("without panels each setting has its own colour and legend entry", {
  g <- plot_survival(list(survival = surv_tbl()), facet = FALSE)
  built <- ggplot2::ggplot_build(g)
  expect_equal(geoms(g), c("GeomRibbon", "GeomLine"))
  # the ends of the covariate colours, teal and deep purple
  expect_equal(unique(built$data[[2]]$colour), c("#7FCDBB", "#3F007D"))
  expect_equal(unique(built$data[[1]]$fill), c("#7FCDBB", "#3F007D"))
  expect_equal(.setting_colours(4), .covariate_colours)
  expect_equal(.setting_labels(two_way <- data.frame(sex = c("F", "F", "M"), mass = c(8, 8, 12)),
                               c("sex", "mass")),
               factor(c("F, 8", "F, 8", "M, 12"), levels = c("F, 8", "M, 12")))
  expect_equal(g$scales$get_scales("colour")$name, "mass")
  expect_equal(levels(g$data$.setting), c("8", "12"))
  # two covariates: the legend names both and each entry gives both values
  two <- surv_tbl(); two <- cbind(sex = rep(c("F", "M"), each = 4), two)
  g2 <- plot_survival(list(survival = two), facet = FALSE)
  expect_equal(g2$scales$get_scales("colour")$name, "sex, mass")
  expect_equal(levels(g2$data$.setting), c("F, 8", "M, 12"))
  # with the record above, it is drawn once
  rec <- tibble::tibble(time = c(0, 10, 20, 30), temp = c(22, 26, 25, 23))
  together <- ggplot2::ggplot_build(
    plot_survival(list(survival = surv_tbl(), record = rec), temperature = TRUE, facet = FALSE))
  expect_equal(nrow(together$layout$layout), 2L)
  expect_equal(nrow(together$data[[1]]), nrow(rec))
  # a single setting keeps the one colour and has no legend to draw
  single <- plot_survival(list(survival = surv_tbl(mass = 8)), facet = FALSE)
  expect_null(single$scales$get_scales("colour"))
  # more settings than colours: said, and drawn in one colour
  expect_message(many <- plot_survival(list(survival = surv_tbl(mass = 1:9)), facet = FALSE),
                 class = "rTFT_plot_colours")
  expect_null(many$scales$get_scales("colour"))
  expect_no_message(plot_survival(list(survival = surv_tbl(mass = 1:8)), facet = FALSE))
})


test_that("plot_survival refuses what it cannot draw", {
  x <- list(survival = surv_tbl())
  expect_error(plot_survival(x, colour_by = "size"), class = "rTFT_bad_colour_by")
  expect_error(plot_survival(x, facet = "some"), class = "rTFT_bad_facet")
  expect_error(plot_survival(list(tdt_curves = surv_tbl())), "predict_survival",
               class = "rTFT_bad_plot_input")
})


test_that("plot_survival draws what predict_survival and average_survival return", {
  fits <- list(exp     = readRDS(test_path("fixtures", "hazard-fit-exp.rds")),
               zeroexp = readRDS(test_path("fixtures", "hazard-fit-zeroexp.rds")))
  rec  <- data.frame(time = seq(0, 30, by = 5), temp = c(22, 23, 25, 26, 25, 24, 23))
  pred <- predict_survival(fits$exp, rec, newdata = data.frame(mass = c(8, 12)), n_draws = 100)
  expect_no_error(ggplot2::ggplot_build(plot_survival(pred)))
  expect_no_error(ggplot2::ggplot_build(plot_survival(pred, temperature = TRUE)))
  w   <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0.5, 0.5))
  avg <- average_survival(fits, w, rec, newdata = data.frame(mass = 10), n_draws = 100)
  expect_no_error(ggplot2::ggplot_build(plot_survival(avg)))
})


# ===========================================================================
# 3. Each period's share of the risk
# ===========================================================================

# A prediction in the layout predict_survival() returns with `by`: two masses, three
# days, the shares of each mass adding to one.
risk_x <- function(mass = c(8, 12), day = 1:3) {
  g <- expand.grid(day = day, mass = mass)
  share <- rep(seq_along(day) / sum(seq_along(day)), length(mass))
  list(survival = surv_tbl(mass),
       risk = tibble::tibble(mass = g$mass, day = g$day, median = share,
                             lower = share * 0.8, upper = pmin(1, share * 1.2),
                             prob = 0.95, n_draws = 100L))
}

test_that("plot_risk draws a bar per period with its interval", {
  g <- plot_risk(risk_x())
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomBlank", "GeomCol", "GeomLinerange"))
  expect_equal(g$labels$x, "day")
  expect_equal(g$labels$y, "Share of the risk")
  built <- ggplot2::ggplot_build(g)
  expect_equal(nrow(built$layout$layout), 2L)                  # a panel per mass
  expect_equal(nrow(built$data[[2]]), 6L)                      # three days by two masses
  expect_equal(built$data[[2]]$y[1:3], (1:3) / 6)
  expect_equal(built$data[[3]]$ymax[1:3], (1:3) / 6 * 1.2)
  expect_equal(geoms(plot_risk(risk_x(), interval = FALSE)), c("GeomBlank", "GeomCol"))
})


test_that("plot_risk draws the settings side by side when there are no panels", {
  g <- plot_risk(risk_x(), facet = FALSE)
  built <- ggplot2::ggplot_build(g)
  expect_equal(nrow(built$layout$layout), 1L)
  expect_equal(unique(built$data[[2]]$fill), c("#7FCDBB", "#3F007D"))
  expect_equal(g$scales$get_scales("fill")$name, "mass")
  # the two masses of a day sit either side of it, not on top of each other
  day_1 <- built$data[[2]][round(built$data[[2]]$x) == 1, ]
  expect_equal(nrow(day_1), 2L)
  expect_true(day_1$xmax[1] <= day_1$xmin[2] + 1e-8)
  expect_error(plot_risk(risk_x(mass = 1:9), facet = FALSE), class = "rTFT_plot_colours")
})


test_that("periods that are not numbers keep their order along the record", {
  x <- risk_x(day = c("Mon", "Tue", "Wed"))
  x$risk$day <- as.character(x$risk$day)
  built <- ggplot2::ggplot_build(plot_risk(x))
  expect_equal(levels(plot_risk(x)$data$day), c("Mon", "Tue", "Wed"))
  expect_equal(nrow(built$data[[2]]), 6L)
  # one setting, no covariates
  bare <- list(survival = surv_tbl(mass = 8)[, -1], risk = risk_x(mass = 8)$risk[, -1])
  expect_equal(nrow(ggplot2::ggplot_build(plot_risk(bare))$layout$layout), 1L)
})


test_that("a period with no risk in any draw has no bar", {
  x <- risk_x()
  x$risk[c("median", "lower", "upper")] <- NA_real_
  x$risk$n_draws <- 0L
  expect_no_warning(built <- ggplot2::ggplot_build(plot_risk(x)))
  expect_equal(nrow(built$data[[2]]), 0L)                      # no bars
  expect_equal(nrow(built$layout$layout), 2L)                  # the panels are still there
  # one period without risk among others: its bar alone is missing
  some <- risk_x()
  some$risk$median[some$risk$day == 2] <- NA_real_
  bars <- ggplot2::ggplot_build(plot_risk(some))$data[[2]]
  expect_equal(nrow(bars), 4L)
  expect_equal(sort(unique(bars$x)), c(1, 3))
})


test_that("plot_risk refuses what it cannot draw", {
  expect_error(plot_risk(list(survival = surv_tbl())), class = "rTFT_no_risk")
  expect_error(plot_risk(list(risk = 1)), class = "rTFT_bad_plot_input")
  expect_error(plot_risk(risk_x(), facet = "some"), class = "rTFT_bad_facet")
})


test_that("plot_risk draws what predict_survival and average_survival return", {
  fits <- list(exp     = readRDS(test_path("fixtures", "hazard-fit-exp.rds")),
               zeroexp = readRDS(test_path("fixtures", "hazard-fit-zeroexp.rds")))
  rec  <- data.frame(time = seq(0, 30, by = 5), temp = c(22, 23, 25, 26, 25, 24, 23),
                     day = c(1, 1, 1, 2, 2, 3, 3))
  pred <- predict_survival(fits$exp, rec, newdata = data.frame(mass = c(8, 12)),
                           by = "day", n_draws = 100)
  built <- ggplot2::ggplot_build(plot_risk(pred))
  expect_equal(nrow(built$data[[2]]), 6L)
  expect_equal(sum(built$data[[2]]$y[built$data[[2]]$PANEL == 1]), 1, tolerance = 0.03)
  w   <- tibble::tibble(model = c("exp", "zeroexp"), weight = c(0.5, 0.5))
  avg <- average_survival(fits, w, rec, newdata = data.frame(mass = 10), by = "day",
                          n_draws = 100)
  expect_no_error(ggplot2::ggplot_build(plot_risk(avg, facet = FALSE)))
})
