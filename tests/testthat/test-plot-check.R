# ===========================================================================
# 1. Helpers
# ===========================================================================

geoms <- function(g) unname(vapply(g$layers, \(l) class(l$geom)[1], character(1)))

# Small tables in the layouts check_survival() and hazard_residuals() return,
# written by hand, so the plots are tested without a fitted model.
check_tbl <- function(path = c("r1", "r2")) {
  g <- expand.grid(time = c(0, 10, 20, 30), path_id = path, stringsAsFactors = FALSE)
  s <- exp(-g$time / 30)
  tibble::tibble(path_id = g$path_id, n = 10L, time = g$time, Ta = 20 + g$time / 4,
                 observed = round(s * 10) / 10, median = s, lower = s * 0.7,
                 upper = pmin(1, s * 1.3), prob = 0.95)
}
res_tbl <- function() {
  tibble::tibble(id = 1:8, path_id = "r1", mass = c(6, 9, 12, 15, 18, 21, 24, 27),
                 trial = rep(c("A", "B"), each = 4), time = c(5, 10, 10, 15, 20, 20, 25, 30),
                 Ta = 21:28, event = c(1, 1, 1, 1, 1, 1, 1, 0),
                 residual = c(0.9, 0.2, 0.7, 0.5, 0.35, 0.6, 0.1, 0.05))
}


# ===========================================================================
# 2. The posterior predictive check of survival
# ===========================================================================

test_that("plot_check draws the simulated band, its median and the observed steps", {
  g <- plot_check(check_tbl(), time_unit = "min")
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomRibbon", "GeomLine", "GeomStep"))
  expect_equal(g$labels$x, "Time (min)")
  expect_equal(g$labels$y, "Survival")
  expect_equal(plot_check(check_tbl())$labels$x, "Time")
  built <- ggplot2::ggplot_build(g)
  expect_equal(nrow(built$layout$layout), 2L)                  # a panel per trajectory
  expect_equal(built$data[[3]]$y[1:4], check_tbl()$observed[1:4])
  # with a further grouping, a panel for each group of each trajectory
  by_trial <- dplyr::bind_rows(dplyr::mutate(check_tbl(), trial = "A", .after = "path_id"),
                               dplyr::mutate(check_tbl(), trial = "B", .after = "path_id"))
  expect_equal(nrow(ggplot2::ggplot_build(plot_check(by_trial))$layout$layout), 4L)
  expect_error(plot_check(res_tbl()), class = "rTFT_bad_plot_input")
  expect_error(plot_check(1), class = "rTFT_bad_plot_input")
})


test_that("plot_check draws what check_survival returns", {
  fit <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))
  set.seed(1)
  expect_no_error(ggplot2::ggplot_build(plot_check(check_survival(fit, n_draws = 30))))
  built <- ggplot2::ggplot_build(plot_check(check_survival(fit, by = "trial", n_draws = 30)))
  expect_equal(nrow(built$layout$layout), 8L)
})


# ===========================================================================
# 3. Residuals by individual
# ===========================================================================

test_that("plot_residuals sets the residuals against the uniform distribution", {
  g <- plot_residuals(res_tbl())
  expect_s3_class(g, "ggplot")
  expect_equal(geoms(g), c("GeomAbline", "GeomPoint"))
  built <- ggplot2::ggplot_build(g)
  expect_equal(built$data[[2]]$y, sort(res_tbl()$residual))
  expect_equal(built$data[[2]]$x, (1:8 - 0.5) / 8)
  expect_equal(g$labels$x, "Expected under the model")
})


test_that("plot_residuals draws the residuals against a variable", {
  by_mass <- plot_residuals(res_tbl(), against = "mass")
  expect_equal(geoms(by_mass), c("GeomHline", "GeomPoint"))
  expect_equal(by_mass$labels$x, "mass")
  built <- ggplot2::ggplot_build(by_mass)
  expect_equal(built$data[[2]]$x, res_tbl()$mass)
  expect_equal(built$data[[1]]$yintercept, 0.5)
  by_trial <- plot_residuals(res_tbl(), against = "trial")
  expect_equal(geoms(by_trial), c("GeomHline", "GeomBoxplot", "GeomPoint"))
  expect_equal(nrow(ggplot2::ggplot_build(by_trial)$data[[2]]), 2L)   # a box per trial
})


test_that("plot_residuals refuses what it cannot draw", {
  expect_error(plot_residuals(check_tbl()), class = "rTFT_bad_plot_input")
  expect_error(plot_residuals(res_tbl(), against = "sex"), class = "rTFT_bad_against")
  expect_error(plot_residuals(res_tbl(), against = "residual"), class = "rTFT_bad_against")
  expect_error(plot_residuals(res_tbl(), against = c("mass", "trial")),
               class = "rTFT_bad_against")
})


test_that("plot_residuals draws what hazard_residuals returns", {
  fit <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))
  set.seed(1)
  res <- hazard_residuals(fit, n_draws = 50)
  expect_no_error(ggplot2::ggplot_build(plot_residuals(res)))
  expect_no_error(ggplot2::ggplot_build(plot_residuals(res, against = "trial")))
})
