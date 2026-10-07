# _precompute.R - fit the models the vignettes report, once
#
# The vignettes show the code that fits each model but do not run it when they are
# built: a fit takes minutes, and a built package cannot depend on a sampler being
# installed. This script fits every model the vignettes report and saves, for each
# vignette, only what it displays (tables of estimates and the tables the plots are
# drawn from), with the package version and the date. Everything else in a vignette
# runs when it is built.
#
# Each fit is checked before anything is saved: rank-normalized split-Rhat below
# 1.01, bulk and tail effective sample sizes above 400, and no divergent
# transitions. A fit that fails stops the script, so a vignette cannot report a
# model that did not converge.
#
# Rerun this script from the package root after a change to the fitting or deriving
# functions, then rebuild the vignettes and read them against the new numbers.
#
# Inputs:
#   data/static_assay.rda
#   data/stream_record.rda
#
# Outputs:
#   vignettes/results-01.rds

# Setup ------------------------------------------------------------------------

library(posterior)
library(pkgload)

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from vignettes/", call. = FALSE)
}
load_all(quiet = TRUE)

stamp <- list(version = as.character(utils::packageVersion("rTFT")),
              date    = format(Sys.Date()))
seed  <- 20261006

# Convergence, checked on every fit --------------------------------------------

# The convergence summary of the population-level parameters and the standard
# deviations of the group-level effects, after checking that the fit converged.
check_fit <- function(fit) {
  draws <- posterior::as_draws_array(fit)
  keep  <- grep("^(b|sd)_", posterior::variables(draws), value = TRUE)
  conv  <- posterior::summarise_draws(posterior::subset_draws(draws, variable = keep),
                                      "median", "sd", "rhat", "ess_bulk", "ess_tail")
  nuts  <- brms::nuts_params(fit)
  n_divergent <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  stopifnot("a fit did not converge" = all(conv$rhat < 1.01),
            "a fit has too few effective draws" = all(conv$ess_bulk > 400, conv$ess_tail > 400),
            "a fit has divergent transitions" = n_divergent == 0)
  list(table = as.data.frame(conv), n_divergent = n_divergent,
       n_draws = posterior::ndraws(draws))
}

# Vignette 1: getting started, on the static assay ------------------------------

static_path <- data.frame(
  path_id = rep(c("T23", "T25", "T27"), each = 2),
  time    = rep(c(0, 7200), times = 3),
  temp    = rep(c(23, 25, 27), each = 2))
static_checks <- expand.grid(
  path_id = c("T23", "T25", "T27"),
  time    = seq(0, 7200, by = 60))

long_static <- build_hazard_data(static_assay, static_path, "interval",
                                 covariates = c("mass", "trial"),
                                 checks = static_checks)

fit_static <- fit_hazard(
  long_static, m_form = "exp", baseline = "weibull",
  dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
  # a smaller step size than the default, as an example of an option passed on to
  # the sampler
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr", refresh = 0)

# The same model with the memoryless baseline, for the comparison of the two
# predictions along the stream record.
fit_static_mem <- fit_hazard(
  long_static, m_form = "exp", baseline = "mem",
  dpar_formulas = list(loga = ~ 1 + I(mass - 15) + (1 | trial)),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = seed,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr", refresh = 0)
invisible(check_fit(fit_static_mem))

# The checks and the prediction draw random numbers, so the seed is set before each.
set.seed(seed)
check_static <- check_survival(fit_static)
set.seed(seed)
residuals_static <- hazard_residuals(fit_static)

results_01 <- list(
  stamp       = stamp,
  convergence = check_fit(fit_static),
  # the posterior probability that the hazard rises with exposure time (k > 1)
  p_shape_above_1 = mean(as.matrix(fit_static, variable = "b_logshape_Intercept") > 0),
  check       = check_static,
  residuals   = residuals_static,
  # thermal tolerance at two body masses, on the default grids
  tolerance   = suppressMessages(
    derive_tdt_ttl(fit_static, newdata = data.frame(mass = c(10, 20)), t_ref = 60)),
  # the CTmax a ramping trial would record, for three ramps from 14 degC
  dctmax      = suppressMessages(
    derive_dctmax(fit_static, rate = c(0.1, 0.2, 0.3), T_start = 14,
                  newdata = data.frame(mass = 15))),
  # survival of a 15 g fish along the stream record, with each day's share of the risk
  prediction  = suppressMessages(
    predict_survival(fit_static, stream_record, newdata = data.frame(mass = 15),
                     by = "day")),
  # the same prediction under the memoryless baseline: its mortality only
  mortality_mem = suppressMessages(
    predict_survival(fit_static_mem, stream_record,
                     newdata = data.frame(mass = 15)))$mortality)

saveRDS(results_01, "vignettes/results-01.rds", compress = "xz")
