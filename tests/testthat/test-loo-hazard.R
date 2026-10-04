# ===========================================================================
# 1. Helpers
# ===========================================================================

# Two fitted-model fixtures: `fit_exp` has no group-level term (40 individuals), and
# `fit_trials` has a random intercept per trial on log a, shared by the six
# individuals of each of eight trials.
fit_exp    <- readRDS(test_path("fixtures", "hazard-fit-exp.rds"))
fit_trials <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))

# A hand-made design for .choose_unit(): six individuals with two rows each, in two
# trials of three, and four tanks nested in the trials (two of them holding two
# individuals, two holding one).
design <- data.frame(id    = rep(1:6, each = 2),
                     trial = rep(c("A", "B"), each = 6),
                     tank  = rep(c("a1", "a1", "a2", "b1", "b1", "b2"), each = 2))

test_that(".grouping_terms reads the groupings of the group-level terms", {
  expect_equal(.grouping_terms(fit_trials), "trial")
  expect_length(.grouping_terms(fit_exp), 0)
  # the syntaxes brms accepts, read as brms reads them, on stand-in fits
  terms_of <- function(...) .grouping_terms(list(formula = list(pforms = list(...))))
  expect_equal(terms_of(loga ~ 1 + (1 || trial)), "trial")
  expect_equal(terms_of(loga ~ 1 + (1 | p | trial), logb ~ 1 + (1 | p | trial)), "trial")
  expect_equal(terms_of(loga ~ 1 + (1 | trial/tank)), c("trial", "trial:tank"))
  expect_equal(terms_of(loga ~ 1 + (1 | trial:tank)), "trial:tank")
  expect_equal(terms_of(loga ~ 1 + (1 | gr(trial, by = sex))), "trial")   # by groups nothing
})


test_that("a grouping of several columns takes its value from their combination", {
  d <- data.frame(trial = c("A", "A", "B", "B"), tank = c("1", "2", "1", "2"))
  expect_equal(.grouping_values(d, "trial"), d$trial)
  expect_equal(.grouping_values(d, "trial:tank"), c("A:1", "A:2", "B:1", "B:2"))
})


test_that("tanks nested in trials are recognised even when tank names are reused", {
  # each trial's tanks are named "1" and "2": as columns, tank and trial look
  # crossed, but the grouping trial:tank is nested in trial
  reused <- data.frame(id    = rep(1:8, each = 2),
                       trial = rep(c("A", "B"), each = 8),
                       tank  = rep(c("1", "1", "2", "2"), each = 2, times = 2))
  expect_message(u <- .choose_unit(reused, reused$id, c("trial", "trial:tank")),
                 class = "rTFT_loo_unit")
  expect_equal(u$name, "trial")
  # leaving out tanks is allowed, and leaks the trial's effect
  expect_warning(u <- .choose_unit(reused, reused$id, c("trial", "trial:tank"),
                                   unit = "trial:tank"), class = "rTFT_loo_leak")
  expect_equal(length(unique(u$values)), 4)              # four tanks, not two
})


test_that(".sum_by_unit sums each unit's rows, draw by draw", {
  ll <- matrix(1:15, nrow = 3)                           # 3 draws, 5 rows
  out <- .sum_by_unit(ll, c("a", "a", "b", "c", "b"))
  expect_equal(colnames(out), c("a", "b", "c"))          # in order of appearance
  expect_equal(out[, "a"], ll[, 1] + ll[, 2])
  expect_equal(out[, "b"], ll[, 3] + ll[, 5])
  expect_equal(out[, "c"], ll[, 4])
})


test_that("with no shared group the unit is the individual, without a message", {
  expect_silent(u <- .choose_unit(design, design$id, character(0)))
  expect_equal(u$name, "id")
  expect_equal(u$values, design$id)
  # a group with one individual per level is not shared
  expect_silent(u <- .choose_unit(design, design$id, "id"))
  expect_equal(u$name, "id")
})


test_that("individuals sharing a group make the coarsest such group the unit", {
  expect_message(u <- .choose_unit(design, design$id, c("tank", "trial")),
                 class = "rTFT_loo_unit")
  expect_equal(u$name, "trial")                          # coarser than tank
  expect_equal(u$values, design$trial)
  expect_message(u <- .choose_unit(design, design$id, "tank"), class = "rTFT_loo_unit")
  expect_equal(u$name, "tank")
})


test_that("a unit finer than a shared group is allowed, with a warning", {
  expect_warning(u <- .choose_unit(design, design$id, "trial", unit = "id"),
                 class = "rTFT_loo_leak")
  expect_equal(u$name, "id")
  expect_equal(u$label, "individual")
  expect_warning(.choose_unit(design, design$id, c("tank", "trial"), unit = "tank"),
                 class = "rTFT_loo_leak")          # tank is finer than trial
  # a coarser unit takes whole tanks with it, so leaves nothing behind
  expect_silent(.choose_unit(design, design$id, c("tank", "trial"), unit = "trial"))
  # with both tank and trial leaking, the fix suggested is the coarsest, the trial,
  # which takes the tanks with it
  w <- tryCatch(.choose_unit(design, design$id, c("tank", "trial"), unit = "id"),
                warning = conditionMessage)
  expect_match(w, 'unit = "trial"', fixed = TRUE)
})


test_that("crossed groupings have no default unit, and either can be chosen", {
  # eight individuals, two rows each, in two trials crossed with two families
  crossed <- data.frame(id     = rep(1:8, each = 2),
                        trial  = rep(c("A", "B"), each = 8),
                        family = rep(c("f1", "f1", "f2", "f2"), each = 2, times = 2))
  expect_error(.choose_unit(crossed, crossed$id, c("trial", "family")),
               class = "rTFT_crossed_groups")
  # leaving out trials keeps the family effects, and leaving out families the trial
  # effects; each is allowed, with the warning naming what stays in the fit
  w <- tryCatch(.choose_unit(crossed, crossed$id, c("trial", "family"), unit = "trial"),
                warning = conditionMessage)
  expect_match(w, "family", fixed = TRUE)
  expect_warning(u <- .choose_unit(crossed, crossed$id, c("trial", "family"),
                                   unit = "family"), class = "rTFT_loo_leak")
  expect_equal(u$values, crossed$family)
})


test_that("a unit must be a column, and must keep each individual in one level", {
  expect_error(.choose_unit(design, design$id, "trial", unit = "pond"),
               class = "rTFT_bad_unit")
  expect_error(.choose_unit(design, design$id, "trial", unit = c("trial", "tank")),
               class = "rTFT_bad_unit")
  split_up <- design
  split_up$day <- rep(1:2, 6)                            # each individual on two days
  expect_error(.choose_unit(split_up, split_up$id, c("trial", "day"), unit = "day"),
               class = "rTFT_bad_unit")
  # a column that is not a grouping variable of the model, even one constant within
  # each individual, is not a unit
  with_mass <- design
  with_mass$mass <- rep(c(8, 8, 10, 10, 12, 12), each = 2)
  expect_error(.choose_unit(with_mass, with_mass$id, "trial", unit = "mass"),
               class = "rTFT_bad_unit")
})


test_that("the Pareto k threshold follows the number of draws", {
  expect_equal(.k_threshold(600), 1 - 1 / log10(600))   # about 0.64
  expect_equal(.k_threshold(1e5), 0.7)                   # capped
})


test_that(".loo_units warns when PSIS cannot be trusted for some units", {
  set.seed(1)
  steady  <- matrix(stats::rnorm(600 * 4, -2, 0.05), ncol = 4)
  chain   <- rep(1:2, each = 300)
  expect_silent(out <- .loo_units(steady, chain, "trial"))
  expect_equal(out$unit, "trial")
  # one unit whose log-likelihood swings widely from draw to draw: a few draws
  # dominate its importance weights
  wild <- steady
  wild[, 2] <- stats::rnorm(600, -10, 6)
  expect_warning(.loo_units(wild, chain, "trial"), class = "rTFT_loo_unreliable")
})


# ===========================================================================
# 2. Leave-one-unit-out for a fitted model (loo_hazard)
# ===========================================================================

test_that("loo_hazard leaves out individuals when they share no group", {
  # none of the package's own messages or warnings: the individual is the honest
  # unit here
  expect_no_message(out <- loo_hazard(fit_exp), class = "rTFT_loo_unit")
  expect_no_warning(loo_hazard(fit_exp), class = "rTFT_loo_leak")
  expect_s3_class(out, "psis_loo")
  expect_equal(out$unit, "id")
  expect_equal(nrow(out$pointwise), 40)
})


test_that("loo_hazard leaves out whole trials when individuals share them", {
  expect_message(out <- loo_hazard(fit_trials), class = "rTFT_loo_unit")
  expect_equal(out$unit, "trial")
  expect_equal(nrow(out$pointwise), 8)
  # the same as PSIS-LOO on the log-likelihood summed by trial by hand
  ll    <- brms::log_lik(fit_trials)
  trial <- as.character(fit_trials$data$trial)
  by_hand <- sapply(unique(trial), \(t) rowSums(ll[, trial == t, drop = FALSE]))
  ref <- loo::loo(by_hand, r_eff = loo::relative_eff(exp(by_hand), chain_id = rep(1:2, each = 300)))
  expect_equal(out$estimates, ref$estimates)
})


test_that("loo_hazard warns that leaving out individuals is optimistic within trials", {
  expect_warning(out <- loo_hazard(fit_trials, unit = "id"), class = "rTFT_loo_leak")
  expect_equal(nrow(out$pointwise), 48)
})


test_that("loo_hazard needs a fit with its settings and a valid unit", {
  bare <- fit_exp
  bare$rtft <- NULL
  expect_error(loo_hazard(bare), class = "rTFT_no_settings")
  # a record of individuals that no longer lines up with the fitted rows
  short <- fit_exp
  short$rtft$id <- short$rtft$id[-1]
  expect_error(loo_hazard(short), class = "rTFT_rows_mismatch")
  expect_error(loo_hazard(fit_trials, unit = "tank"), class = "rTFT_bad_unit")
  expect_error(loo_hazard(fit_exp, unit = "mass"), class = "rTFT_bad_unit")   # a covariate
})
