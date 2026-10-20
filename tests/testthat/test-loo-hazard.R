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


test_that(".loo_units keeps its relative efficiency where a unit's likelihood underflows", {
  set.seed(1)
  steady  <- matrix(stats::rnorm(600 * 4, -2, 0.05), ncol = 4)
  chain   <- rep(1:2, each = 300)
  # one unit's log-likelihood shifted below -745, where exp() is zero
  shifted <- steady
  shifted[, 2] <- shifted[, 2] - 800
  expect_true(all(exp(shifted[, 2]) == 0))
  a <- .loo_units(steady, chain, "trial")
  b <- .loo_units(shifted, chain, "trial")
  # a constant shift of one unit changes its log predictive density by that constant
  # and nothing else: the Pareto k of every unit is the same
  expect_equal(b$diagnostics$pareto_k, a$diagnostics$pareto_k)
  expect_equal(b$pointwise[, "elpd_loo"] + c(0, 800, 0, 0), a$pointwise[, "elpd_loo"])
})


test_that(".loo_units asks loo for one core, whatever mc.cores says", {
  set.seed(1)
  steady <- matrix(stats::rnorm(600 * 4, -2, 0.05), ncol = 4)
  chain  <- rep(1:2, each = 300)
  # With more than one core loo copies the process with forks, which hung for good
  # on a cluster once a sampler's threads were running. The cores asked for in the
  # first call to each of loo's functions, the one rTFT makes, are recorded on the
  # way; loo's own later calls pass through.
  asked <- new.env()
  relative_eff <- loo::relative_eff
  psis_loo     <- loo::loo
  local_mocked_bindings(
    relative_eff = function(x, ..., cores = getOption("mc.cores", 1)) {
      if (is.null(asked$relative_eff)) asked$relative_eff <- cores
      relative_eff(x, ..., cores = cores)
    },
    loo = function(x, ..., cores = getOption("mc.cores", 1)) {
      if (is.null(asked$loo)) asked$loo <- cores
      psis_loo(x, ..., cores = cores)
    },
    .package = "loo")
  withr::local_options(mc.cores = 4)
  out <- .loo_units(steady, chain, "trial")
  expect_equal(asked$relative_eff, 1)
  expect_equal(asked$loo, 1)
  expect_s3_class(out, "loo")
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
  # With eight trials, leaving one out changes the fit too much for PSIS to be
  # trusted for all of them, and the function says so.
  expect_warning(
    expect_message(out <- loo_hazard(fit_trials), class = "rTFT_loo_unit"),
    class = "rTFT_loo_unreliable"
  )
  expect_equal(out$unit, "trial")
  expect_equal(nrow(out$pointwise), 8)
  # the same as PSIS-LOO on the log-likelihood summed by trial by hand, which gives
  # the warning of the loo package about the same trial
  ll    <- brms::log_lik(fit_trials)
  trial <- as.character(fit_trials$data$trial)
  by_hand <- sapply(unique(trial), \(t) rowSums(ll[, trial == t, drop = FALSE]))
  ref <- suppressWarnings(
    loo::loo(by_hand, r_eff = loo::relative_eff(exp(by_hand), chain_id = rep(1:2, each = 300)))
  )
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


# ===========================================================================
# 3. Leave-one-unit-out by refitting (kfold_hazard)
# ===========================================================================

# Refitting takes minutes, so these tests replace the refit with a stand-in that
# records what each refit would have been given and returns a known value per unit
# left out: minus the unit's position among all the units. The fold logic and the
# assembly are then checked exactly; the refit itself is the brms update() and
# log_lik() already relied on, checked by hand and in the vignettes.
stand_in_refit <- function(record) {
  function(fit, data, train, unit_values, ...) {
    record$calls[[length(record$calls) + 1]] <- list(
      train = train, held_out = unique(as.character(unit_values[!train])),
      args = list(...), from = fit$from)
    held <- unique(as.character(unit_values[!train]))
    # the "refit" is a stand-in object recording which fold made it
    list(elpd = stats::setNames(-match(held, unique(as.character(unit_values))), held),
         fit = list(from = length(record$calls)))
  }
}

test_that(".log_mean_exp is the log of the mean, even for very negative values", {
  x <- cbind(c(-1, -2, -3), c(-0.5, -0.5, -0.5))
  expect_equal(.log_mean_exp(x), log(colMeans(exp(x))))
  far <- cbind(c(-1000, -1001))                          # exp() of these is 0 in R
  expect_equal(.log_mean_exp(far), -1000 + log(mean(exp(c(0, -1)))))
})


test_that(".assign_folds gives each unit a fold, as evenly as possible", {
  units <- paste0("T", 1:8)
  expect_equal(unname(.assign_folds(units, 8)), 1:8)     # one per unit
  set.seed(3)
  f <- .assign_folds(units, 3)
  expect_equal(names(f), units)
  expect_equal(sort(as.vector(table(f))), c(2, 3, 3))   # as even as possible
  set.seed(3)
  expect_equal(.assign_folds(units, 3), f)               # the same with the same seed
})


test_that("kfold_hazard refits once per trial, leaving each trial out in turn", {
  record <- new.env()
  record$calls <- list()
  local_mocked_bindings(.refit_fold = stand_in_refit(record))
  out <- suppressMessages(kfold_hazard(fit_trials))
  expect_length(record$calls, 8)                         # one refit per trial
  held <- vapply(record$calls, \(cl) cl$held_out, character(1))
  expect_setequal(held, paste0("T", 1:8))                # each trial left out once
  # each refit is given every row except the left-out trial's
  trial <- as.character(fit_trials$data$trial)
  for (cl in record$calls) expect_equal(cl$train, trial != cl$held_out)
  # the first fold is refitted from the original fit, each later one from the refit
  # before it, which saves compiling the model again
  expect_null(record$calls[[1]]$from)
  expect_equal(vapply(record$calls[-1], \(cl) cl$from, numeric(1)), 1:7)
  # the result: one value per trial, their sum, and the standard error over trials
  expect_s3_class(out, c("kfold", "loo"))
  expect_equal(out$unit, "trial")
  expect_equal(unname(out$pointwise[, "elpd_kfold"]), -(1:8))
  expect_equal(unname(out$estimates["elpd_kfold", ]), c(-36, sqrt(8 * var(-(1:8)))))
})


test_that("a unit that no refit scored is missing from the result, not counted as zero", {
  # a stand-in refit that returns no score for trial T5
  local_mocked_bindings(.refit_fold = function(fit, data, train, unit_values, ...) {
    held <- setdiff(unique(as.character(unit_values[!train])), "T5")
    list(elpd = stats::setNames(rep(-1, length(held)), held), fit = fit)
  })
  out <- suppressMessages(kfold_hazard(fit_trials))
  scores <- out$pointwise[, "elpd_kfold"]
  expect_equal(names(scores)[is.na(scores)], "T5")
  expect_equal(unname(scores[names(scores) != "T5"]), rep(-1, 7))
  expect_true(is.na(out$estimates["elpd_kfold", "Estimate"]))
})


test_that("kfold_hazard with K folds leaves out every unit exactly once", {
  record <- new.env()
  record$calls <- list()
  local_mocked_bindings(.refit_fold = stand_in_refit(record))
  set.seed(5)
  out <- suppressMessages(kfold_hazard(fit_trials, K = 3))
  expect_length(record$calls, 3)
  held <- unlist(lapply(record$calls, \(cl) cl$held_out))
  expect_setequal(held, paste0("T", 1:8))
  expect_length(held, 8)                                 # none twice
  expect_equal(sort(as.vector(table(out$folds))), c(2, 3, 3))
  expect_equal(unname(out$pointwise[, "elpd_kfold"]), -(1:8))  # in the units' order
})


test_that("kfold_hazard uses loo_hazard's choice of unit, and passes options to the refits", {
  record <- new.env()
  record$calls <- list()
  local_mocked_bindings(.refit_fold = stand_in_refit(record))
  # each call gives several messages; the one checked is caught and the rest kept
  # out of the test output
  suppressMessages(
    expect_message(kfold_hazard(fit_trials, K = 2, chains = 1), class = "rTFT_loo_unit"))
  expect_equal(record$calls[[1]]$args, list(chains = 1))
  suppressMessages(
    expect_message(kfold_hazard(fit_trials, K = 2), class = "rTFT_kfold_progress"))
  expect_warning(suppressMessages(kfold_hazard(fit_trials, unit = "id", K = 2)),
                 class = "rTFT_loo_leak")
})


test_that("kfold_hazard needs a sensible number of folds", {
  for (bad in list(1, 9, 2.5, "four", c(2, 3)))
    expect_error(suppressMessages(kfold_hazard(fit_trials, K = bad)), class = "rTFT_bad_K")
})


test_that("two kfold_hazard results can be compared with loo_compare", {
  record <- new.env()
  record$calls <- list()
  local_mocked_bindings(.refit_fold = stand_in_refit(record))
  a <- suppressMessages(kfold_hazard(fit_trials, K = 2))
  b <- a
  b$pointwise[, "elpd_kfold"] <- b$pointwise[, "elpd_kfold"] - 0.5
  b$estimates["elpd_kfold", "Estimate"] <- sum(b$pointwise[, "elpd_kfold"])
  comparison <- loo::loo_compare(list(a = a, b = b))
  # The models are in a `model` column from loo 2.10.0 and in the row names of a
  # matrix before it.
  models <- if ("model" %in% colnames(comparison)) comparison$model else rownames(comparison)
  expect_equal(models, c("a", "b"))
  expect_equal(as.numeric(comparison[models == "b", "elpd_diff"]), -4)
})


test_that(".kfold_setup takes folds that are given, in the order of the model's units", {
  own <- suppressMessages(.kfold_setup(fit_trials, K = 4))
  expect_equal(own$K, 4)
  expect_named(own$folds, paste0("T", 1:8))
  given <- stats::setNames(rep(1:2, 4), paste0("T", 8:1))
  taken <- suppressMessages(.kfold_setup(fit_trials, folds = given))
  expect_equal(taken$K, 2)
  expect_equal(taken$folds, given[paste0("T", 1:8)])
  expect_equal(taken$units, paste0("T", 1:8))
  expect_error(suppressMessages(.kfold_setup(fit_trials, folds = given[-1])),
               class = "rTFT_bad_folds")
  # folds are run through by number, so none may be skipped and the first is 1
  skipped   <- stats::setNames(rep(c(1, 3), 4), paste0("T", 1:8))
  from_zero <- stats::setNames(rep(0:1, 4), paste0("T", 1:8))
  expect_error(suppressMessages(.kfold_setup(fit_trials, folds = skipped)),
               class = "rTFT_bad_folds")
  expect_error(suppressMessages(.kfold_setup(fit_trials, folds = from_zero)),
               class = "rTFT_bad_folds")
})


# ===========================================================================
# 4. A fitted model of static and ramping trials together
# ===========================================================================

# A small fitted model saved as a fixture (made by fixtures/make-hazard-fit-mixed.R):
# three static trials and two ramps of ten individuals each, with no group-level
# term, so the unit is the individual.
fit_mixed <- readRDS(test_path("fixtures", "hazard-fit-mixed.rds"))

test_that("loo_hazard leaves out whole individuals of either kind of assay", {
  expect_no_message(out <- loo_hazard(fit_mixed), class = "rTFT_loo_unit")
  expect_s3_class(out, "psis_loo")
  expect_equal(out$unit, "id")
  expect_equal(nrow(out$pointwise), 50)
  expect_true(all(is.finite(out$pointwise[, "elpd_loo"])))
  # an individual's log-likelihood is the sum of its rows, whether they end in a
  # failure between two checks, in a failure at a known time or in neither
  ll <- brms::log_lik(fit_mixed)
  id <- fit_mixed$rtft$id
  by_hand <- vapply(unique(id), \(i) rowSums(ll[, id == i, drop = FALSE]),
                    numeric(nrow(ll)))
  expect_equal(unname(.sum_by_unit(ll, id)), unname(by_hand))
})
