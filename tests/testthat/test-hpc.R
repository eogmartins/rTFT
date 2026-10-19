# ===========================================================================
# 1. Helpers
# ===========================================================================

# A small ramping assay in two trials of ten individuals, and two models of it with
# a random intercept per trial. Writing a bundle fits nothing, so the table and the
# models are all these tests need.
ramp   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
checks <- data.frame(path_id = "A", time = seq(0, 60, by = 5))
ind    <- data.frame(id = 1:20, path_id = "A", event = 1,
                     trial = rep(c("a", "b"), each = 10),
                     t_upper = rep(seq(15, 60, 5), 2))
long   <- suppressMessages(build_hazard_data(ind, ramp, "interval", covariates = "trial",
                                             checks = checks, grid_dT = Inf))
models <- build_hazard_models("exp", c("mem", "weibull"),
                              dpar_formulas = list(loga = ~ 1 + (1 | trial)))

# Where the cluster installs rTFT from, given to every bundle written here: the
# rTFT these tests load is a development version, whose source cannot be told.
src <- "eogmartins/rTFT@dev-0.2.0"

# Write a bundle into `dir`, or a new temporary folder, without its messages, and
# return the folder.
write_quietly <- function(..., dir = NULL) {
  if (is.null(dir))
    dir <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "run")
  suppressMessages(write_hpc_bundle(models, long, dir, account = "def-test",
                                    rtft_source = src, ...))
  dir
}

# The scripts a bundle holds.
scripts <- c("settings.sh", "install.sh", "install.R", "check.sh", "README.txt")

# The files a bundle holds, below its folder.
bundle_files <- function(dir) sort(list.files(dir, recursive = TRUE))

# Whether this machine's bash can read the files R writes. The bash that Windows
# offers through its Linux subsystem, for one, cannot read Windows paths; the
# scripts are tested where their users run them, on Linux, and on macOS.
bash_reads_files <- function() {
  if (!nzchar(Sys.which("bash"))) return(FALSE)
  probe <- withr::local_tempfile(fileext = ".sh")
  writeLines("true", probe)
  status <- suppressWarnings(system2("bash", c("-n", shQuote(probe)),
                                     stdout = FALSE, stderr = FALSE))
  identical(status, 0L)
}

# The cluster's software, as check_hpc_bundle() finds it, made present for the tests
# that are not about it: these machines need not have CmdStan.
with_cmdstan <- function(env = parent.frame()) {
  local_mocked_bindings(.has_cmdstanr = function() TRUE, .cmdstan_found = function() TRUE,
                        .env = env)
}


# A fitted model of eight trials of six individuals, with a random intercept per
# trial on log a (two chains, 600 draws), the table it was fitted to, and two models
# of that table. The tasks of a bundle of them are run with stand-ins for the fit and
# for the refits, as in the tests of fit_hazard_models(), since fitting takes minutes.
fit_trials   <- readRDS(test_path("fixtures", "hazard-fit-trials.rds"))
long_trials  <- cbind(fit_trials$data, id = fit_trials$rtft$id)
trial_models <- build_hazard_models("exp", c("mem", "weibull"),
                                    dpar_formulas = list(loga = ~ 1 + (1 | trial)))

new_record <- function() {
  record <- new.env()
  record$calls <- list()
  record$fail  <- FALSE
  record
}

# The stand-in for fit_hazard(): it records what it was given and returns the stored
# fit, or stops, as a sampler can, when `record$fail` is TRUE.
stand_in_fit <- function(record) {
  function(long, m_form, baseline, dpar_formulas, observation, costressors, ...) {
    record$calls[[length(record$calls) + 1]] <- list(
      m_form = m_form, baseline = baseline, args = list(...))
    if (record$fail) stop("the sampler stopped")
    fit_trials
  }
}

# The stand-in for a refit: it records what it was given and returns, for each unit
# left out, minus the unit's position among all the units.
stand_in_refit <- function(record) {
  function(fit, data, train, unit_values, ...) {
    held <- unique(as.character(unit_values[!train]))
    record$calls[[length(record$calls) + 1]] <- list(held = held, args = list(...))
    list(elpd = stats::setNames(-match(held, unique(as.character(unit_values))), held),
         fit = fit)
  }
}

# Write a bundle of the two models of the trials, with refits in three folds (six
# tasks), into a new temporary folder, and return the folder. Two bundles written
# alike are the same bundle, so `seed` is changed to write another; `models` gives
# other models of the trials.
trial_bundle <- function(kfold = TRUE, seed = 1, models = trial_models,
                         env = parent.frame()) {
  dir <- file.path(withr::local_tempdir(.local_envir = env), "run")
  set.seed(3)
  suppressMessages(write_hpc_bundle(
    models, long_trials, dir, account = "def-test", rtft_source = src,
    kfold = kfold, K = if (kfold) 3, chains = 2, seed = seed,
    sd_prior = "exponential(1)"))
  dir
}

# Run a task, or every task of a bundle, without the messages of the tasks.
run_quietly <- function(...) suppressMessages(run_hpc_task(...))
run_all <- function(dir, kfold = 1:6) {
  for (i in 1:2) run_quietly(dir, "fit", i)
  for (i in kfold) run_quietly(dir, "kfold", i)
}


# ===========================================================================
# 2. Writing the bundle (write_hpc_bundle)
# ===========================================================================

test_that("write_hpc_bundle writes the scripts, the inputs and the empty folders", {
  dir <- write_quietly(chains = 2, seed = 1)
  expect_equal(bundle_files(dir),
               c("README.txt", "bundle/checksums.txt", "bundle/data.rds",
                 "bundle/models.rds", "bundle/settings.rds", "bundle/tasks-fit.csv",
                 "check.sh", "install.R", "install.sh", "settings.sh"))
  expect_true(all(dir.exists(file.path(dir, c("stan", "results", "logs")))))
  # without cross-validation by refitting there is no table of refits
  expect_false(file.exists(file.path(dir, "bundle", "tasks-kfold.csv")))
})


test_that("the bundle holds the data, the models and the settings of the run", {
  dir <- write_quietly(chains = 2, iter = 1000, seed = 1, time = "02:30:00",
                       mem_per_cpu = "2000M", email = "me@example.org")
  expect_identical(readRDS(file.path(dir, "bundle", "data.rds")), long)
  # the formulas come back with environments of their own, so only their content is
  # compared
  expect_equal(readRDS(file.path(dir, "bundle", "models.rds")), models,
               ignore_formula_env = TRUE)
  s <- readRDS(file.path(dir, "bundle", "settings.rds"))
  expect_equal(s$rtft_version, as.character(utils::packageVersion("rTFT")))
  expect_equal(s$brms_version, as.character(utils::packageVersion("brms")))
  expect_equal(s$observation, "interval")
  # the unit is settled now, as fit_hazard_models() settles it
  expect_equal(s$unit, "trial")
  expect_equal(s$units, c("a", "b"))
  expect_false(s$kfold)
  expect_null(s$K)
  expect_null(s$folds)
  expect_equal(s$rule, list(rhat = 1.01, ess = 400, divergent = 0.001))
  expect_equal(s$chains, 2)
  # the arguments of the fits, with the backend the cluster uses
  expect_equal(s$dots, list(chains = 2, iter = 1000, seed = 1, backend = "cmdstanr"))
  expect_equal(s$slurm, list(account = "def-test", time = "02:30:00",
                             mem_per_cpu = "2000M", cpus_per_task = 2,
                             email = "me@example.org",
                             modules = c("StdEnv/2023", "r/4.5.0"), r_library = NULL))
  tasks <- utils::read.csv(file.path(dir, "bundle", "tasks-fit.csv"))
  expect_equal(tasks, data.frame(task = 1:2, model = c("exp_mem", "exp_weibull"),
                                 m_form = "exp", baseline = c("mem", "weibull")))
})


test_that("the checksums of the bundle are those of its files", {
  dir <- write_quietly(kfold = TRUE)
  sums <- utils::read.table(file.path(dir, "bundle", "checksums.txt"),
                            col.names = c("md5", "file"), colClasses = "character")
  expect_equal(sums$file, paste0("bundle/", c("data.rds", "models.rds", "settings.rds",
                                              "tasks-fit.csv", "tasks-kfold.csv")))
  expect_equal(sums$md5, unname(tools::md5sum(file.path(dir, sums$file))))
})


test_that("with kfold = TRUE the folds are assigned once and each refit is a task", {
  set.seed(4)
  dir <- write_quietly(kfold = TRUE)
  set.seed(4)
  expected <- .assign_folds(c("a", "b"), 2)
  s <- readRDS(file.path(dir, "bundle", "settings.rds"))
  expect_true(s$kfold)
  expect_equal(s$K, 2)
  expect_equal(s$folds, expected)
  tasks <- utils::read.csv(file.path(dir, "bundle", "tasks-kfold.csv"))
  expect_equal(tasks$task, 1:4)
  expect_equal(tasks$model, rep(c("exp_mem", "exp_weibull"), each = 2))
  expect_equal(tasks$fold, rep(1:2, 2))
  expect_equal(tasks$units, rep(names(expected)[order(expected)], 2))
})


test_that("the scripts are filled in from the template, with nothing left to fill", {
  dir <- write_quietly(chains = 3, time = "1-00:00:00", modules = c("StdEnv/2023", "r/4.4.0"),
                       r_library = "~/R/rtft",
                       remote = "someone@nibi.alliancecan.ca:scratch/run")
  settings <- readLines(file.path(dir, "settings.sh"))
  expect_true(all(c('ACCOUNT="def-test"', 'TIME="1-00:00:00"', 'MEM_PER_CPU="4G"',
                    'CPUS_PER_TASK="3"', 'MODULES="StdEnv/2023 r/4.4.0"',
                    'R_LIBRARY="~/R/rtft"', 'EMAIL=""') %in% settings))
  readme <- paste(readLines(file.path(dir, "README.txt")), collapse = "\n")
  expect_match(readme, "ssh someone@nibi.alliancecan.ca", fixed = TRUE)
  expect_match(readme, "cd scratch/run", fixed = TRUE)
  expect_match(readme, "someone@nibi.alliancecan.ca:scratch/run/", fixed = TRUE)
  # the copy makes the whole path of the folder on the cluster first
  expect_match(readme, 'rsync -av --rsync-path="mkdir -p scratch/run && rsync"', fixed = TRUE)
  for (file in scripts)
    expect_false(any(grepl("{{", readLines(file.path(dir, file)), fixed = TRUE)))
  # the version to install, and where from
  install <- readLines(file.path(dir, "install.R"))
  expect_true(any(grepl(paste0('^wanted +<- "', utils::packageVersion("rTFT"), '"$'),
                        install)))
  expect_true(any(grepl('remotes::install_github("eogmartins/rTFT@dev-0.2.0"', install,
                        fixed = TRUE)))
  # what the bundle holds, said with its plurals
  expect_match(readme, "It holds 2 models to fit, as 2 tasks, without cross-validation by",
               fixed = TRUE)
  with_refits <- readLines(file.path(write_quietly(kfold = TRUE), "README.txt"))
  expect_match(paste(with_refits, collapse = " "),
               "2 models to fit, as 2 tasks, and 4 refits for cross-validation", fixed = TRUE)
  # without a remote, the commands keep placeholders a reader fills in
  plain <- paste(readLines(file.path(write_quietly(), "README.txt")), collapse = "\n")
  expect_match(plain, "ssh <user>@<cluster>", fixed = TRUE)
})


test_that("every text file of a bundle has line feeds alone, on every system", {
  dir <- write_quietly(kfold = TRUE)
  text_files <- c(scripts, "bundle/checksums.txt", "bundle/tasks-fit.csv",
                  "bundle/tasks-kfold.csv")
  for (file in text_files) {
    bytes <- readBin(file.path(dir, file), "raw", file.size(file.path(dir, file)))
    expect_false(any(bytes == as.raw(13)), label = file)       # no carriage return
  }
})


test_that("the scripts are valid shell scripts, and the R script valid R", {
  dir <- write_quietly()
  expect_no_error(parse(file.path(dir, "install.R")))
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  for (file in c("settings.sh", "install.sh", "check.sh"))
    expect_equal(system2("bash", c("-n", shQuote(file.path(dir, file)))), 0L)
})


test_that("the install job stops at once when it cannot reach the internet", {
  skip_if(!bash_reads_files() || !nzchar(Sys.which("curl")),
          "this machine's bash cannot read the files R writes, or it has no curl")
  dir <- write_quietly()
  # run as Slurm would run the job, with an address that cannot be reached
  out <- withr::with_envvar(
    c(SLURM_JOB_NAME = "rtft-install", SLURM_SUBMIT_DIR = dir,
      RTFT_CHECK_URL = "https://rtft-check.invalid"),
    suppressWarnings(system2("bash", shQuote(file.path(dir, "install.sh")),
                             stdout = TRUE, stderr = TRUE)))
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = " "), "cannot reach the internet", fixed = TRUE)
})


test_that("the source of rTFT on the cluster is read from how it was installed here", {
  from_github <- list(Version = "0.1.0.9000", RemoteType = "github",
                      RemoteUsername = "eogmartins", RemoteRepo = "rTFT",
                      RemoteSha = "3c690c9abc")
  expect_equal(.rtft_source(from_github), "eogmartins/rTFT@3c690c9abc")
  expect_equal(.rtft_source(list(Version = "0.1.0")), "eogmartins/rTFT@v0.1.0")
  expect_true(is.na(.rtft_source(list(Version = "0.1.0.9000"))))
  # and the commit it came from, which only an install from GitHub records
  expect_equal(.installed_rtft_commit(from_github), "3c690c9abc")
  expect_true(is.na(.installed_rtft_commit(list(Version = "0.1.0"))))
  expect_true(is.na(.installed_rtft_commit(list(Version = "0.1.0.9000",
                                                RemoteType = "local"))))
  # a development version from a local copy must be told where to come from
  dir <- file.path(withr::local_tempdir(), "run")
  local_mocked_bindings(.rtft_source = function(...) NA_character_)
  expect_error(write_hpc_bundle(models, long, dir, account = "def-test"),
               class = "rTFT_no_rtft_source")
  expect_error(write_hpc_bundle(models, long, dir, account = "def-test", rtft_source = ""),
               class = "rTFT_no_rtft_source")
  expect_false(dir.exists(dir))
  local_mocked_bindings(.rtft_source = function(...) "eogmartins/rTFT@v9.9.9")
  suppressMessages(write_hpc_bundle(models, long, dir, account = "def-test"))
  expect_true(any(grepl("eogmartins/rTFT@v9.9.9", readLines(file.path(dir, "install.R")),
                        fixed = TRUE)))
})


test_that("the copy is one rsync command over three lines, or ssh and scp on Windows", {
  copy <- .copy_command("/home/me/run", "me@host", "scratch/run", windows = FALSE)
  expect_equal(copy$command,
               c('rsync -av --rsync-path="mkdir -p scratch/run && rsync" \\',
                 '  "/home/me/run/" \\',
                 "  me@host:scratch/run/"))
  expect_match(copy$note, "one command over three lines", fixed = TRUE)
  on_windows <- .copy_command("C:/Users/me/run", "me@host", "scratch/projects/run",
                              windows = TRUE)
  expect_equal(on_windows$command,
               c('ssh me@host "mkdir -p scratch/projects"',
                 'scp -r "C:/Users/me/run" me@host:scratch/projects/run'))
  expect_match(on_windows$note, "WinSCP", fixed = TRUE)
  # the three lines of rsync are read by a shell as one command
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  lines <- withr::local_tempfile(fileext = ".sh")
  writeLines(copy$command, lines)
  expect_equal(system2("bash", c("-n", shQuote(lines))), 0L)
})


test_that("README.txt gives the copy command of the system the bundle is written on", {
  remote <- "someone@nibi.alliancecan.ca:scratch/run"
  unix <- readLines(file.path(write_quietly(remote = remote), "README.txt"))
  first <- grep("rsync -av", unix)
  expect_length(first, 1)
  expect_match(unix[first], "\\\\$")                     # continued on the next line
  expect_match(unix[first + 2], "someone@nibi.alliancecan.ca:scratch/run/", fixed = TRUE)
  local_mocked_bindings(.on_windows = function() TRUE)
  windows <- paste(readLines(file.path(write_quietly(remote = remote), "README.txt")),
                   collapse = "\n")
  expect_match(windows, 'ssh someone@nibi.alliancecan.ca "mkdir -p scratch"', fixed = TRUE)
  expect_match(windows, "scp -r ", fixed = TRUE)
  expect_no_match(windows, "rsync", fixed = TRUE)
})


test_that("write_hpc_bundle says what it wrote and what to do next, and returns the folder", {
  dir <- file.path(withr::local_tempdir(), "run")
  suppressMessages(
    expect_message(out <- write_hpc_bundle(models, long, dir, account = "def-test",
                                           rtft_source = "eogmartins/rTFT@dev-0.2.0"),
                   class = "rTFT_bundle_written"),
    classes = "rTFT_loo_unit")
  expect_equal(out, normalizePath(dir))
  # it points to the command in README.txt, rather than printing it wrapped
  said <- suppressMessages(
    tryCatch(write_hpc_bundle(models, long, file.path(withr::local_tempdir(), "run"),
                              account = "def-test", rtft_source = src),
             rTFT_bundle_written = function(m) conditionMessage(m)),
    classes = "rTFT_loo_unit")
  expect_match(said, "README.txt", fixed = TRUE)
  expect_no_match(said, "rsync", fixed = TRUE)
})


test_that("a bundle is never written over a folder that holds something", {
  dir <- withr::local_tempdir()
  writeLines("mine", file.path(dir, "notes.txt"))
  expect_error(write_hpc_bundle(models, long, dir, account = "def-test",
                                rtft_source = src),
               class = "rTFT_dir_not_empty")
  expect_equal(list.files(dir), "notes.txt")
  expect_error(write_hpc_bundle(models, long, c("a", "b"), account = "def-test",
                                rtft_source = src),
               class = "rTFT_bad_dir")
})


test_that("settings of the cluster that Slurm would refuse are refused here", {
  dir <- file.path(withr::local_tempdir(), "run")
  bad <- function(...) expect_error(
    suppressMessages(write_hpc_bundle(models, long, dir, ...)),
    class = "rTFT_bad_hpc_setting")
  bad()                                                   # no account
  bad(account = "")
  bad(account = c("def-a", "def-b"))
  bad(account = "def-test", time = "12h")
  bad(account = "def-test", time = "12:75:00")
  bad(account = "def-test", mem_per_cpu = "4 GB")
  bad(account = "def-test", email = "nobody")
  bad(account = "def-test", modules = character(0))
  bad(account = "def-test", r_library = c("a", "b"))
  bad(account = "def-test", remote = "nibi.alliancecan.ca")
  bad(account = "def-test", remote = "me@nibi.alliancecan.ca:scratch/my run")
  expect_false(dir.exists(dir))                           # nothing was written
})


test_that("the cores and the backend are set for the cluster, not given", {
  dir <- file.path(withr::local_tempdir(), "run")
  expect_error(write_hpc_bundle(models, long, dir, account = "def-test", rtft_source = src,
                                cores = 4),
               class = "rTFT_set_in_bundle")
  expect_error(write_hpc_bundle(models, long, dir, account = "def-test", rtft_source = src,
                                backend = "rstan"),
               class = "rTFT_set_in_bundle")
  expect_no_error(suppressMessages(
    write_hpc_bundle(models, long, dir, account = "def-test", rtft_source = src,
                     backend = "cmdstanr")))
})


test_that("the run is checked as fit_hazard_models() checks it, before anything is written", {
  dir <- file.path(withr::local_tempdir(), "run")
  mixed_units <- rbind(models[1, ], build_hazard_models("exp", "llogis"))
  expect_error(suppressMessages(write_hpc_bundle(mixed_units, long, dir,
                                                 account = "def-test", rtft_source = src)),
               class = "rTFT_units_differ")
  expect_error(suppressMessages(write_hpc_bundle(models, long, dir, account = "def-test",
                                                 rtft_source = src,
                                                 kfold = TRUE, K = 3)),
               class = "rTFT_bad_K")
  expect_error(write_hpc_bundle(models, long, dir, account = "def-test", rtft_source = src,
                                m_form = "exp"),
               class = "rTFT_set_in_models")
  expect_false(dir.exists(dir))
})


test_that("two models that differ only in their names are refused", {
  dir <- file.path(withr::local_tempdir(), "run")
  again <- models[1, ]
  again$model <- "exp_mem_again"
  expect_error(suppressMessages(write_hpc_bundle(rbind(models, again), long, dir,
                                                 account = "def-test",
                                                 rtft_source = src)),
               class = "rTFT_duplicate_models")
  said <- tryCatch(suppressMessages(write_hpc_bundle(
    rbind(models, again), long, dir, account = "def-test", rtft_source = src)),
    error = function(e) cli::ansi_strip(conditionMessage(e)))
  expect_match(said, '"exp_mem" and "exp_mem_again"', fixed = TRUE)
  expect_false(dir.exists(dir))
  # the same function and baseline with another formula, or with a co-stressor, is
  # another model
  other <- again
  other$dpar_formulas <- list(list(loga = ~ 1 + mass + (1 | trial)))
  expect_true(.check_distinct_models(rbind(models, other)))
  stressed <- again
  stressed$costressors <- list(list(costressor("DO", column = "do", ref = 8.5, scale = 1,
                                               direction = -1)))
  expect_true(.check_distinct_models(rbind(models, stressed)))
})


test_that("a template of one's own is used, and one that cannot be filled is refused", {
  own <- withr::local_tempdir()
  for (file in scripts)
    writeLines(c("# mine", "ACCOUNT={{account}}"), file.path(own, file))
  dir <- write_quietly(template = own)
  expect_equal(readLines(file.path(dir, "settings.sh")), c("# mine", "ACCOUNT=def-test"))
  # a placeholder the bundle has no value for
  writeLines("PARTITION={{partition}}", file.path(own, "check.sh"))
  half <- file.path(withr::local_tempdir(), "run")
  expect_error(write_quietly(template = own, dir = half), class = "rTFT_template_unfilled")
  expect_false(dir.exists(half))                          # nothing was written
  # a folder without the three files, or something that is not a folder
  expect_error(write_quietly(template = withr::local_tempdir()), class = "rTFT_bad_template")
  expect_error(write_quietly(template = 1), class = "rTFT_bad_template")
})


# ===========================================================================
# 3. Checking the bundle (check_hpc_bundle)
# ===========================================================================

test_that("check_hpc_bundle passes a bundle as it was written, and says what it checked", {
  with_cmdstan()
  dir <- write_quietly(kfold = TRUE)
  expect_message(out <- check_hpc_bundle(dir), class = "rTFT_bundle_checked")
  expect_true(out)
  # from the folder itself, as check.sh runs it
  withr::local_dir(dir)
  expect_true(suppressMessages(check_hpc_bundle()))
})


test_that("a bundle with a file missing, changed or added is refused", {
  with_cmdstan()
  dir <- write_quietly(kfold = TRUE)
  # missing
  file.rename(file.path(dir, "bundle", "models.rds"), file.path(dir, "models.rds"))
  expect_error(check_hpc_bundle(dir), class = "rTFT_bad_bundle")
  file.rename(file.path(dir, "models.rds"), file.path(dir, "bundle", "models.rds"))
  expect_true(suppressMessages(check_hpc_bundle(dir)))
  # changed: a copy cut short
  data_file <- file.path(dir, "bundle", "data.rds")
  whole <- readBin(data_file, "raw", file.size(data_file))
  writeBin(whole[seq_len(length(whole) %/% 2)], data_file)
  expect_error(check_hpc_bundle(dir), class = "rTFT_bundle_changed")
  writeBin(whole, data_file)
  # gone, though recorded
  file.rename(file.path(dir, "bundle", "tasks-kfold.csv"), file.path(dir, "tasks.csv"))
  expect_error(check_hpc_bundle(dir), class = "rTFT_bundle_changed")
  file.rename(file.path(dir, "tasks.csv"), file.path(dir, "bundle", "tasks-kfold.csv"))
  # added since
  writeLines("x", file.path(dir, "bundle", "extra.csv"))
  expect_error(check_hpc_bundle(dir), class = "rTFT_bundle_changed")
  # and a folder that is not a bundle at all
  expect_error(check_hpc_bundle(withr::local_tempdir()), class = "rTFT_bad_bundle")
})


test_that("check_hpc_bundle refuses another version of rTFT, and missing software", {
  dir <- write_quietly()
  local_mocked_bindings(.has_cmdstanr = function() TRUE, .cmdstan_found = function() TRUE,
                        .installed_rtft_version = function() "0.0.1")
  expect_error(check_hpc_bundle(dir), class = "rTFT_version_differs")
  local_mocked_bindings(.installed_rtft_version = function()
    as.character(utils::packageVersion("rTFT")))
  local_mocked_bindings(.has_cmdstanr = function() FALSE)
  expect_error(check_hpc_bundle(dir), class = "rTFT_needs_cmdstanr")
  local_mocked_bindings(.has_cmdstanr = function() TRUE, .cmdstan_found = function() FALSE)
  expect_error(check_hpc_bundle(dir), class = "rTFT_no_cmdstan")
})


test_that("copies of one version are told apart by their commits, where both have one", {
  with_cmdstan()
  sha_a <- "3c690c9abc0000000000000000000000000000aa"
  sha_b <- "e15b42aef00000000000000000000000000000bb"
  local_mocked_bindings(.installed_rtft_commit = function(...) sha_a)
  dir <- write_quietly()
  expect_equal(readRDS(file.path(dir, "bundle", "settings.rds"))$rtft_commit, sha_a)
  # the install reinstalls when the commit differs, and checks what it installed
  install <- readLines(file.path(dir, "install.R"))
  expect_true(any(grepl(paste0('^wanted_commit <- "', sha_a, '"$'), install)))
  expect_no_error(parse(text = install))
  expect_true(suppressMessages(check_hpc_bundle(dir)))
  # another commit of the same version is refused, here and at the start of a task
  local_mocked_bindings(.installed_rtft_commit = function(...) sha_b)
  expect_error(check_hpc_bundle(dir), class = "rTFT_version_differs")
  expect_match(tryCatch(check_hpc_bundle(dir), error = conditionMessage), "e15b42a",
               fixed = TRUE)
  expect_error(suppressMessages(run_hpc_task(dir, "fit", 1)), class = "rTFT_task_failed")
  expect_match(.read_record(dir, "fit", 1)$message, "from commit \"e15b42a\"", fixed = TRUE)
  # a copy without a commit, installed from a folder, is checked by its number alone
  local_mocked_bindings(.installed_rtft_commit = function(...) NA_character_)
  expect_true(suppressMessages(check_hpc_bundle(dir)))
  # and so is a bundle that recorded none
  plain <- write_quietly()
  expect_true(is.na(readRDS(file.path(plain, "bundle", "settings.rds"))$rtft_commit))
  expect_true('wanted_commit <- ""' %in% readLines(file.path(plain, "install.R")))
  local_mocked_bindings(.installed_rtft_commit = function(...) sha_b)
  expect_true(suppressMessages(check_hpc_bundle(plain)))
})


# ===========================================================================
# 4. The helpers
# ===========================================================================

test_that("the settings Slurm takes are accepted in each of their forms", {
  ok <- function(...) expect_true(.check_hpc_settings(...))
  ok("def-a", "12:00:00", "4G", NULL, "r/4.5.0", NULL, NULL)
  ok("def-a", "1-00:00:00", "4000M", "me@example.org", c("StdEnv/2023", "r/4.5.0"),
     "~/R/lib", "me@host:scratch/run")
  ok("def-a", "100:00:00", "512", NULL, "r", NULL, NULL)
})


test_that(".fill_template fills every placeholder it has a value for", {
  expect_equal(.fill_template(c("A={{a}}", "B={{b}} and {{a}}"), list(a = "1", b = 2)),
               c("A=1", "B=2 and 1"))
  expect_error(.fill_template("C={{c}}", list(a = "1")), class = "rTFT_template_unfilled")
  # on any line of the script, not only the first
  expect_error(.fill_template(c("A={{a}}", "C={{c}}"), list(a = "1")),
               class = "rTFT_template_unfilled")
})


test_that("what is installed is read from the machine itself", {
  expect_equal(.installed_rtft_version(), as.character(utils::packageVersion("rTFT")))
  expect_type(.has_cmdstanr(), "logical")
  expect_length(.has_cmdstanr(), 1)
  skip_if_not_installed("cmdstanr")
  expect_type(.cmdstan_found(), "logical")
  expect_length(.cmdstan_found(), 1)
})

test_that("a file is written whole or not at all", {
  dir  <- withr::local_tempdir()
  path <- file.path(dir, "a.rds")
  expect_equal(.write_whole(1:3, path, saveRDS), path)
  expect_equal(readRDS(path), 1:3)
  expect_false(file.exists(paste0(path, ".part")))
  # a file cannot be renamed onto a folder that holds something
  taken <- file.path(dir, "taken")
  dir.create(taken)
  writeLines("x", file.path(taken, "x"))
  expect_error(.write_whole(1:3, taken, saveRDS), class = "rTFT_not_written")
})


test_that("a task is named the same way in its own messages and in the status", {
  fits <- data.frame(task = 2L, model = "exp_weibull")
  expect_equal(.task_label("fit", fits), "fit 2 (exp_weibull)")
  mixed <- data.frame(task = c(1L, 3L), model = "exp_mem", fold = c(NA, 3L))
  expect_equal(.task_label(c("fit", "kfold"), mixed),
               c("fit 1 (exp_mem)", "refit 3 (exp_mem, fold 3)"))
})


test_that("a reason is kept as one line of ASCII text, and a time in its unit", {
  made_by_cli <- tryCatch(cli::cli_abort(c("Only 1 of the {2} chains finished.",
                                           "i" = "See {.file logs/}.")),
                          error = identity)
  expect_equal(.one_line(made_by_cli), "Only 1 of the 2 chains finished. i See 'logs/'.")
  # base R quotes names with curly quotes where R runs in UTF-8
  quoted <- simpleError("there is no package called \u2018brms\u2019 \u201cnow\u201d")
  expect_equal(.one_line(quoted), "there is no package called 'brms' \"now\"")
  # the unit is chosen after rounding
  expect_equal(vapply(c(14, 119.6, 840, 7199, 41000), .duration, ""),
               c("14 s", "2 min", "14 min", "2 h", "11.4 h"))
})



# ===========================================================================
# 5. Running a task (run_hpc_task)
# ===========================================================================

test_that("run_hpc_task fits the model of its task as fit_hazard_models() does, and keeps it", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  before <- getOption("cmdstanr_write_stan_file_dir")
  expect_message(expect_message(path <- run_hpc_task(dir, "fit", 2),
                                class = "rTFT_task_progress"),
                 class = "rTFT_task_done")
  expect_equal(path, .task_file(dir, "fit", 2, "rds"))
  # the model of the task, with the arguments of the bundle and a core per chain
  expect_length(record$calls, 1)
  expect_equal(record$calls[[1]]$baseline, "weibull")
  expect_mapequal(record$calls[[1]]$args,
                  list(chains = 2, seed = 1, sd_prior = "exponential(1)",
                       backend = "cmdstanr", cores = 2))
  expect_identical(getOption("cmdstanr_write_stan_file_dir"), before)
  # the result, stamped with the bundle
  res <- readRDS(path)
  expect_s3_class(res$fit, "brmsfit")
  expect_s3_class(res$loo, "loo")
  expect_equal(res$loo$unit, "trial")
  expect_equal(res$model, "exp_weibull")
  expect_equal(res$bundle, .bundle_id(dir))
  expect_true(all(c("convergence", "high_k", "seconds") %in% names(res)))
  # and the records of its start and of its end
  ended <- .read_record(dir, "fit", 2)
  expect_equal(ended$status, "done")
  expect_equal(ended$model, "exp_weibull")
  expect_equal(ended$bundle, .bundle_id(dir))
  expect_true(is.na(ended$message))
  expect_true(file.exists(.task_file(dir, "fit", 2, "started")))
  expect_length(list.files(file.path(dir, "results"), pattern = "[.]part$"), 0)
})


test_that("a task that has ended well is not run again, and one that failed is", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  run_quietly(dir, "fit", 1)
  expect_message(run_hpc_task(dir, "fit", 1), class = "rTFT_task_skipped")
  expect_length(record$calls, 1)
  # a task that fails says so, and keeps its reason
  record$fail <- TRUE
  expect_error(run_quietly(dir, "fit", 2), class = "rTFT_task_failed")
  ended <- .read_record(dir, "fit", 2)
  expect_equal(ended$status, "failed")
  expect_match(ended$message, "the sampler stopped", fixed = TRUE)
  expect_false(file.exists(.task_file(dir, "fit", 2, "rds")))
  # and is run again
  record$fail <- FALSE
  run_quietly(dir, "fit", 2)
  expect_length(record$calls, 3)
  expect_equal(.read_record(dir, "fit", 2)$status, "done")
})


test_that("the number of the task is read from Slurm", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  withr::local_envvar(SLURM_ARRAY_TASK_ID = "2", SLURM_ARRAY_JOB_ID = "77")
  run_quietly(dir)
  expect_equal(record$calls[[1]]$baseline, "weibull")
  started <- utils::read.csv(.task_file(dir, "fit", 2, "started"), colClasses = "character")
  expect_equal(started$job, "77_2")
  expect_equal(.read_record(dir, "fit", 2)$job, "77_2")
})


test_that("a model keeps its name through the tables of the tasks, whatever it is", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record),
                        .refit_fold = stand_in_refit(new_record()))
  odd <- build_hazard_models("exp", c("mem", "weibull"),
                             dpar_formulas = list(loga = ~ 1 + (1 | trial)),
                             model_names = c("NA", "1"))
  dir <- trial_bundle(models = odd)
  tasks <- .read_tasks(dir, "kfold")
  expect_identical(tasks$model, rep(c("NA", "1"), each = 3))
  expect_type(tasks$task, "integer")
  expect_type(tasks$fold, "integer")
  run_all(dir)
  # each task fitted its own model, and kept its name
  expect_equal(vapply(record$calls, `[[`, "", "baseline"), c("mem", "weibull"))
  expect_identical(readRDS(.task_file(dir, "fit", 1, "rds"))$model, "NA")
})


test_that("a refit starts from its model's fit and leaves out the units of its fold", {
  refits <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(refits))
  dir <- trial_bundle()
  settings <- readRDS(file.path(dir, "bundle", "settings.rds"))
  tasks <- .read_tasks(dir, "kfold")
  # before its fit has ended, a refit cannot start
  expect_error(.kfold_task(dir, settings, tasks[1, ], .bundle_id(dir)),
               class = "rTFT_fit_not_ready")
  expect_error(run_quietly(dir, "kfold", 1), class = "rTFT_task_failed")
  expect_match(.read_record(dir, "kfold", 1)$message, "start from its fit", fixed = TRUE)
  expect_length(refits$calls, 0)
  # once it has, each refit leaves out its own fold, and is given what brms takes
  run_quietly(dir, "fit", 1)
  for (i in 1:3) run_quietly(dir, "kfold", i)
  expect_length(refits$calls, 3)
  for (k in 1:3) {
    held <- names(settings$folds)[settings$folds == k]
    expect_setequal(refits$calls[[k]]$held, held)
    expect_mapequal(refits$calls[[k]]$args,
                    list(chains = 2, seed = 1, backend = "cmdstanr", cores = 2))
    res <- readRDS(.task_file(dir, "kfold", k, "rds"))
    expect_setequal(names(res$elpd), held)
    expect_equal(res$fold, k)
    expect_equal(res$bundle, .bundle_id(dir))
  }
})


test_that("run_hpc_task refuses a task that is not one, before fitting anything", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  expect_error(run_hpc_task(dir, "refit", 1), class = "rTFT_bad_stage")
  expect_error(run_hpc_task(dir, c("fit", "kfold"), 1), class = "rTFT_bad_stage")
  for (bad in list(0, 3, "", "a", c(1, 2)))
    expect_error(run_hpc_task(dir, "fit", bad), class = "rTFT_no_task")
  expect_error(run_hpc_task(dir, "kfold", 7), class = "rTFT_no_task")
  withr::local_envvar(SLURM_ARRAY_TASK_ID = NA)
  expect_error(run_hpc_task(dir, "fit"), class = "rTFT_no_task")
  # outside Slurm, with no task given, the message says to give one
  expect_match(tryCatch(run_hpc_task(dir, "kfold"), error = conditionMessage),
               'run_hpc_task(stage = "kfold", task = 1)', fixed = TRUE)
  # a bundle without refits has no task for them
  expect_error(run_hpc_task(trial_bundle(kfold = FALSE), "kfold", 1),
               class = "rTFT_no_task")
  expect_length(record$calls, 0)
  expect_false(file.exists(.task_file(dir, "fit", 1, "started")))
})


test_that("run_hpc_task refuses a bundle that changed, or another version of rTFT", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  writeLines("x", file.path(dir, "bundle", "extra.csv"))
  expect_error(run_hpc_task(dir, "fit", 1), class = "rTFT_bundle_changed")
  file.remove(file.path(dir, "bundle", "extra.csv"))
  expect_false(file.exists(.task_file(dir, "fit", 1, "started")))
  # another version stops the task before it fits, and the reason is recorded
  local_mocked_bindings(.installed_rtft_version = function() "0.0.1")
  expect_error(run_quietly(dir, "fit", 1), class = "rTFT_task_failed")
  expect_length(record$calls, 0)
  ended <- .read_record(dir, "fit", 1)
  expect_equal(ended$status, "failed")
  expect_equal(ended$rtft_version, "0.0.1")
  expect_match(ended$message, "The installed rTFT is version \"0.0.1\"", fixed = TRUE)
})


test_that("a task run again keeps the record of its earlier end aside", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  record$fail <- TRUE
  try(run_quietly(dir, "fit", 2), silent = TRUE)
  expect_equal(.read_record(dir, "fit", 2)$status, "failed")
  # while the second attempt runs (seen from inside the fit), the task has no record
  # of an end
  record$fail <- FALSE
  during <- NULL
  local_mocked_bindings(fit_hazard = function(...) {
    during <<- .read_record(dir, "fit", 2)
    fit_trials
  })
  run_quietly(dir, "fit", 2)
  expect_null(during)
  # the earlier end is beside it, with its reason
  previous <- utils::read.csv(.task_file(dir, "fit", 2, "previous.csv"),
                              colClasses = "character")
  expect_equal(previous$status, "failed")
  expect_match(previous$message, "the sampler stopped", fixed = TRUE)
  expect_equal(.read_record(dir, "fit", 2)$status, "done")
})
