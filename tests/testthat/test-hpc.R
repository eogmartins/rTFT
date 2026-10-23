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
scripts <- c("settings.sh", "install.sh", "install.R", "check.sh", "submit.sh", "task.sh",
             "status.sh", "README.txt")

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

# Write the record that a task started, as run_hpc_task() writes it, as the Slurm job
# `job`, and the log Slurm wrote for it, holding `log`.
fake_start <- function(dir, stage, task, job, log = NULL) {
  .write_csv_lf(data.frame(stage = stage, task = task, job = job, host = "nc1",
                           started = "2026-10-17 10:00:00 PDT"),
                .task_file(dir, stage, task, "started"))
  if (!is.null(log))
    writeLines(log, file.path(dir, "logs", sprintf("%s-%s.out", stage, job)))
}

# The status of each task, as hpc_status() gives it, without its message.
status_quietly <- function(dir) suppressMessages(hpc_status(dir))

# The message of a condition of class `class` that an expression gives, as plain text.
message_of <- function(expr, class) {
  said <- NULL
  suppressWarnings(suppressMessages(withCallingHandlers(expr, condition = function(cnd) {
    if (inherits(cnd, class)) said <<- cli::ansi_strip(conditionMessage(cnd))
  })))
  said
}

# Check that `said` gives `command` to run on the cluster as plain text, ready to paste:
# not between backticks, which cli puts around code where the terminal has no colours
# and bash would run, and not followed by the full stop of a sentence, which bash would
# read as part of the last argument.
expect_plain_command <- function(said, command) {
  expect_match(said, command, fixed = TRUE)
  expect_no_match(said, paste0("`", command), fixed = TRUE)
  expect_no_match(said, paste0(command, "."), fixed = TRUE)
}

# The message of an error, as plain text.
said_by <- function(err) cli::ansi_strip(conditionMessage(err))


# ===========================================================================
# 2. Writing the bundle (write_hpc_bundle)
# ===========================================================================

test_that("write_hpc_bundle writes the scripts, the inputs and the empty folders", {
  dir <- write_quietly(chains = 2, seed = 1)
  expect_equal(bundle_files(dir),
               c("README.txt", "check.sh", "inputs/checksums.txt", "inputs/data.rds",
                 "inputs/models.rds", "inputs/settings.rds", "inputs/tasks-fit.csv",
                 "install.R", "install.sh", "settings.sh", "status.sh", "submit.sh",
                 "task.sh"))
  expect_true(all(dir.exists(file.path(dir, c("stan", "results", "logs")))))
  # without cross-validation by refitting there is no table of refits
  expect_false(file.exists(file.path(dir, "inputs", "tasks-kfold.csv")))
})


test_that("the bundle holds the data, the models and the settings of the run", {
  dir <- write_quietly(chains = 2, iter = 1000, seed = 1, time = "02:30:00",
                       mem_per_cpu = "2000M", email = "me@example.org")
  expect_identical(readRDS(file.path(dir, "inputs", "data.rds")), long)
  # the formulas come back with environments of their own, so only their content is
  # compared
  expect_equal(readRDS(file.path(dir, "inputs", "models.rds")), models,
               ignore_formula_env = TRUE)
  s <- readRDS(file.path(dir, "inputs", "settings.rds"))
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
  tasks <- utils::read.csv(file.path(dir, "inputs", "tasks-fit.csv"))
  expect_equal(tasks, data.frame(task = 1:2, model = c("exp_mem", "exp_weibull"),
                                 m_form = "exp", baseline = c("mem", "weibull")))
})


test_that("the checksums of the bundle are those of its files", {
  dir <- write_quietly(kfold = TRUE)
  sums <- utils::read.table(file.path(dir, "inputs", "checksums.txt"),
                            col.names = c("md5", "file"), colClasses = "character")
  expect_equal(sums$file, paste0("inputs/", c("data.rds", "models.rds", "settings.rds",
                                              "tasks-fit.csv", "tasks-kfold.csv")))
  expect_equal(sums$md5, unname(tools::md5sum(file.path(dir, sums$file))))
})


test_that("with kfold = TRUE the folds are assigned once and each refit is a task", {
  set.seed(4)
  dir <- write_quietly(kfold = TRUE)
  set.seed(4)
  expected <- .assign_folds(c("a", "b"), 2)
  s <- readRDS(file.path(dir, "inputs", "settings.rds"))
  expect_true(s$kfold)
  expect_equal(s$K, 2)
  expect_equal(s$folds, expected)
  tasks <- utils::read.csv(file.path(dir, "inputs", "tasks-kfold.csv"))
  expect_equal(tasks$task, 1:4)
  expect_equal(tasks$model, rep(c("exp_mem", "exp_weibull"), each = 2))
  expect_equal(tasks$fold, rep(1:2, 2))
  expect_equal(tasks$units, rep(names(expected)[order(expected)], 2))
})


test_that("the scripts are filled in from the template, with nothing left to fill", {
  # the commands checked here are those written on macOS and Linux; Windows has its
  # own, checked below
  local_mocked_bindings(.on_windows = function() FALSE)
  dir <- write_quietly(chains = 3, time = "1-00:00:00", modules = c("StdEnv/2023", "r/4.4.0"),
                       r_library = "~/R/rtft",
                       remote = "someone@nibi.alliancecan.ca:scratch/run")
  settings <- readLines(file.path(dir, "settings.sh"))
  expect_true(all(c('ACCOUNT="def-test"', 'TIME="1-00:00:00"', 'MEM_PER_CPU="4000M"',
                    'CPUS_PER_TASK="3"', 'MODULES="StdEnv/2023 r/4.4.0"',
                    'R_LIBRARY="~/R/rtft"', 'EMAIL=""') %in% settings))
  readme <- paste(readLines(file.path(dir, "README.txt")), collapse = "\n")
  expect_match(readme, "ssh someone@nibi.alliancecan.ca", fixed = TRUE)
  expect_match(readme, "cd scratch/run", fixed = TRUE)
  expect_match(readme, "someone@nibi.alliancecan.ca:scratch/run/", fixed = TRUE)
  # the copy makes the whole path of the folder on the cluster first
  expect_match(readme, 'rsync -av --rsync-path="mkdir -p scratch/run && rsync"', fixed = TRUE)
  # and the results come back into the bundle, to be assembled there
  expect_match(readme, "someone@nibi.alliancecan.ca:scratch/run/results/ \\", fixed = TRUE)
  # the paths are written with forward slashes on every system
  expect_match(readme, paste0('"', normalizePath(dir, winslash = "/"), '/results/"'),
               fixed = TRUE)
  expect_match(readme, paste0('rTFT::collect_hpc("', normalizePath(dir, winslash = "/"),
                              '")'), fixed = TRUE)
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
  text_files <- c(scripts, "inputs/checksums.txt", "inputs/tasks-fit.csv",
                  "inputs/tasks-kfold.csv")
  for (file in text_files) {
    bytes <- readBin(file.path(dir, file), "raw", file.size(file.path(dir, file)))
    expect_false(any(bytes == as.raw(13)), label = file)       # no carriage return
  }
})


test_that("the scripts are valid shell scripts, and the R script valid R", {
  dir <- write_quietly()
  expect_no_error(parse(file.path(dir, "install.R")))
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  for (file in c("settings.sh", "install.sh", "check.sh", "submit.sh", "task.sh",
                 "status.sh"))
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
  local_mocked_bindings(.on_windows = function() FALSE)
  unix <- readLines(file.path(write_quietly(remote = remote), "README.txt"))
  first <- grep("rsync -av --rsync-path", unix, fixed = TRUE)
  expect_length(first, 1)
  expect_match(unix[first], "\\\\$")                     # continued on the next line
  expect_match(unix[first + 2], "someone@nibi.alliancecan.ca:scratch/run/", fixed = TRUE)
  local_mocked_bindings(.on_windows = function() TRUE)
  windows <- paste(readLines(file.path(write_quietly(remote = remote), "README.txt")),
                   collapse = "\n")
  expect_match(windows, 'ssh someone@nibi.alliancecan.ca "mkdir -p scratch"', fixed = TRUE)
  expect_match(windows, "scp -r ", fixed = TRUE)
  expect_no_match(windows, "rsync", fixed = TRUE)
  expect_match(windows, "scp -r someone@nibi.alliancecan.ca:scratch/run/results ", fixed = TRUE)
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
  expect_plain_command(said, "bash install.sh")
  expect_plain_command(said, "bash check.sh")
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


# Write the next round of the bundle in `dir`, without its messages, into a new
# temporary folder, and return the folder.
round_quietly <- function(dir, models, ..., env = parent.frame()) {
  new_dir <- file.path(withr::local_tempdir(.local_envir = env), "round")
  suppressMessages(next_hpc_round(dir, models, new_dir, rtft_source = src, ...))
  new_dir
}


test_that("next_hpc_round copies what makes the rounds comparable, and changes the sampler", {
  dir <- trial_bundle()
  first <- readRDS(file.path(dir, "inputs", "settings.rds"))
  expect_message(
    next_hpc_round(dir, "exp_weibull", file.path(withr::local_tempdir(), "round"),
                   iter = 4000, rtft_source = src),
    class = "rTFT_round_written")
  round <- round_quietly(dir, "exp_weibull", iter = 4000,
                         control = list(adapt_delta = 0.99), time = "24:00:00")
  # a bundle of its own, laid out as the first
  expect_equal(bundle_files(round), bundle_files(dir)[!grepl("^results/|^logs/",
                                                              bundle_files(dir))])
  # the data file itself, not written again
  expect_equal(unname(tools::md5sum(file.path(round, "inputs", "data.rds"))),
               unname(tools::md5sum(file.path(dir, "inputs", "data.rds"))))
  # the whole table of models, and tasks for the models named, and only those
  expect_equal(readRDS(file.path(round, "inputs", "models.rds"))$model,
               c("exp_mem", "exp_weibull"))
  expect_equal(.read_tasks(round, "fit")$model, "exp_weibull")
  # copied: the unit, the folds, the rule, the observation, the account and modules
  settings <- readRDS(file.path(round, "inputs", "settings.rds"))
  for (what in c("observation", "unit", "units", "kfold", "K", "folds", "rule",
                 "rtft_version"))
    expect_identical(settings[[what]], first[[what]])
  expect_identical(settings$slurm[c("account", "modules", "r_library", "mem_per_cpu")],
                   first$slurm[c("account", "modules", "r_library", "mem_per_cpu")])
  # changed: the sampler and the time; the prior of the first bundle is kept
  expect_equal(settings$dots$iter, 4000)
  expect_equal(settings$dots$control$adapt_delta, 0.99)
  expect_equal(settings$dots$seed, first$dots$seed)
  expect_equal(settings$dots$sd_prior, "exponential(1)")
  expect_equal(settings$slurm$time, "24:00:00")
  # its scripts say what wrote them
  here <- normalizePath(dir, winslash = "/")
  readme <- readLines(file.path(round, "README.txt"))
  expect_match(readme[4], "Written by rTFT::next_hpc_round() on", fixed = TRUE)
  expect_match(paste(readme[5:8], collapse = " "),
               paste0("It is round 2 of the bundle in ", here, ". It holds"),
               fixed = TRUE)
  expect_match(readLines(file.path(dir, "README.txt"))[4], "write_hpc_bundle()",
               fixed = TRUE)
  # its place in the chain of rounds
  expect_equal(settings$round, list(number = 2L, first = .bundle_id(dir),
                                    previous = .bundle_id(dir), first_dir = here,
                                    previous_dir = here))
  # the refits of the model, on the folds of the first bundle
  tasks <- .read_tasks(round, "kfold")
  of_first <- .read_tasks(dir, "kfold")
  expect_equal(tasks$units, of_first$units[of_first$model == "exp_weibull"])
  expect_true(all(tasks$model == "exp_weibull"))
})


test_that("a round merges control into the bundle's, and changes the chains with the cores", {
  dir <- trial_bundle()
  first <- readRDS(file.path(dir, "inputs", "settings.rds"))
  first$dots$control <- list(max_treedepth = 12)
  saveRDS(first, file.path(dir, "inputs", "settings.rds"))
  sums <- .bundle_checksums(dir)
  .write_lines_lf(paste(sums, names(sums), sep = "  "),
                  file.path(dir, "inputs", "checksums.txt"))
  round <- round_quietly(dir, "exp_mem", chains = 4, control = list(adapt_delta = 0.99))
  settings <- readRDS(file.path(round, "inputs", "settings.rds"))
  expect_equal(settings$dots$control, list(max_treedepth = 12, adapt_delta = 0.99))
  expect_equal(settings$chains, 4)
  expect_equal(settings$slurm$cpus_per_task, 4)
  expect_match(readLines(file.path(round, "settings.sh")), 'CPUS_PER_TASK="4"',
               fixed = TRUE, all = FALSE)
})


test_that("each round follows the one before, and counts from the first bundle", {
  dir <- trial_bundle()
  two   <- round_quietly(dir, c("exp_mem", "exp_weibull"), iter = 3000)
  three <- round_quietly(two, "exp_weibull", iter = 6000)
  settings <- readRDS(file.path(three, "inputs", "settings.rds"))
  expect_equal(settings$round, list(
    number = 3L, first = .bundle_id(dir), previous = .bundle_id(two),
    first_dir = normalizePath(dir, winslash = "/"),
    previous_dir = normalizePath(two, winslash = "/")))
  expect_equal(settings$dots$iter, 6000)
  readme <- paste(readLines(file.path(three, "README.txt"))[5:9], collapse = " ")
  expect_match(readme, paste0("It is round 3 of the bundle in ",
                              normalizePath(dir, winslash = "/"), ", after ",
                              normalizePath(two, winslash = "/"), "."), fixed = TRUE)
  # a model left out of a round can be fitted again in the next
  only_one <- round_quietly(dir, "exp_weibull")
  back <- round_quietly(only_one, "exp_mem")
  expect_equal(.read_tasks(back, "fit")$model, "exp_mem")
  expect_equal(readRDS(file.path(back, "inputs", "settings.rds"))$round$number, 3L)
})


test_that("a round is written beside the bundle by default, named after the first", {
  dir <- trial_bundle()
  two <- suppressMessages(next_hpc_round(dir, "exp_mem", rtft_source = src))
  expect_equal(two, normalizePath(file.path(dirname(dir), "run-round-2")))
  three <- suppressMessages(next_hpc_round(two, "exp_mem", rtft_source = src))
  expect_equal(three, normalizePath(file.path(dirname(dir), "run-round-3")))
  # the folder of a round is refused when it holds something, as any other
  expect_error(suppressMessages(next_hpc_round(dir, "exp_mem", rtft_source = src)),
               class = "rTFT_dir_not_empty")
  # a folder that is not a string is refused by the name of its argument
  said <- tryCatch(next_hpc_round(dir, "exp_mem", 1, rtft_source = src),
                   rTFT_bad_dir = function(e) cli::ansi_strip(conditionMessage(e)))
  expect_match(said, "`new_dir` must be", fixed = TRUE)
  expect_error(next_hpc_round(NULL, "exp_mem"), class = "rTFT_bad_dir")
})


test_that("next_hpc_round refuses what would make another model or another run", {
  dir <- trial_bundle()
  # anything but the sampler
  expect_error(round_quietly(dir, "exp_mem", sd_prior = "exponential(2)"),
               class = "rTFT_not_sampler_setting")
  expect_error(round_quietly(dir, "exp_mem", k_sw = 5),
               class = "rTFT_not_sampler_setting")
  expect_error(round_quietly(dir, "exp_mem", backend = "rstan"),
               class = "rTFT_not_sampler_setting")
  expect_error(round_quietly(dir, "exp_mem", 4000), class = "rTFT_not_sampler_setting")
  # models that are not the bundle's, each once
  for (bad in list("exp_gompertz", character(0), c("exp_mem", "exp_mem"), NA, 1))
    expect_error(round_quietly(dir, bad), class = "rTFT_bad_round_models")
  # another version of rTFT
  local({
    local_mocked_bindings(.installed_rtft_version = function() "0.0.1")
    expect_error(round_quietly(dir, "exp_mem"), class = "rTFT_round_version_differs")
  })
  # settings of the jobs that Slurm cannot take
  expect_error(round_quietly(dir, "exp_mem", time = "a day"),
               class = "rTFT_bad_hpc_setting")
  # a folder that holds something already
  full <- withr::local_tempdir()
  writeLines("x", file.path(full, "note.txt"))
  expect_error(
    suppressMessages(next_hpc_round(dir, "exp_mem", full, rtft_source = src)),
    class = "rTFT_dir_not_empty")
  # a bundle changed since it was written
  writeLines("x", file.path(dir, "inputs", "extra.csv"))
  expect_error(round_quietly(dir, "exp_mem"), class = "rTFT_bundle_changed")
})


test_that("the tasks of a round fit with its own settings, and refit from its own fit", {
  record <- new_record()
  refits <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record),
                        .refit_fold = stand_in_refit(refits))
  with_cmdstan()
  dir   <- trial_bundle()
  round <- round_quietly(dir, "exp_weibull", iter = 4000)
  expect_true(suppressMessages(check_hpc_bundle(round)))
  # it assembles the model the round fits, not the whole table it holds
  expect_match(message_of(check_hpc_bundle(round), "rTFT_bundle_checked"),
               "1 model assembled", fixed = TRUE)
  run_quietly(round, "fit", 1)
  expect_equal(record$calls[[1]]$baseline, "weibull")
  expect_equal(record$calls[[1]]$args$iter, 4000)
  expect_equal(record$calls[[1]]$args$sd_prior, "exponential(1)")
  run_quietly(round, "kfold", 1)
  expect_equal(refits$calls[[1]]$args$iter, 4000)
  first_fold <- .read_tasks(dir, "kfold")$units[1]
  expect_setequal(refits$calls[[1]]$held, strsplit(first_fold, ";")[[1]])
})


# ===========================================================================
# 3. Checking the bundle (check_hpc_bundle)
# ===========================================================================

test_that("check_hpc_bundle passes a bundle as it was written, and says what it checked", {
  with_cmdstan()
  dir <- write_quietly(kfold = TRUE)
  expect_message(out <- check_hpc_bundle(dir), class = "rTFT_bundle_checked")
  expect_true(out)
  expect_plain_command(message_of(check_hpc_bundle(dir), "rTFT_bundle_checked"),
                       "bash submit.sh")
  # from the folder itself, as check.sh runs it
  withr::local_dir(dir)
  expect_true(suppressMessages(check_hpc_bundle()))
})


test_that("a bundle with a file missing, changed or added is refused", {
  with_cmdstan()
  dir <- write_quietly(kfold = TRUE)
  # missing
  file.rename(file.path(dir, "inputs", "models.rds"), file.path(dir, "models.rds"))
  expect_error(check_hpc_bundle(dir), class = "rTFT_bad_bundle")
  file.rename(file.path(dir, "models.rds"), file.path(dir, "inputs", "models.rds"))
  expect_true(suppressMessages(check_hpc_bundle(dir)))
  # changed: a copy cut short
  data_file <- file.path(dir, "inputs", "data.rds")
  whole <- readBin(data_file, "raw", file.size(data_file))
  writeBin(whole[seq_len(length(whole) %/% 2)], data_file)
  expect_error(check_hpc_bundle(dir), class = "rTFT_bundle_changed")
  writeBin(whole, data_file)
  # gone, though recorded
  file.rename(file.path(dir, "inputs", "tasks-kfold.csv"), file.path(dir, "tasks.csv"))
  expect_error(check_hpc_bundle(dir), class = "rTFT_bundle_changed")
  file.rename(file.path(dir, "tasks.csv"), file.path(dir, "inputs", "tasks-kfold.csv"))
  # added since
  writeLines("x", file.path(dir, "inputs", "extra.csv"))
  expect_error(check_hpc_bundle(dir), class = "rTFT_bundle_changed")
  # and a folder that is not a bundle at all
  expect_error(check_hpc_bundle(withr::local_tempdir()), class = "rTFT_bad_bundle")
})


test_that("check_hpc_bundle refuses another version of rTFT, and missing software", {
  dir <- write_quietly()
  local_mocked_bindings(.has_cmdstanr = function() TRUE, .cmdstan_found = function() TRUE,
                        .installed_rtft_version = function() "0.0.1")
  err <- expect_error(check_hpc_bundle(dir), class = "rTFT_version_differs")
  expect_plain_command(said_by(err), "bash install.sh")
  local_mocked_bindings(.installed_rtft_version = function()
    as.character(utils::packageVersion("rTFT")))
  local_mocked_bindings(.has_cmdstanr = function() FALSE)
  err <- expect_error(check_hpc_bundle(dir), class = "rTFT_needs_cmdstanr")
  expect_plain_command(said_by(err), "bash install.sh")
  local_mocked_bindings(.has_cmdstanr = function() TRUE, .cmdstan_found = function() FALSE)
  err <- expect_error(check_hpc_bundle(dir), class = "rTFT_no_cmdstan")
  expect_plain_command(said_by(err), "bash install.sh --job")
})


test_that("copies of one version are told apart by their commits, where both have one", {
  with_cmdstan()
  sha_a <- "3c690c9abc0000000000000000000000000000aa"
  sha_b <- "e15b42aef00000000000000000000000000000bb"
  local_mocked_bindings(.installed_rtft_commit = function(...) sha_a)
  dir <- write_quietly()
  expect_equal(readRDS(file.path(dir, "inputs", "settings.rds"))$rtft_commit, sha_a)
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
  expect_true(is.na(readRDS(file.path(plain, "inputs", "settings.rds"))$rtft_commit))
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
  expect_equal(.one_line(made_by_cli), "Only 1 of the 2 chains finished. See 'logs/'.")
  # base R quotes names with curly quotes where R runs in UTF-8
  quoted <- simpleError("there is no package called \u2018brms\u2019 \u201cnow\u201d")
  expect_equal(.one_line(quoted), "there is no package called 'brms' \"now\"")
  # the unit is chosen after rounding
  expect_equal(vapply(c(14, 119.6, 840, 7199, 41000), .duration, ""),
               c("14 s", "2 min", "14 min", "2 h", "11.4 h"))
})


test_that("the task numbers are written as Slurm's --array takes them, and read back", {
  expect_equal(.array_ranges(c(1, 2, 3, 5, 7, 8)), "1-3,5,7-8")
  expect_equal(.array_ranges(c(4, 2)), "2,4")
  expect_equal(.array_ranges(6), "6")
  expect_equal(.expand_ranges("1-3,5,7-8"), c(1:3, 5L, 7:8))
  expect_equal(.expand_ranges("6"), 6L)
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


test_that("a task leaves mc.cores as it found it", {
  # Set, mc.cores makes other functions copy the process with forks, which hung
  # for good on Nibi once the sampler's threads were running.
  seen <- NULL
  local_mocked_bindings(fit_hazard = function(...) {
    seen <<- getOption("mc.cores")
    fit_trials
  })
  dir <- trial_bundle()
  withr::local_options(mc.cores = 7)
  run_quietly(dir, "fit", 1)
  expect_identical(seen, 7)                  # not the number of chains, 2
  expect_identical(getOption("mc.cores"), 7)
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
  expect_identical(status_quietly(dir)$model, c("NA", "1", rep(c("NA", "1"), each = 3)))
  out <- suppressMessages(suppressWarnings(collect_hpc(dir)))
  expect_named(out$fits, c("NA", "1"))
  expect_named(out$kfold, c("NA", "1"))
  expect_identical(out$summary$model, c("NA", "1"))
})


test_that("a refit starts from its model's fit and leaves out the units of its fold", {
  refits <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(refits))
  dir <- trial_bundle()
  settings <- readRDS(file.path(dir, "inputs", "settings.rds"))
  tasks <- .read_tasks(dir, "kfold")
  # before its fit has ended, a refit cannot start
  err <- expect_error(.kfold_task(dir, settings, tasks[1, ], .bundle_id(dir)),
                      class = "rTFT_fit_not_ready")
  expect_plain_command(said_by(err), "bash submit.sh --fit 1 --kfold 1")
  err <- expect_error(run_quietly(dir, "kfold", 1), class = "rTFT_task_failed")
  expect_plain_command(said_by(err), "bash status.sh")
  expect_match(.read_record(dir, "kfold", 1)$message, "start from its fit", fixed = TRUE)
  # and status.sh says why as a sentence, with no bullet of cli left among its words
  said <- message_of(hpc_status(dir), "rTFT_hpc_status")
  expect_match(said, "refit 1 (exp_mem, fold 1): failed. The refits of", fixed = TRUE)
  expect_match(said, "not ended well. Run that task first;", fixed = TRUE)
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


test_that("a refit run again draws the same effects for the units it leaves out", {
  # a stand-in refit whose result is a draw of R's generator, as the effects of the
  # units left out are
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = function(...) list(elpd = c(u = stats::runif(1))))
  dir <- trial_bundle()
  settings <- readRDS(file.path(dir, "inputs", "settings.rds"))
  tasks <- .read_tasks(dir, "kfold")
  run_quietly(dir, "fit", 1)
  # each in a generator seeded at random, as a new R on the cluster is
  draw <- function(k, settings) {
    set.seed(NULL)
    .kfold_task(dir, settings, tasks[k, ], .bundle_id(dir))$elpd
  }
  expect_identical(draw(1, settings), draw(1, settings))
  expect_false(identical(draw(1, settings), draw(2, settings)))
  # without a seed, the generator is left alone
  settings$dots$seed <- NULL
  expect_false(identical(draw(1, settings), draw(1, settings)))
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
  writeLines("x", file.path(dir, "inputs", "extra.csv"))
  expect_error(run_hpc_task(dir, "fit", 1), class = "rTFT_bundle_changed")
  file.remove(file.path(dir, "inputs", "extra.csv"))
  expect_false(file.exists(.task_file(dir, "fit", 1, "started")))
  # another version stops the task before it fits, and the reason is recorded, so
  # that hpc_status() can say it
  local_mocked_bindings(.installed_rtft_version = function() "0.0.1")
  expect_error(run_quietly(dir, "fit", 1), class = "rTFT_task_failed")
  expect_length(record$calls, 0)
  ended <- .read_record(dir, "fit", 1)
  expect_equal(ended$status, "failed")
  expect_equal(ended$rtft_version, "0.0.1")
  expect_match(ended$message, "The installed rTFT is version \"0.0.1\"", fixed = TRUE)
  states <- status_quietly(dir)
  expect_match(states$reason[1], "The installed rTFT is version \"0.0.1\"", fixed = TRUE)
})


test_that("a task run again is seen as started, and its earlier end is kept aside", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  record$fail <- TRUE
  try(run_quietly(dir, "fit", 2), silent = TRUE)
  expect_equal(status_quietly(dir)$status[2], "failed")
  # while the second attempt runs (seen from inside the fit), the task is started,
  # and is not among those to submit again
  record$fail <- FALSE
  during <- seen <- said <- NULL
  local_mocked_bindings(fit_hazard = function(...) {
    during <<- .read_record(dir, "fit", 2)
    seen   <<- status_quietly(dir)
    said   <<- message_of(hpc_status(dir), "rTFT_hpc_status")
    fit_trials
  })
  run_quietly(dir, "fit", 2)
  expect_null(during)                                  # no record of an end yet
  expect_equal(seen$status[2], "started")
  expect_no_match(said, "--fit 1-2", fixed = TRUE)
  expect_match(said, "--fit 1 ", fixed = TRUE)
  # the earlier end is beside it, with its reason
  previous <- utils::read.csv(.task_file(dir, "fit", 2, "previous.csv"),
                              colClasses = "character")
  expect_equal(previous$status, "failed")
  expect_match(previous$message, "the sampler stopped", fixed = TRUE)
  expect_equal(.read_record(dir, "fit", 2)$status, "done")
})


# ===========================================================================
# 6. Following the run (hpc_status)
# ===========================================================================

test_that("hpc_status says that nothing has started, and how to start it", {
  dir <- trial_bundle()
  expect_message(states <- hpc_status(dir), class = "rTFT_hpc_status")
  expect_equal(nrow(states), 8)
  expect_equal(states$stage, rep(c("fit", "kfold"), c(2, 6)))
  expect_equal(states$task, c(1:2, 1:6))
  expect_equal(states$fold, c(NA, NA, rep(1:3, 2)))
  expect_equal(unique(states$status), "not started")
  said <- message_of(hpc_status(dir), "rTFT_hpc_status")
  expect_match(said, "Fits: 0 of 2 done; 2 not started.", fixed = TRUE)
  expect_match(said, "bash submit.sh --fit 1-2 --kfold 1-6", fixed = TRUE)
})


test_that("hpc_status gives each task's state, and why it is not done", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_quietly(dir, "fit", 1)
  record$fail <- TRUE
  try(run_quietly(dir, "fit", 2), silent = TRUE)
  run_quietly(dir, "kfold", 1)
  fake_start(dir, "kfold", 2, "80_2",
             "slurmstepd: error: *** JOB 82 ON nc1 CANCELLED AT 2026-10-17T22:00:00 DUE TO TIME LIMIT ***")
  fake_start(dir, "kfold", 3, "80_3",
             "slurmstepd: error: Detected 1 oom_kill event in StepId=83.batch.")
  fake_start(dir, "kfold", 4, "80_4", "Running task refit 4")
  fake_start(dir, "kfold", 5, "80_5")
  states <- status_quietly(dir)
  expect_equal(states$status, c("done", "failed", "done", "started", "started",
                                "started", "started", "not started"))
  expect_match(states$reason[2], "the sampler stopped", fixed = TRUE)
  expect_match(states$reason[4], "ran out of time", fixed = TRUE)
  expect_match(states$reason[5], "ran out of memory", fixed = TRUE)
  expect_match(states$reason[6], "may still be running", fixed = TRUE)
  expect_match(states$reason[6], "see logs/kfold-80_4.out", fixed = TRUE)
  expect_match(states$reason[7], "may still be running", fixed = TRUE)
  expect_true(all(is.na(states$reason[c(1, 3)])))
  expect_false(anyNA(states$seconds[1:3]))
  # what to submit again: what failed, never started, or was stopped by Slurm, but
  # not what may still be running
  said <- message_of(hpc_status(dir), "rTFT_hpc_status")
  expect_match(said, "bash submit.sh --fit 2 --kfold 2-3,6", fixed = TRUE)
  expect_match(said, "Refits for cross-validation: 1 of 6 done; 4 started, 1 not started.",
               fixed = TRUE)
  expect_plain_command(said, "the command sq lists")
})


test_that("a task stopped by Slurm is told apart by the line Slurm wrote in its log", {
  dir <- trial_bundle()
  dir.create(file.path(dir, "logs"), showWarnings = FALSE)
  reason <- function(line) {
    writeLines(c("Running task refit 1", line), file.path(dir, "logs", "kfold-90_1.out"))
    .log_reason(dir, "kfold", "90_1")
  }
  expect_match(reason("slurmstepd: error: *** JOB 91 ON nc1 CANCELLED AT 2026-10-19T01:00:00 DUE TO NODE FAILURE, SEE SLURMCTLD LOG FOR DETAILS ***"),
               "^it was stopped by Slurm, due to node failure; see logs/kfold-90_1.out$")
  expect_match(reason("slurmstepd: error: *** JOB 91 ON nc1 CANCELLED AT 2026-10-19T01:00:00 DUE TO PREEMPTION ***"),
               "due to preemption;", fixed = TRUE)
  # cancelled by hand, with scancel, Slurm gives no reason
  expect_match(reason("slurmstepd: error: *** JOB 91 ON nc1 CANCELLED AT 2026-10-19T01:00:00 ***"),
               "^it was cancelled; see")
  # a task stopped by Slurm is among those to submit again; one that may still be
  # running is not
  fake_start(dir, "kfold", 1, "90_1")
  reason("slurmstepd: error: *** JOB 91 ON nc1 CANCELLED AT 2026-10-19T01:00:00 DUE TO NODE FAILURE ***")
  expect_match(message_of(hpc_status(dir), "rTFT_hpc_status"), "--kfold 1-6", fixed = TRUE)
  reason("Running task refit 1, still sampling")
  expect_match(message_of(hpc_status(dir), "rTFT_hpc_status"), "--kfold 2-6", fixed = TRUE)
})


test_that("a task submitted again waits in the queue, and is not given to submit again", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle()
  run_quietly(dir, "fit", 1)
  record$fail <- TRUE
  try(run_quietly(dir, "fit", 2), silent = TRUE)
  # every file of the tasks written an hour ago; the record of submit.sh says fit 2
  # and every refit were submitted since, and fit 1 too (it stays done)
  files <- list.files(file.path(dir, "results"), full.names = TRUE)
  Sys.setFileTime(files, Sys.time() - 3600)
  now <- round(as.numeric(Sys.time()))
  writeLines(c('"time","stage","tasks","job"',
               paste0(now - 7200, ',"fit","1-2","90"'),
               paste0(now, ',"fit","1-2","95"'),
               paste0(now, ',"kfold","1-6","96"')),
             file.path(dir, "results", "submitted.csv"))
  sent <- .read_submissions(dir)
  expect_equal(nrow(sent), 10)
  expect_equal(sent$task[sent$job == "96"], 1:6)
  states <- status_quietly(dir)
  expect_equal(states$status, c("done", "queued", rep("queued", 6)))
  expect_match(states$reason[2], "as job 95, waiting to start", fixed = TRUE)
  said <- message_of(hpc_status(dir), "rTFT_hpc_status")
  expect_no_match(said, "To submit again", fixed = TRUE)
  expect_match(said, "bash submit.sh --fit 2 --kfold 1-6", fixed = TRUE)   # if gone
  expect_plain_command(said, "the command sq no longer lists")
  # a submission older than a task's files is an earlier one
  writeLines(c('"time","stage","tasks","job"', paste0(now - 7200, ',"fit","2","90"')),
             file.path(dir, "results", "submitted.csv"))
  states <- status_quietly(dir)
  expect_equal(states$status[1:3], c("done", "failed", "not started"))
  # and without a record of submissions nothing is queued
  file.remove(file.path(dir, "results", "submitted.csv"))
  expect_false(any(status_quietly(dir)$status == "queued"))
})


test_that("each command to submit again is printed on one line, whatever the width", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  dir <- trial_bundle()
  run_quietly(dir, "fit", 1)
  # cli wraps a message to the width in force when it is printed; testthat does not
  # wrap at all, so what is checked is the width in force as the message is given
  withr::local_options(cli.width = 20)
  width <- NULL
  said  <- NULL
  suppressMessages(withCallingHandlers(hpc_status(dir), rTFT_hpc_status = function(m) {
    width <<- getOption("cli.width")
    said  <<- conditionMessage(m)
  }))
  expect_identical(width, Inf)
  expect_match(said, "bash submit.sh --fit 2 --kfold 1-6", fixed = TRUE)
  expect_no_match(said, "`bash submit.sh", fixed = TRUE)     # plain text to paste
  expect_identical(getOption("cli.width"), 20)                 # and the width is back
})


test_that("hpc_status says when every task is done, and what comes next", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_all(dir)
  states <- status_quietly(dir)
  expect_equal(unique(states$status), "done")
  said <- message_of(hpc_status(dir), "rTFT_hpc_status")
  expect_match(said, "Every task is done", fixed = TRUE)
  # and no line without words, such as an empty list of tasks would leave
  expect_true(all(grepl("[[:alpha:]]", strsplit(said, "\n")[[1]])))
  expect_match(said, "collect_hpc()", fixed = TRUE)
  expect_no_match(said, "submit.sh", fixed = TRUE)
})


test_that("a result written for another bundle, or that went missing, is not counted done", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  dir   <- trial_bundle()
  other <- trial_bundle(kfold = FALSE, seed = 2)
  run_quietly(other, "fit", 1)
  file.copy(file.path(other, "results", c("fit-001.csv", "fit-001.rds")),
            file.path(dir, "results"))
  run_quietly(dir, "fit", 2)
  file.remove(.task_file(dir, "fit", 2, "rds"))
  states <- status_quietly(dir)
  expect_equal(states$status[1:2], c("failed", "failed"))
  expect_match(states$reason[1], "another bundle", fixed = TRUE)
  expect_match(states$reason[2], "missing", fixed = TRUE)
})


# ===========================================================================
# 7. Collecting the results (collect_hpc)
# ===========================================================================

test_that("collect_hpc returns what fit_hazard_models() returns for the same run", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_all(dir)
  expect_message(
    out <- suppressWarnings(collect_hpc(dir)), class = "rTFT_hpc_collected")
  set.seed(3)
  ref <- suppressMessages(suppressWarnings(fit_hazard_models(
    trial_models, long_trials, kfold = TRUE, K = 3, chains = 2, seed = 1,
    sd_prior = "exponential(1)")))
  expect_named(out, c("fits", "loo", "kfold", "summary"))
  expect_identical(names(out$fits), c("exp_mem", "exp_weibull"))
  expect_equal(out$fits$exp_mem, ref$fits$exp_mem)
  expect_equal(out$loo, ref$loo)
  # the same summary, with the round each model comes from: here, the first bundle
  expect_equal(names(out$summary), append(names(ref$summary), "round", after = 1))
  expect_equal(out$summary$round, c(1L, 1L))
  expect_equal(dplyr::select(out$summary, -"seconds", -"round"),
               dplyr::select(ref$summary, -"seconds"))
  expect_false(anyNA(out$summary$seconds))
  # the refits are scored on the folds of the bundle, unit by unit
  settings <- readRDS(file.path(dir, "inputs", "settings.rds"))
  for (m in names(out$kfold)) {
    expect_s3_class(out$kfold[[m]], c("kfold", "loo"))
    expect_equal(out$kfold[[m]]$estimates, ref$kfold[[m]]$estimates)
    expect_equal(out$kfold[[m]]$pointwise, ref$kfold[[m]]$pointwise)
    expect_equal(out$kfold[[m]]$folds, settings$folds)
    expect_equal(out$kfold[[m]]$K, 3)
  }
  expect_no_error(suppressMessages(hazard_weights(out$kfold)))
})


test_that("models and folds without a result are named, and the rest is collected", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_quietly(dir, "fit", 1)
  for (i in 1:2) run_quietly(dir, "kfold", i)
  expect_warning(expect_warning(
    out <- suppressMessages(suppressWarnings(
      collect_hpc(dir), classes = c("rTFT_not_converged", "rTFT_loo_unreliable"))),
    class = "rTFT_tasks_missing"), class = "rTFT_models_failed")
  expect_equal(out$summary$status, c("fitted", "missing"))
  expect_named(out$fits, "exp_mem")
  expect_length(out$kfold, 0)
  expect_match(out$summary$message[1], "fold 3 has no result", fixed = TRUE)
  expect_match(out$summary$message[2], "No result", fixed = TRUE)
  expect_true(all(is.na(out$summary$elpd_kfold)))
  said <- message_of(collect_hpc(dir), "rTFT_tasks_missing")
  expect_match(said, "1 of the 2 models has no result yet", fixed = TRUE)
  expect_plain_command(said, "bash status.sh")
})


test_that("a fit that failed is reported with its reason, and its refits are not sought", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_quietly(dir, "fit", 1)
  for (i in 1:3) run_quietly(dir, "kfold", i)
  record$fail <- TRUE
  try(run_quietly(dir, "fit", 2), silent = TRUE)
  out <- suppressMessages(suppressWarnings(collect_hpc(dir)))
  expect_equal(out$summary$status, c("fitted", "failed"))
  expect_match(out$summary$message[2], "the sampler stopped", fixed = TRUE)
  expect_false(is.na(out$summary$seconds[2]))
  expect_named(out$kfold, "exp_mem")
})


test_that("a model that failed is told apart from one that has no result yet", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  three <- build_hazard_models("exp", c("mem", "weibull", "llogis"),
                               dpar_formulas = list(loga = ~ 1 + (1 | trial)))
  dir <- trial_bundle(kfold = FALSE, models = three)
  run_quietly(dir, "fit", 1)
  record$fail <- TRUE
  try(run_quietly(dir, "fit", 2), silent = TRUE)
  # the third has not run: it is missing, not failed, and the others are not all
  # said to have been fitted
  failed <- message_of(collect_hpc(dir), "rTFT_models_failed")
  expect_match(failed, "1 of the 3 models could not be fitted: \"exp_weibull\"")
  expect_no_match(failed, "The other models were fitted", fixed = TRUE)
  expect_match(message_of(collect_hpc(dir), "rTFT_tasks_missing"),
               "1 of the 3 models has no result yet: \"exp_llogis\"")
  # once it has run, every model is fitted or failed, and that is said
  record$fail <- FALSE
  run_quietly(dir, "fit", 3)
  expect_match(message_of(collect_hpc(dir), "rTFT_models_failed"),
               "The other models were fitted", fixed = TRUE)
})


test_that("a refit whose score is undefined is reported, and one of -Inf is kept", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_all(dir)
  # model 1, fold 2: an undefined score
  file <- .task_file(dir, "kfold", 2, "rds")
  res  <- readRDS(file)
  res$elpd[1] <- NaN
  saveRDS(res, file)
  # model 2, fold 2: a unit given no chance at all
  file <- .task_file(dir, "kfold", 5, "rds")
  res  <- readRDS(file)
  res$elpd[1] <- -Inf
  saveRDS(res, file)
  out <- suppressMessages(suppressWarnings(collect_hpc(dir)))
  expect_named(out$kfold, "exp_weibull")
  expect_match(out$summary$message[1],
               paste("undefined for", names(res$elpd)[1]), fixed = TRUE)
  expect_true(is.na(out$summary$elpd_kfold[1]))
  # kept as kfold_hazard() keeps it: the estimate is -Inf
  expect_equal(out$kfold$exp_weibull$estimates["elpd_kfold", "Estimate"], -Inf)
  expect_equal(out$summary$elpd_kfold[2], -Inf)
})


test_that("a result or a record cut short in the copy back is named, and not read", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  dir <- trial_bundle(kfold = FALSE)
  for (i in 1:2) run_quietly(dir, "fit", i)
  cut_in_half <- function(file) {
    whole <- readBin(file, "raw", file.size(file))
    writeBin(whole[seq_len(length(whole) %/% 2)], file)
  }
  cut_in_half(.task_file(dir, "fit", 2, "rds"))
  expect_error(suppressWarnings(collect_hpc(dir)), class = "rTFT_result_unreadable")
  expect_match(tryCatch(suppressWarnings(collect_hpc(dir)), error = conditionMessage),
               "fit-002.rds", fixed = TRUE)
  # a record cut short is not taken for one written by another bundle
  file.copy(.task_file(dir, "fit", 1, "csv"), .task_file(dir, "fit", 2, "csv"),
            overwrite = TRUE)
  cut_in_half(.task_file(dir, "fit", 1, "csv"))
  expect_error(suppressWarnings(collect_hpc(dir)), class = "rTFT_result_unreadable")
  expect_match(tryCatch(suppressWarnings(collect_hpc(dir)), error = conditionMessage),
               "fit-001.csv", fixed = TRUE)
  writeBin(raw(0), .task_file(dir, "fit", 1, "csv"))                    # and empty
  expect_error(suppressWarnings(collect_hpc(dir)), class = "rTFT_result_unreadable")
})


test_that("collecting with another rTFT than the one that fitted says so", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  sha_a <- "3c690c9abc0000000000000000000000000000aa"
  local_mocked_bindings(.installed_rtft_commit = function(...) sha_a)
  dir <- trial_bundle(kfold = FALSE)
  for (i in 1:2) run_quietly(dir, "fit", i)
  quiet <- function(expr) suppressWarnings(expr)
  # the same rTFT: nothing is said
  expect_no_message(quiet(collect_hpc(dir)), class = "rTFT_collected_by_other_version")
  # another commit of the same version
  local_mocked_bindings(.installed_rtft_commit = function(...) "e15b42aef0")
  said <- message_of(collect_hpc(dir), "rTFT_collected_by_other_version")
  expect_match(said, "made by rTFT .* \\(commit 3c690c9\\); this is rTFT .* \\(commit e15b42a\\)")
  # another version, installed from a folder (no commit): it still collects
  local_mocked_bindings(.installed_rtft_version = function() "9.9.9",
                        .installed_rtft_commit = function(...) NA_character_)
  expect_message(out <- quiet(collect_hpc(dir)), class = "rTFT_collected_by_other_version")
  expect_named(out$fits, c("exp_mem", "exp_weibull"))
})


test_that("results are never assembled with a bundle other than their own", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  dir   <- trial_bundle(kfold = FALSE)
  other <- trial_bundle(kfold = FALSE, seed = 2)
  for (i in 1:2) run_quietly(other, "fit", i)
  run_quietly(dir, "fit", 1)
  file.copy(file.path(other, "results", c("fit-002.csv", "fit-002.rds")),
            file.path(dir, "results"))
  expect_error(suppressWarnings(collect_hpc(dir)), class = "rTFT_result_not_of_bundle")
  # nor with a bundle whose files changed after the run
  writeLines("x", file.path(other, "inputs", "extra.csv"))
  expect_error(collect_hpc(other), class = "rTFT_bundle_changed")
})


test_that("collect_hpc takes each model from its latest round", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()),
                        .refit_fold = stand_in_refit(new_record()))
  dir <- trial_bundle()
  run_all(dir)
  two <- round_quietly(dir, "exp_weibull", iter = 4000)
  # the refits of the first bundle for exp_weibull are taken out: those of its
  # round must be the ones collected
  of_weibull <- .read_tasks(dir, "kfold")$task[.read_tasks(dir, "kfold")$model ==
                                                  "exp_weibull"]
  for (j in of_weibull)
    file.remove(.task_file(dir, "kfold", j, c("csv", "rds")))
  for (stage in c("fit", "kfold"))
    for (j in .read_tasks(two, stage)$task) run_quietly(two, stage, j)
  expect_message(out <- suppressWarnings(collect_hpc(c(dir, two))),
                 "from 2 rounds", class = "rTFT_hpc_collected")
  expect_equal(out$summary$model, c("exp_mem", "exp_weibull"))
  expect_equal(out$summary$round, c(1L, 2L))
  expect_equal(out$summary$status, c("fitted", "fitted"))
  expect_named(out$kfold, c("exp_mem", "exp_weibull"))
  expect_false(anyNA(out$summary$elpd_kfold))
  # the folds of every round are the first bundle's, so the models can be weighted
  expect_equal(out$kfold$exp_weibull$folds, out$kfold$exp_mem$folds)
  expect_no_error(suppressMessages(hazard_weights(out$kfold)))
  # in any order
  again <- suppressMessages(suppressWarnings(collect_hpc(c(two, dir))))
  expect_equal(dplyr::select(again$summary, -"seconds"),
               dplyr::select(out$summary, -"seconds"))
  # a round alone holds only its own models
  alone <- suppressMessages(suppressWarnings(collect_hpc(two)))
  expect_equal(alone$summary$model, "exp_weibull")
  expect_equal(alone$summary$round, 2L)
})


test_that("a model left out of a round is collected from the round that fitted it last", {
  local_mocked_bindings(fit_hazard = stand_in_fit(new_record()))
  dir <- trial_bundle(kfold = FALSE)
  for (i in 1:2) run_quietly(dir, "fit", i)
  two <- round_quietly(dir, "exp_weibull", iter = 4000)
  run_quietly(two, "fit", 1)
  three <- round_quietly(two, "exp_mem", iter = 4000)
  run_quietly(three, "fit", 1)
  out <- suppressMessages(suppressWarnings(collect_hpc(c(dir, two, three))))
  expect_equal(out$summary$model, c("exp_mem", "exp_weibull"))
  expect_equal(out$summary$round, c(3L, 2L))
  expect_equal(out$summary$status, c("fitted", "fitted"))
  # a round collected alone holds only the model it fits, though its table holds both
  alone <- suppressMessages(suppressWarnings(collect_hpc(three)))
  expect_equal(alone$summary$model, "exp_mem")
  expect_equal(alone$summary$round, 3L)
})


test_that("a model's latest round counts even when it failed or has no result", {
  record <- new_record()
  local_mocked_bindings(fit_hazard = stand_in_fit(record))
  dir <- trial_bundle(kfold = FALSE)
  for (i in 1:2) run_quietly(dir, "fit", i)
  two <- round_quietly(dir, "exp_weibull", iter = 4000)
  # the round has not run: the model is missing, not taken from the first bundle
  expect_warning(
    out <- suppressMessages(suppressWarnings(collect_hpc(c(dir, two)),
                                             classes = "rTFT_models_failed")),
    class = "rTFT_tasks_missing")
  expect_equal(out$summary$status, c("fitted", "missing"))
  expect_equal(out$summary$round, c(1L, 2L))
  expect_named(out$fits, "exp_mem")
  # the round failed: the model failed, with its reason
  record$fail <- TRUE
  expect_error(run_quietly(two, "fit", 1), class = "rTFT_task_failed")
  out <- suppressMessages(suppressWarnings(collect_hpc(c(dir, two))))
  expect_equal(out$summary$status, c("fitted", "failed"))
  expect_match(out$summary$message[2], "the sampler stopped", fixed = TRUE)
  expect_named(out$fits, "exp_mem")
})


test_that("collect_hpc refuses folders that are not a bundle and its rounds", {
  dir   <- trial_bundle(kfold = FALSE)
  two   <- round_quietly(dir, "exp_weibull")
  three <- round_quietly(two, "exp_weibull")
  other <- trial_bundle(kfold = FALSE, seed = 2)
  of_other <- round_quietly(other, "exp_mem")
  also_two <- round_quietly(dir, "exp_mem")
  # refused, with the reason given in the message
  not_chain <- function(dirs, why) {
    said <- tryCatch(suppressWarnings(collect_hpc(dirs)),
                     rTFT_rounds_not_a_chain = function(e) conditionMessage(e))
    expect_match(said, why, fixed = TRUE)
  }
  not_chain(c(dir, dir), "given twice")
  not_chain(c(two, three), "first bundle is not among them")
  not_chain(c(dir, three), "the bundle round 3 follows is not among them")
  not_chain(c(dir, two, of_other), "more than one first bundle")
  not_chain(c(dir, other), "more than one first bundle")
  not_chain(c(dir, two, also_two), "two of them are the same round")
  # a round whose data were changed, with its checksums written again
  saveRDS(long_trials[-1, ], file.path(two, "inputs", "data.rds"))
  sums <- .bundle_checksums(two)
  .write_lines_lf(paste(sums, names(sums), sep = "  "),
                  file.path(two, "inputs", "checksums.txt"))
  settings <- readRDS(file.path(three, "inputs", "settings.rds"))
  settings$round$previous <- .bundle_id(two)
  saveRDS(settings, file.path(three, "inputs", "settings.rds"))
  sums <- .bundle_checksums(three)
  .write_lines_lf(paste(sums, names(sums), sep = "  "),
                  file.path(three, "inputs", "checksums.txt"))
  not_chain(c(dir, two, three), "do not hold the same data")
  # and what is not a folder
  for (bad in list(character(0), NA_character_, 1))
    expect_error(collect_hpc(bad), class = "rTFT_bad_dir")
})


# ===========================================================================
# 8. The scripts that submit and run the tasks
# ===========================================================================

# Put stand-ins for the commands of a cluster first on the path: `sbatch`, which
# records how it was called and prints a job number as Slurm prints it with
# --parsable; `module`, which does nothing; and `Rscript`, which records how it was
# called. Return the file of the records.
fake_cluster <- function(env = parent.frame()) {
  bin <- withr::local_tempdir(.local_envir = env)
  log <- file.path(bin, "calls.txt")
  writeLines(c("#!/bin/bash", 'echo "sbatch $*" >> "$CALLS"',
               'echo "$(( 1000 + $(wc -l < "$CALLS") ));nibi"'),
             file.path(bin, "sbatch"))
  writeLines(c("#!/bin/bash", "exit 0"), file.path(bin, "module"))
  writeLines(c("#!/bin/bash", 'echo "Rscript $*" >> "$CALLS"',
               'echo "TMPDIR=$TMPDIR" >> "$CALLS"'), file.path(bin, "Rscript"))
  Sys.chmod(file.path(bin, c("sbatch", "module", "Rscript")), "0755")
  file.create(log)
  withr::local_envvar(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"), CALLS = log,
                      .local_envir = env)
  log
}

# Run a script of the bundle in `dir` with bash, with `args`, and return what it
# printed, with its exit status as the attribute "status" when it is not 0.
run_script <- function(dir, script, args = character(0)) {
  suppressWarnings(system2("bash", shQuote(c(file.path(dir, script), args)),
                           stdout = TRUE, stderr = TRUE))
}


test_that("submit.sh submits the fits, then the refits to start when the fits end", {
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  calls <- fake_cluster()
  dir <- trial_bundle()
  out <- run_script(dir, "submit.sh")
  expect_null(attr(out, "status"))
  sent <- readLines(calls)
  expect_length(sent, 2)
  expect_match(sent[1], "--account=def-test --job-name=rtft-fit --array=1-2 --time=12:00:00 --mem-per-cpu=4000M --cpus-per-task=2 --output=logs/fit-%A_%a.out task.sh fit",
               fixed = TRUE)
  expect_match(sent[2], "--array=1-6", fixed = TRUE)
  expect_match(sent[2], "--dependency=afterany:1001 task.sh kfold", fixed = TRUE)
  expect_match(paste(out, collapse = "\n"), "as job 1001.", fixed = TRUE)
  # each submission is recorded, for hpc_status() to read
  sent <- utils::read.csv(file.path(dir, "results", "submitted.csv"),
                          colClasses = "character")
  expect_equal(sent$stage, c("fit", "kfold"))
  expect_equal(sent$tasks, c("1-2", "1-6"))
  expect_equal(sent$job, c("1001", "1002"))
  expect_true(all(abs(as.numeric(sent$time) - as.numeric(Sys.time())) < 120))
  expect_equal(unique(status_quietly(dir)$status), "queued")
  expect_true(all(dir.exists(file.path(dir, c("logs", "results", "stan")))))
})


test_that("submit.sh submits only the tasks it is given, and refuses what it cannot submit", {
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  calls <- fake_cluster()
  dir <- trial_bundle()
  run_script(dir, "submit.sh", c("--kfold", "2-3,5"))
  sent <- readLines(calls)
  expect_length(sent, 1)
  expect_match(sent, "--array=2-3,5", fixed = TRUE)
  expect_no_match(sent, "dependency", fixed = TRUE)
  run_script(dir, "submit.sh", c("--fit", "2", "--kfold", "4"))
  sent <- readLines(calls)[2:3]
  expect_match(sent[1], "--array=2 ", fixed = TRUE)
  expect_match(sent[2], "--dependency=afterany:1002", fixed = TRUE)
  # what cannot be submitted is not
  for (bad in list("--fit", c("--fit", "x"), c("--fit", "1;2"), c("--all", "1"),
                   c("--fit", "0"), c("--fit", "1-3"), c("--kfold", "2,7")))
    expect_equal(attr(run_script(dir, "submit.sh", bad), "status"), 1L)
  expect_match(paste(run_script(dir, "submit.sh", c("--fit", "1-3")), collapse = " "),
               "This bundle has 2 fits; fit 3 is not one of them.", fixed = TRUE)
  expect_length(readLines(calls), 3)
  # a bundle without refits submits its fits alone, and has no refits to submit
  plain <- trial_bundle(kfold = FALSE)
  run_script(plain, "submit.sh")
  expect_length(readLines(calls), 4)
  expect_match(readLines(calls)[4], "task.sh fit", fixed = TRUE)
  expect_equal(attr(run_script(plain, "submit.sh", c("--kfold", "1")), "status"), 1L)
})


test_that("task.sh loads R and runs the task of its stage", {
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  calls <- fake_cluster()
  dir <- trial_bundle()
  node_disk <- withr::local_tempdir()
  withr::local_envvar(SLURM_SUBMIT_DIR = dir, SLURM_ARRAY_TASK_ID = "2",
                      SLURM_TMPDIR = node_disk)
  out <- run_script(dir, "task.sh", "kfold")
  expect_null(attr(out, "status"))
  # R is run with its temporary files on the node's own disk
  expect_equal(readLines(calls), c("Rscript -e rTFT::run_hpc_task('.', 'kfold')",
                                   paste0("TMPDIR=", node_disk)))
  expect_equal(attr(run_script(dir, "task.sh", "all"), "status"), 1L)
  expect_length(readLines(calls), 2)
})


test_that("install.sh keeps its temporary files out of the login node's /tmp", {
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  calls <- fake_cluster()
  dir <- trial_bundle()
  tmp_of <- function(line) normalizePath(sub("^TMPDIR=", "", line), mustWork = FALSE)
  # on the login node, in a folder of the bundle
  withr::local_envvar(SLURM_JOB_NAME = NA, SLURM_TMPDIR = NA)
  expect_null(attr(run_script(dir, "install.sh"), "status"))
  expect_equal(tmp_of(readLines(calls)[2]), normalizePath(file.path(dir, "tmp")))
  expect_true(dir.exists(file.path(dir, "tmp")))
  # in a Slurm job, on the node's own disk
  node_disk <- withr::local_tempdir()
  withr::local_envvar(SLURM_TMPDIR = node_disk)
  run_script(dir, "install.sh")
  expect_equal(tmp_of(readLines(calls)[4]), normalizePath(node_disk))
})


test_that("status.sh loads R and reports on the bundle", {
  skip_if(!bash_reads_files(), "this machine's bash cannot read the files R writes")
  calls <- fake_cluster()
  dir <- trial_bundle()
  out <- run_script(dir, "status.sh")
  expect_null(attr(out, "status"))
  expect_equal(readLines(calls)[1], "Rscript -e invisible(rTFT::hpc_status(\".\"))")
})
