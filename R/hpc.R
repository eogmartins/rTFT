# hpc.R

# Fitting a table of models on a cluster that schedules its jobs with Slurm. A run
# that would take hours or days on one machine is written into a folder, the
# bundle, which holds the data, the models and the settings of the run with their
# checksums, and the short scripts that load R and call the functions of this
# package on the cluster. The bundle is copied to the cluster, checked there before
# anything is submitted, and its results are brought back and assembled into what
# fit_hazard_models() returns. Every value the scripts need, the account above all,
# is written once, in the bundle's settings.sh.

# ===========================================================================
# 1. Helpers
# ===========================================================================

# .has_cmdstanr() / .cmdstan_found() / .on_windows() / .installed_rtft_version()
# .installed_rtft_commit() / .rtft_source()

# Whether the cmdstanr package is installed. It is the interface to CmdStan that the
# cluster's fits use, and a suggested package, since it is not on CRAN.
.has_cmdstanr <- function() requireNamespace("cmdstanr", quietly = TRUE)

# Whether cmdstanr can find a CmdStan installation, the compiler and libraries that
# turn a model into a program. FALSE when cmdstanr has no path to it, or the path
# leads nowhere.
.cmdstan_found <- function() {
  path <- tryCatch(suppressWarnings(suppressMessages(cmdstanr::cmdstan_path())),
                   error = function(e) NULL)
  !is.null(path) && dir.exists(path)
}

# Whether this machine runs Windows, which has ssh and scp but no rsync.
.on_windows <- function() .Platform$OS.type == "windows"

# The version of rTFT that is installed, as a character string.
.installed_rtft_version <- function() as.character(utils::packageVersion("rTFT"))

# The commit the installed rTFT was installed from, read from `d`, its description: the
# full identifier GitHub gives it, when it came from GitHub, and NA otherwise (a
# release from an archive, or a copy on this machine).
.installed_rtft_commit <- function(d = utils::packageDescription("rTFT")) {
  if (identical(d$RemoteType, "github") && !is.null(d$RemoteSha)) return(d$RemoteSha)
  NA_character_
}

# Where the cluster can install the rTFT that is installed here, as
# remotes::install_github() takes it, read from `d`, its description, or NA when that
# cannot be told: the commit it was installed from, when it came from GitHub; the
# tag of its release, for a released version; and NA for a development version
# installed from a copy on this machine, whose code may be on no branch at all.
.rtft_source <- function(d = utils::packageDescription("rTFT")) {
  if (identical(d$RemoteType, "github") && !is.null(d$RemoteSha))
    return(paste0(d$RemoteUsername, "/", d$RemoteRepo, "@", d$RemoteSha))
  if (length(strsplit(d$Version, ".", fixed = TRUE)[[1]]) == 3L)
    return(paste0("eogmartins/rTFT@v", d$Version))
  NA_character_
}


# .check_hpc_settings() / .template_dir() / .fill_template() -------------------

# Check the settings of the cluster given to write_hpc_bundle(): `account` and
# `time` as Slurm takes them, `mem_per_cpu` as an amount with its unit, `email`,
# `r_library` and `remote` as single strings or NULL (`remote` as a host and a
# folder, without spaces, which the copy command could not carry), and `modules`
# as a character vector of module names. Return TRUE invisibly.
.check_hpc_settings <- function(account, time, mem_per_cpu, email, modules, r_library,
                                remote) {
  one <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  bad <- c(
    account     = !one(account),
    time        = !one(time) ||
                    !grepl("^([0-9]+-)?[0-9]+:[0-5][0-9]:[0-5][0-9]$", time),
    mem_per_cpu = !one(mem_per_cpu) || !grepl("^[0-9]+[KMGT]?$", mem_per_cpu),
    email       = !is.null(email) && (!one(email) || !grepl("@", email, fixed = TRUE)),
    modules     = !is.character(modules) || !length(modules) || anyNA(modules) ||
                    any(!nzchar(modules)),
    r_library   = !is.null(r_library) && !one(r_library),
    remote      = !is.null(remote) &&
                    (!one(remote) || !grepl(":", remote, fixed = TRUE) ||
                       grepl("[[:space:]]", remote)))
  if (any(bad))
    cli::cli_abort(c(
      "Some settings of the cluster cannot be used: {.arg {names(bad)[bad]}}.",
      "i" = "{.arg account} is the name of your allocation (such as {.val def-someone}), {.arg time} is written as {.val 12:00:00} or {.val 1-00:00:00}, and {.arg mem_per_cpu} as {.val 4000M} or {.val 4G}.",
      "i" = "{.arg remote} is the place on the cluster to copy the bundle to, as {.val user@nibi.alliancecan.ca:scratch/my-run}."
    ), class = "rTFT_bad_hpc_setting")
  invisible(TRUE)
}

# The files of a template, which write_hpc_bundle() fills in and writes into the
# bundle.
.script_files <- c("settings.sh", "install.sh", "install.R", "check.sh", "submit.sh",
                   "task.sh", "status.sh", "README.txt")

# The folder that holds the scripts of `template`: "alliance" for the one written
# for the clusters of the Digital Research Alliance of Canada, or the path of a
# folder of your own holding the same files.
.template_dir <- function(template) {
  dir <- if (identical(template, "alliance"))
    system.file("templates", "alliance", package = "rTFT") else template
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir) ||
      !all(file.exists(file.path(dir, .script_files))))
    cli::cli_abort(c(
      "{.arg template} must be {.val alliance}, or a folder that holds {.file {(.script_files)}}.",
      "x" = "You supplied {.val {template}}."
    ), class = "rTFT_bad_template")
  dir
}

# Fill the placeholders of a template's `text`, written as {{name}}, with the named
# `values`, and return the text. A placeholder left unfilled is an error, so a
# script never reaches the cluster with one in it.
.fill_template <- function(text, values) {
  for (name in names(values))
    text <- gsub(paste0("{{", name, "}}"), values[[name]], text, fixed = TRUE)
  left <- unique(unlist(regmatches(text, gregexpr("\\{\\{[A-Za-z_]+\\}\\}", text))))
  if (length(left))
    cli::cli_abort(c(
      "A script of the template has placeholders that were not filled: {.val {left}}.",
      "i" = "The placeholders a template may use are {.val {paste0('{{', names(values), '}}')}}."
    ), class = "rTFT_template_unfilled")
  text
}


# .write_lines_lf() / .write_csv_lf() ------------------------------------------

# Write the lines `text` to the file `path` with line feeds alone, whatever the
# system. Written on Windows in the usual way, each line would also end with a
# carriage return, which bash on the cluster reads as part of the line.
.write_lines_lf <- function(text, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(text, con, sep = "\n")
}

# Write the data frame `x` to the file `path` as a table of comma-separated values
# without row names, with line feeds alone, whatever the system.
.write_csv_lf <- function(x, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  utils::write.csv(x, con, row.names = FALSE)
}


# .copy_command() / .copy_back_command() ---------------------------------------

# The command that copies the bundle at `local_path` to the folder `remote_folder` of
# the cluster `ssh_target`, as lines to write into README.txt, and a note on it,
# returned as a list. Elsewhere it is one rsync command, split over lines a shell
# reads as one, which makes the whole path of the folder on the cluster's end of the
# same connection, since rsync makes only the last folder of its destination. On
# Windows, which has ssh and scp but no rsync, with `windows`, it is an ssh command
# that makes the folders above the destination and an scp command that copies the
# folder there; scp names the copy after the destination when it does not exist yet.
.copy_command <- function(local_path, ssh_target, remote_folder, windows) {
  if (!windows)
    return(list(
      command = c(
        sprintf('rsync -av --rsync-path="mkdir -p %s && rsync" \\', remote_folder),
        sprintf('  "%s/" \\', local_path),
        sprintf("  %s:%s/", ssh_target, remote_folder)),
      note = paste(
        "It is one command over three lines; copy all three. The --rsync-path part",
        "makes the folder on the cluster first, with any folders above it that do",
        "not exist yet; rsync alone makes only the last.")))
  parent <- dirname(remote_folder)
  list(
    command = c(
      sprintf('ssh %s "mkdir -p %s"', ssh_target, parent),
      sprintf('scp -r "%s" %s:%s', local_path, ssh_target, remote_folder)),
    note = paste(
      "Two commands, for the PowerShell of Windows: the first makes the folders",
      "above the one on the cluster, the second copies the bundle into it, which",
      "must not exist yet. Graphical tools do the same: Globus, WinSCP or MobaXterm,",
      "which the Alliance documents."))
}


# The command that copies the folder results/ of the bundle from the folder
# `remote_folder` of the cluster `ssh_target` back into the bundle at `local_path`,
# as lines to write into README.txt: one rsync command over three lines, or on
# Windows, with `windows`, one scp command, which copies the folder into the bundle.
.copy_back_command <- function(local_path, ssh_target, remote_folder, windows) {
  if (!windows)
    return(c("rsync -av \\",
             sprintf("  %s:%s/results/ \\", ssh_target, remote_folder),
             sprintf('  "%s/results/"', local_path)))
  sprintf('scp -r %s:%s/results "%s"', ssh_target, remote_folder, local_path)
}


# .check_rtft_version() --------------------------------------------------------

# Check that the installed rTFT is the version that wrote the bundle, as recorded in
# `settings`, the bundle's settings, and return the installed version. Copies of a
# development version share a number and can differ in their code, so where both
# the installed rTFT and the one that wrote the bundle came from GitHub, their
# commits must be the same as well.
.check_rtft_version <- function(settings) {
  installed <- .installed_rtft_version()
  fix <- "Run {.code bash install.sh}, which installs the version that wrote the bundle, so that the models fitted here are the ones written there."
  if (!identical(installed, settings$rtft_version))
    cli::cli_abort(c(
      "The installed rTFT is version {.val {installed}}, and the bundle was written by version {.val {settings$rtft_version}}.",
      "i" = fix
    ), class = "rTFT_version_differs")
  # a bundle written before the commit was recorded has none, and is not checked
  # by it
  wrote <- if (is.null(settings$rtft_commit)) NA_character_ else settings$rtft_commit
  have  <- .installed_rtft_commit()
  if (!is.na(wrote) && !is.na(have) && !identical(have, wrote))
    cli::cli_abort(c(
      "The installed rTFT {installed} is from commit {.val {substr(have, 1, 7)}}, and the bundle was written by the rTFT {installed} of commit {.val {substr(wrote, 1, 7)}}.",
      "i" = fix
    ), class = "rTFT_version_differs")
  installed
}


# .model_key() / .check_distinct_models() ---------------------------------------

# The model in row `i` of `models` as one string, to compare models by: its
# failure-rate function, its baseline, the formula of each parameter, by parameter,
# and its co-stressors, everything but its name.
.model_key <- function(models, i) {
  as_text <- function(x) paste(deparse(x), collapse = " ")
  formulas <- models$dpar_formulas[[i]]
  given    <- sort(names(formulas))
  paste(c(models$m_form[i], models$baseline[i],
          paste(given, vapply(formulas[given], as_text, character(1))),
          as_text(models$costressors[[i]])),
        collapse = "\r")
}

# Check that no two models of `models` are the same model under different names, and
# return TRUE invisibly. On the cluster their tasks would run at the same time and
# compile one program, whose file is named by its code, into the same file.
.check_distinct_models <- function(models) {
  keys <- vapply(seq_len(nrow(models)), \(i) .model_key(models, i), character(1))
  same <- keys %in% keys[duplicated(keys)]
  if (any(same)) {
    sets <- vapply(split(models$model[same], keys[same]),
                   \(m) paste0('"', m, '"', collapse = " and "), character(1))
    cli::cli_abort(c(
      "Some models of the table are the same model under different names: {sets}.",
      "i" = "Fitted at the same time on the cluster, they would compile one program into the same file. Keep one of each."
    ), class = "rTFT_duplicate_models")
  }
  invisible(TRUE)
}


# .bundle_files() / .bundle_checksums() / .verify_bundle() / .bundle_id() ------

# The files of the folder `bundle/` of the bundle in `dir` that are checked for
# changes, with paths relative to `dir`.
.bundle_files <- function(dir) {
  files <- setdiff(list.files(file.path(dir, "bundle")), "checksums.txt")
  file.path("bundle", sort(files))
}

# The MD5 checksum of each of those files, named by its path within `dir`.
.bundle_checksums <- function(dir) {
  files <- .bundle_files(dir)
  stats::setNames(unname(tools::md5sum(file.path(dir, files))), files)
}

# Check that `dir` holds a complete bundle whose files are all as they were written,
# by their checksums, and return the checksums of the files found, invisibly. A file
# missing, changed since or added since is an error, so a copy that was cut short or
# a file edited by hand is found before anything is fitted or assembled.
.verify_bundle <- function(dir) {
  needed <- file.path("bundle", c("data.rds", "models.rds", "settings.rds",
                                  "tasks-fit.csv", "checksums.txt"))
  absent <- needed[!file.exists(file.path(dir, needed))]
  if (length(absent))
    cli::cli_abort(c(
      "{.path {dir}} is not a complete bundle: {.file {absent}} {?is/are} missing.",
      "i" = "Copy the whole folder written by {.fn write_hpc_bundle} again."
    ), class = "rTFT_bad_bundle")
  recorded <- utils::read.table(file.path(dir, "bundle", "checksums.txt"),
                                col.names = c("md5", "file"), colClasses = "character")
  found <- .bundle_checksums(dir)
  now   <- unname(found[recorded$file])          # missing for a file no longer there
  changed <- c(recorded$file[is.na(now) | now != recorded$md5],
               setdiff(names(found), recorded$file))       # and a file added since
  if (length(changed))
    cli::cli_abort(c(
      "{cli::qty(length(changed))}{?This file/These files} of the bundle {?is/are} not as {?it was/they were} written: {.file {changed}}.",
      "i" = "The copy may have been cut short, or a file edited. Copy the folder again; to change the run, write a new bundle."
    ), class = "rTFT_bundle_changed")
  invisible(found)
}

# The identity of the bundle in `dir`: the checksum of its list of checksums, which
# changes when any file of the bundle does. Every result of a task is stamped with
# it, so that no result is assembled with a bundle other than the one it came from.
.bundle_id <- function(dir) {
  unname(tools::md5sum(file.path(dir, "bundle", "checksums.txt")))
}


# ===========================================================================
# 2. Writing the bundle
# ===========================================================================

# write_hpc_bundle() -----------------------------------------------------------

#' Write a table of models into a folder to fit on a cluster
#'
#' Writes what is needed to fit the models of a table built with
#' [build_hazard_models()] on a cluster that schedules its jobs with Slurm, such as
#' those of the Digital Research Alliance of Canada, into one folder, the bundle.
#' The bundle is then copied to the cluster, where the scripts it contains install what
#' it needs (`install.sh`) and check it (`check.sh`, which calls
#' [check_hpc_bundle()]).
#'
#' Everything about the run is settled on the user's own machine, before anything
#' is copied: the arguments are checked as [fit_hazard_models()] checks them, the
#' unit left out in cross-validation is chosen and, with `kfold = TRUE`, the units
#' are assigned to folds once for all the models. A mistake therefore is quickly
#' found here and not on the cluster hours later.
#'
#' The bundle holds:
#'
#' * `settings.sh`, the account, time, memory, cores and modules, written once and
#'   read by every script of the bundle. It is the one file to edit when a task
#'   needs more time or memory.
#' * `install.sh` and `install.R`, which install what the bundle needs into the R of
#'   the cluster: rTFT at the version that wrote the bundle, from `rtft_source`,
#'   cmdstanr, and CmdStan. It is run once on each cluster, on the login node or, with
#'   `bash install.sh --job`, as a Slurm job.
#' * `check.sh`, which checks the bundle and the R of the cluster on the login node.
#' * `submit.sh`, which submits the tasks to Slurm as job arrays, one for the fits
#'   and one for the refits, queued to start when the fits have ended; `task.sh`,
#'   which Slurm runs for each task and which calls [run_hpc_task()]; and
#'   `status.sh`, which calls [hpc_status()] to say what each task has done.
#' * `README.txt`, the steps to take, with the commands filled in.
#' * `bundle/`, the data (`data.rds`), the models (`models.rds`), the settings of the
#'   run (`settings.rds`), the tasks as tables to read (`tasks-fit.csv`, and
#'   `tasks-kfold.csv` with `kfold = TRUE`), and `checksums.txt`, the checksum of
#'   each of these files, against which the cluster checks what it received.
#' * `stan/`, `results/` and `logs/`, empty, for the compiled models, the results and
#'   what each task prints.
#'
#' The fits on the cluster use the cmdstanr backend, and each task runs its chains on
#' as many cores as there are chains, so `backend` and `cores` cannot be given in
#' `...`.
#'
#' @param models A data frame, the table of models from [build_hazard_models()].
#' @param long A data frame, the table of intervals at risk from
#'   [build_hazard_data()], to which every model is fitted.
#' @param dir A character string, the folder to write the bundle into. It is created
#'   if it does not exist, and must be empty if it does.
#' @param account A character string, the allocation the jobs are charged to, as
#'   Slurm's `--account` takes it (for example `"def-someone"`).
#' @param observation A character string, how the failure times were observed, as in
#'   [fit_hazard()].
#' @param unit An optional character string, the unit to leave out in
#'   cross-validation, as in [loo_hazard()]. `NULL` (the default) chooses it as
#'   [loo_hazard()] does.
#' @param kfold A single logical value. `TRUE` also cross-validates each model by
#'   refitting, one task per model and fold. Defaults to `FALSE`.
#' @param K An optional single whole number, the number of folds when `kfold` is
#'   `TRUE`, as in [kfold_hazard()]. `NULL` (the default) gives one fold per unit.
#'   With fewer folds than units the units are assigned to folds at random, once;
#'   call [set.seed()] first to get the same folds again.
#' @param rhat,ess,divergent The limits of the convergence rule, as in
#'   [fit_hazard_models()].
#' @param time A character string, the time each task may run, as Slurm's `--time`
#'   takes it: hours, minutes and seconds (`"12:00:00"`, the default), or days
#'   first (`"1-12:00:00"`).
#' @param mem_per_cpu A character string, the memory for each core, with its unit, as
#'   Slurm's `--mem-per-cpu` takes it. Defaults to `"4000M"`. Slurm reads `"4G"` as
#'   4096 M, which the clusters of the Alliance point out at each submission.
#' @param email An optional character string, the address Slurm writes to when a task
#'   fails. `NULL` (the default) sends no mail.
#' @param modules A character vector, the modules that give R on the cluster. The
#'   default, `c("StdEnv/2023", "r/4.5.0")`, is for the clusters of the Alliance;
#'   `module avail r` on the cluster lists the versions it offers.
#' @param r_library An optional character string, the R library on the cluster that
#'   holds rTFT, brms and cmdstanr. `NULL` (the default) uses R's own personal
#'   library.
#' @param remote An optional character string, where on the cluster to copy the
#'   bundle, as `"user@host:folder"` (for example
#'   `"someone@nibi.alliancecan.ca:scratch/my-run"`), without spaces. A folder
#'   written without a leading `/` is in your home folder on the cluster. It is used
#'   only to fill in the commands of `README.txt` and of the message, and nothing is
#'   sent to the cluster; `NULL` (the default) leaves placeholders there. The copy
#'   command creates the folder, with any folders above it that do not exist yet.
#' @param rtft_source An optional character string, where the cluster installs rTFT
#'   from, as `remotes::install_github()` takes it. `NULL` (the default) uses the
#'   commit the installed rTFT came from, when it was installed from GitHub, or the
#'   tag of its release; a development version installed from a copy on your machine
#'   must be given, such as `"eogmartins/rTFT@dev-0.2.0"`.
#' @param template A character string, `"alliance"` (the default) for the scripts
#'   written for the clusters of the Alliance, or the path of a folder of your own
#'   holding `settings.sh`, `install.sh`, `install.R`, `check.sh`, `submit.sh`,
#'   `task.sh`, `status.sh` and `README.txt`, with the same placeholders.
#'   [hpc_status()] finds the log of a task as `logs/<stage>-<array job>_<task>.out`;
#'   a template that names its logs otherwise loses only the reason given for a task
#'   stopped by Slurm.
#' @param ... Further arguments for every fit, as in [fit_hazard_models()], such as
#'   `chains`, `iter`, `warmup`, `seed` or `control`, or the priors that
#'   [fit_hazard()] takes.
#' @return The path of the bundle, invisibly. A message says what was written and
#'   what to do next.
#' @seealso [check_hpc_bundle()], which checks the bundle on the cluster,
#'   [hpc_status()], which says what its tasks have done, and [collect_hpc()], which
#'   assembles their results.
#' @examples
#' # a small ramping assay in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1,
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = "trial",
#'                             checks = checks)
#' models <- build_hazard_models("exp", c("mem", "weibull"),
#'                               dpar_formulas = list(loga = ~ 1 + (1 | trial)))
#'
#' # the bundle, written into a temporary folder here; rtft_source says the cluster
#' # installs rTFT from its repository on GitHub (by default it is read from how rTFT
#' # was installed on this machine)
#' bundle <- file.path(tempdir(), "my-run")
#' write_hpc_bundle(models, long, bundle, account = "def-someone", kfold = TRUE,
#'                  remote = "someone@nibi.alliancecan.ca:scratch/my-run",
#'                  rtft_source = "eogmartins/rTFT", chains = 4, iter = 2000, seed = 1)
#' list.files(bundle, recursive = TRUE)
#' @export
write_hpc_bundle <- function(models, long, dir, account,
                             observation = c("interval", "continuous", "mixed"),
                             unit = NULL, kfold = FALSE, K = NULL,
                             rhat = 1.01, ess = 400, divergent = 0.001,
                             time = "12:00:00", mem_per_cpu = "4000M", email = NULL,
                             modules = c("StdEnv/2023", "r/4.5.0"), r_library = NULL,
                             remote = NULL, rtft_source = NULL, template = "alliance",
                             ...) {
  dots <- list(...)
  other_backend <- !is.null(dots$backend) && !identical(dots$backend, "cmdstanr")
  if ("cores" %in% names(dots) || other_backend)
    cli::cli_abort(c(
      "{.arg cores} and {.arg backend} are set for the cluster, not given here.",
      "i" = "Each task runs its chains on as many cores as there are chains, with the {.val cmdstanr} backend."
    ), class = "rTFT_set_in_bundle")
  if (missing(account)) account <- NULL
  .check_hpc_settings(account, time, mem_per_cpu, email, modules, r_library, remote)
  template_dir <- .template_dir(template)
  installed <- .installed_rtft_version()
  given     <- !is.null(rtft_source)
  if (!given) rtft_source <- .rtft_source()
  if (!is.character(rtft_source) || length(rtft_source) != 1L || is.na(rtft_source) ||
      !nzchar(rtft_source))
    cli::cli_abort(c(
      "The cluster cannot be told where to install this version of rTFT from.",
      "x" = if (given) "You supplied {.val {rtft_source}}."
            else "It is the development version {.val {installed}}, installed from a copy on this machine.",
      "i" = "Give {.arg rtft_source}, as {.fn remotes::install_github} takes it, such as {.val eogmartins/rTFT@dev-0.2.0}. The branch must be on GitHub."
    ), class = "rTFT_no_rtft_source")
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir))
    cli::cli_abort("{.arg dir} must be a single character string, a folder.",
                   class = "rTFT_bad_dir")
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE)))
    cli::cli_abort(c(
      "{.path {dir}} is not empty.",
      "i" = "Give a new folder, or an empty one: a bundle is never written over another."
    ), class = "rTFT_dir_not_empty")
  # the run itself, checked last, since it also chooses the unit and says so
  plan <- .models_plan(models, long, observation, unit, kfold, K, rhat, ess, divergent,
                       dots)
  .check_distinct_models(models)

  # the folds, assigned once for all the models
  folds <- NULL
  if (kfold) {
    K <- .check_K(K, length(plan$units))
    folds <- .assign_folds(plan$units, K)
  }
  dots$backend <- "cmdstanr"
  settings <- list(
    rtft_version = installed, rtft_commit = .installed_rtft_commit(),
    brms_version = as.character(utils::packageVersion("brms")),
    created      = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    observation  = plan$observation, unit = plan$unit, units = plan$units,
    kfold = kfold, K = if (kfold) K else NULL, folds = folds,
    rule = list(rhat = rhat, ess = ess, divergent = divergent),
    chains = plan$chains, dots = dots,
    slurm = list(account = account, time = time, mem_per_cpu = mem_per_cpu,
                 cpus_per_task = plan$chains, email = email, modules = modules,
                 r_library = r_library))

  # the tasks, as tables to read: one per model, and one per model and fold
  tasks_fit <- data.frame(task = seq_len(nrow(models)), model = models$model,
                          m_form = models$m_form, baseline = models$baseline)
  tasks_kfold <- if (kfold) {
    grid <- expand.grid(fold = seq_len(K), model = models$model,
                        stringsAsFactors = FALSE)
    of_fold <- \(k) paste(names(folds)[folds == k], collapse = ";")
    data.frame(task = seq_len(nrow(grid)), model = grid$model, fold = grid$fold,
               units = vapply(grid$fold, of_fold, character(1)))
  }

  # the scripts, filled in from the template
  ssh_target    <- if (is.null(remote)) "<user>@<cluster>" else sub(":.*$", "", remote)
  remote_folder <- if (is.null(remote)) "<folder>" else sub("^[^:]*:", "", remote)
  # the folder does not exist yet, so its full path is read from the folder above it
  local_path    <- file.path(normalizePath(dirname(dir), winslash = "/",
                                           mustWork = FALSE), basename(dir))
  copy <- .copy_command(local_path, ssh_target, remote_folder, .on_windows())
  copy_back <- .copy_back_command(local_path, ssh_target, remote_folder, .on_windows())
  # what the bundle holds, as a sentence with its plurals, for README.txt
  plural   <- function(n, word) paste(n, if (n == 1) word else paste0(word, "s"))
  n_refits <- if (kfold) nrow(tasks_kfold) else 0
  contents <- paste0(
    "It holds ", plural(nrow(models), "model"), " to fit, as ",
    plural(nrow(tasks_fit), "task"), ", ",
    if (kfold)
      paste("and", plural(n_refits, "refit"), "for cross-validation by refitting")
    else "without cross-validation by refitting", ".")
  values <- list(
    created = settings$created, account = account,
    email = if (is.null(email)) "" else email,
    time = time, mem_per_cpu = mem_per_cpu, cpus_per_task = plan$chains,
    modules = paste(modules, collapse = " "),
    r_library = if (is.null(r_library)) "" else r_library,
    rtft_version = settings$rtft_version, rtft_source = rtft_source,
    rtft_commit = if (is.na(settings$rtft_commit)) "" else settings$rtft_commit,
    contents = paste(strwrap(contents, width = 80), collapse = "\n"),
    n_models = nrow(models), n_fit_tasks = nrow(tasks_fit),
    n_kfold_tasks = n_refits,
    copy_command = paste(copy$command, collapse = "\n       "),
    copy_note = paste(strwrap(copy$note, width = 74), collapse = "\n     "),
    copy_back_command = paste(copy_back, collapse = "\n       "),
    local_path = local_path,
    ssh_target = ssh_target, remote_folder = remote_folder)
  # Filled before anything is written, so that a template that cannot be filled
  # leaves no half-written bundle behind. A plain loop, so that its error reaches
  # the user as it was raised and not wrapped in the error of a mapping function.
  scripts <- list()
  for (file in .script_files)
    scripts[[file]] <- .fill_template(readLines(file.path(template_dir, file)), values)

  # the bundle itself
  for (folder in c("bundle", "stan", "results", "logs"))
    dir.create(file.path(dir, folder), recursive = TRUE, showWarnings = FALSE)
  saveRDS(long, file.path(dir, "bundle", "data.rds"))
  saveRDS(models, file.path(dir, "bundle", "models.rds"))
  saveRDS(settings, file.path(dir, "bundle", "settings.rds"))
  .write_csv_lf(tasks_fit, file.path(dir, "bundle", "tasks-fit.csv"))
  if (kfold)
    .write_csv_lf(tasks_kfold, file.path(dir, "bundle", "tasks-kfold.csv"))
  sums <- .bundle_checksums(dir)
  .write_lines_lf(paste(sums, names(sums), sep = "  "),
                  file.path(dir, "bundle", "checksums.txt"))

  for (file in .script_files) .write_lines_lf(scripts[[file]], file.path(dir, file))

  refits <- if (kfold)
    paste(" and", nrow(tasks_kfold), "refits for cross-validation") else ""
  label  <- if (plan$unit == "id") "individual" else plan$unit
  cli::cli_inform(c(
    "v" = "Wrote a bundle of {nrow(models)} model{?s} to {.path {dir}}: {nrow(tasks_fit)} fit{?s}{refits}, leaving out one {label} at a time.",
    "i" = "Copy it to the cluster with the command in step 1 of {.file README.txt}, in that folder.",
    "i" = "Then, on the cluster, in that folder: {.code bash install.sh} (once on each cluster), then {.code bash check.sh}."
  ), class = "rTFT_bundle_written")
  invisible(normalizePath(dir))
}


# ===========================================================================
# 3. Checking the bundle on the cluster
# ===========================================================================

# check_hpc_bundle() -----------------------------------------------------------

#' Check a bundle and the R that will fit it, before anything is submitted
#'
#' Checks a bundle written by [write_hpc_bundle()], where it is to be run: on the
#' login node of the cluster, from the bundle's script `check.sh`, which loads R as
#' the bundle's settings say and calls this function. It takes seconds, fits nothing
#' and installs nothing.
#'
#' It checks, in order, that the files of the bundle are all there and unchanged
#' since they were written (by their checksums, so a copy that was cut short or a
#' file edited by hand is found); that the installed rTFT is the version that wrote
#' the bundle; that cmdstanr is installed and finds CmdStan; and that every model of
#' the bundle can be assembled on its data. The first problem stops the check, with
#' an error that says what to do about it.
#'
#' The version is compared by its number and, when both the rTFT that wrote the
#' bundle and the one installed came from GitHub, by the commit each was installed
#' from, since copies of a development version share a number and can differ in
#' their code. A copy installed from a folder on your own machine has no commit to
#' compare, so for it only the number is checked. [run_hpc_task()] makes the same
#' check at the start of each task, since the library can change while tasks wait
#' in the queue.
#'
#' @param dir A character string, the folder of the bundle. Defaults to the working
#'   folder, which is where `check.sh` runs it.
#' @return `TRUE`, invisibly, when every check passes. A message lists what was
#'   checked.
#' @seealso [write_hpc_bundle()], which writes the bundle.
#' @examples
#' \dontrun{
#' # on the cluster, from the folder of the bundle; check.sh runs the same
#' check_hpc_bundle()
#' }
#' @export
check_hpc_bundle <- function(dir = ".") {
  # every file as it was written
  found <- .verify_bundle(dir)
  settings <- readRDS(file.path(dir, "bundle", "settings.rds"))
  models   <- readRDS(file.path(dir, "bundle", "models.rds"))

  # the software that will fit the models
  installed <- .check_rtft_version(settings)
  if (!.has_cmdstanr())
    cli::cli_abort(c(
      "The {.pkg cmdstanr} package is not installed in this R.",
      "i" = "Run {.code bash install.sh}, which installs cmdstanr and then CmdStan, which cmdstanr needs."
    ), class = "rTFT_needs_cmdstanr")
  if (!.cmdstan_found())
    cli::cli_abort(c(
      "{.pkg cmdstanr} cannot find CmdStan, which compiles the models.",
      "i" = "Run {.code bash install.sh}, which compiles CmdStan; it can take 10 to 20 minutes. If the login node stops it, run {.code bash install.sh --job}."
    ), class = "rTFT_no_cmdstan")

  # every model as it will be fitted
  for (i in seq_len(nrow(models)))
    build_hazard_spec(models$m_form[i], models$baseline[i], models$dpar_formulas[[i]],
                      settings$observation, costressors = models$costressors[[i]])

  brms   <- as.character(utils::packageVersion("brms"))
  refits <- if (isTRUE(settings$kfold))
    paste(", with", nrow(models) * settings$K, "refits for cross-validation") else ""
  cli::cli_inform(c(
    "v" = "Every file of the bundle is as it was written ({length(found)} file{?s}).",
    "v" = "rTFT {installed}, the version that wrote the bundle; brms {brms}; CmdStan found.",
    "v" = "{nrow(models)} model{?s} assembled{refits}.",
    "i" = "Submit the tasks with {.code bash submit.sh}."
  ), class = "rTFT_bundle_checked")
  invisible(TRUE)
}


# ===========================================================================
# 4. Running one task on the cluster
# ===========================================================================

# .check_stage() / .read_tasks() / .check_task() / .task_label() ---------------

# Check that `stage` is "fit" or "kfold", and return it.
.check_stage <- function(stage) {
  if (!is.character(stage) || length(stage) != 1L || is.na(stage) ||
      !stage %in% c("fit", "kfold"))
    cli::cli_abort(c(
      "{.arg stage} must be {.val fit} or {.val kfold}.",
      "x" = "You supplied {.val {stage}}."
    ), class = "rTFT_bad_stage")
  stage
}

# The table of the tasks of `stage` in the bundle in `dir`, as a data frame, or NULL
# for the refits of a bundle written without them. Every column is read as text and
# none as missing, so that a model named "NA", "1" or "TRUE" keeps its name; the
# numbers of the tasks and the folds are then made whole numbers.
.read_tasks <- function(dir, stage) {
  file <- file.path(dir, "bundle", paste0("tasks-", stage, ".csv"))
  if (!file.exists(file)) return(NULL)
  tasks <- utils::read.csv(file, colClasses = "character", na.strings = character(0))
  for (column in intersect(c("task", "fold"), names(tasks)))
    tasks[[column]] <- as.integer(tasks[[column]])
  tasks
}

# The row of `tasks`, the table of the tasks of `stage`, that holds task `task`, given
# as a number or as text holding one, as Slurm gives it. A task the table does not
# have is an error, so an array that asks for more tasks than there are fails loudly.
.check_task <- function(task, tasks, stage) {
  if (is.null(tasks))
    cli::cli_abort(c(
      "This bundle has no refits for cross-validation.",
      "i" = "It was written with {.code kfold = FALSE}; a bundle written with {.code kfold = TRUE} has them."
    ), class = "rTFT_no_task")
  # empty: run by hand, outside a Slurm job, with no task given
  if (identical(task, ""))
    cli::cli_abort(c(
      "No task was given, and {.envvar SLURM_ARRAY_TASK_ID}, which Slurm sets for each task of a job array, is not set.",
      "i" = "Give {.arg task}, such as {.code run_hpc_task(stage = \"{stage}\", task = 1)}."
    ), class = "rTFT_no_task")
  number <- suppressWarnings(as.numeric(task))
  if (length(number) != 1L || is.na(number) || !number %in% tasks$task)
    cli::cli_abort(c(
      "There is no task {.val {task}} among the {nrow(tasks)} task{?s} of the {.val {stage}} stage of this bundle.",
      "i" = "The tasks are numbered from 1 to {nrow(tasks)}. On Slurm the number is that of the array task, so an {.code --array} beyond it ends here."
    ), class = "rTFT_no_task")
  match(number, tasks$task)
}

# How messages name the tasks in the rows of `row`, of the stages `stage` (one, or
# one for each row): "fit 2 (exp_mem)", or "refit 3 (exp_mem, fold 1)".
.task_label <- function(stage, row) {
  ifelse(rep_len(stage, nrow(row)) == "fit",
         sprintf("fit %d (%s)", row$task, row$model),
         sprintf("refit %d (%s, fold %d)", row$task, row$model, row$fold))
}


# .task_file() / .write_whole() / .read_record() / .slurm_job() ----------------
# .one_line() / .duration()

# The path of a file of task `task` of `stage` in the folder results/ of `dir`, with
# the extension `ext`: "rds" for its result, "csv" for the record of how it ended,
# and "started" for the record of when and where it began.
.task_file <- function(dir, stage, task, ext) {
  file.path(dir, "results", sprintf("%s-%03d.%s", stage, as.integer(task), ext))
}

# Write `x` to the file `path` with `write`, a function of an object and a path,
# under another name first and then renamed, so that a task stopped while it writes
# (e.g., at the end of its wall time) leaves no file that looks complete and is not.
# Return `path`, invisibly.
.write_whole <- function(x, path, write) {
  part <- paste0(path, ".part")
  write(x, part)
  # file.rename() warns as well as failing; the error below says it once
  if (!suppressWarnings(file.rename(part, path)))
    cli::cli_abort("{.path {path}} could not be written.", class = "rTFT_not_written")
  invisible(path)
}

# The record of how task `task` of `stage` ended, as a list of one character value
# per column, or NULL when it has none: it has not ended, or not started. Every column
# is read as text, so that a checksum made only of digits stays as it was written.
.read_record <- function(dir, stage, task) {
  file <- .task_file(dir, stage, task, "csv")
  if (!file.exists(file)) return(NULL)
  as.list(utils::read.csv(file, colClasses = "character"))
}

# The job and the array task that run this task, as Slurm names them in the names of
# its logs ("1234_5"), or NA outside Slurm.
.slurm_job <- function() {
  job   <- Sys.getenv("SLURM_ARRAY_JOB_ID")
  index <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  if (!nzchar(job) || !nzchar(index)) return(NA_character_)
  paste0(job, "_", index)
}

# The message of the condition `cnd` as one line of plain ASCII text, to keep in a
# table that may be written where R does not run in UTF-8: the bullets of cli in
# their ASCII form, and the curly quotes of base R's messages made straight.
.one_line <- function(cnd) {
  old <- options(cli.unicode = FALSE)
  on.exit(options(old))
  x <- cli::ansi_strip(conditionMessage(cnd))
  x <- gsub("[\u2018\u2019]", "'", x)
  x <- gsub("[\u201c\u201d]", "\"", x)
  gsub("[[:space:]]+", " ", trimws(x))
}

# A time of `seconds` as text, in seconds, minutes or hours, the unit chosen after
# rounding, so that 119.6 seconds is "2 min" and not "120 s".
.duration <- function(seconds) {
  if (round(seconds) < 120) return(paste(round(seconds), "s"))
  if (round(seconds / 60) < 120) return(paste(round(seconds / 60), "min"))
  paste(round(seconds / 3600, 1), "h")
}


# .fit_task() / .kfold_task() --------------------------------------------------

# Fit the model of `row`, a task of the fits, to the data of the bundle in `dir`
# with the bundle's `settings`, as fit_hazard_models() fits each model, with the
# chains on as many cores as there are chains, and return what .fit_one_model()
# returns: the fit, its cross-validation by importance sampling, its convergence
# diagnostics and the number of units for which importance sampling cannot be
# trusted.
.fit_task <- function(dir, settings, row) {
  models <- readRDS(file.path(dir, "bundle", "models.rds"))
  long   <- readRDS(file.path(dir, "bundle", "data.rds"))
  .fit_one_model(models[models$model == row$model, ], long, settings$observation,
                 settings$unit, settings$chains,
                 c(settings$dots, list(cores = settings$chains)))
}

# Refit the model of `row`, a task of the refits, without the units of its fold, and
# return, as `elpd` in a list, the log predictive density of each unit of the fold,
# named by unit, as kfold_hazard() computes it. The refit starts from the model's
# fit, read from the result of its task in `dir`, which must have ended well for the
# bundle whose identity is `id`; the bundle's `settings` give the unit, the folds
# and the arguments of the refit.
.kfold_task <- function(dir, settings, row, id) {
  fits   <- .read_tasks(dir, "fit")
  j      <- fits$task[fits$model == row$model]
  record <- .read_record(dir, "fit", j)
  file   <- .task_file(dir, "fit", j, "rds")
  if (is.null(record) || !identical(record$status, "done") ||
      !identical(record$bundle, id) || !file.exists(file))
    cli::cli_abort(c(
      "The refits of {.val {row$model}} start from its fit, task {j} of the fits, which has not ended well.",
      "i" = "Run that task first; {.code bash submit.sh --fit {j} --kfold {row$task}} submits the two, the refit to start when the fit has ended."
    ), class = "rTFT_fit_not_ready")
  fit   <- readRDS(file)$fit
  setup <- withCallingHandlers(
    .kfold_setup(fit, settings$unit, folds = settings$folds),
    rTFT_loo_leak = function(w) invokeRestart("muffleWarning"))
  row_fold <- unname(setup$folds[as.character(setup$unit$values)])
  fold <- do.call(.refit_fold, c(
    list(fit, setup$data, row_fold != row$fold, setup$unit$values),
    .refit_args(settings$dots), list(cores = settings$chains)))
  list(elpd = fold$elpd)
}


# run_hpc_task() ---------------------------------------------------------------

#' Run one task of a bundle on the cluster
#'
#' Runs one task of a bundle written by [write_hpc_bundle()]: the fit of one model,
#' or one refit of a model for cross-validation, without the units of one fold. It is
#' called on a compute node of the cluster by the bundle's script `task.sh`, which
#' Slurm runs once for each task of a job array, and is seldom called by hand.
#'
#' A fit is made as [fit_hazard_models()] makes it: with [fit_hazard()], refused if it
#' lost a chain, and cross-validated by importance sampling with [loo_hazard()],
#' leaving out the unit settled when the bundle was written. A refit starts from the
#' model's fit, read from the result of its task, so it runs only once that fit has
#' ended well, and it predicts the units left out as [kfold_hazard()] does. Each
#' model is compiled once, by its fit, into the folder `stan/` of the bundle, where
#' its refits find the compiled program.
#'
#' Before anything is fitted, the files of the bundle are checked against their
#' checksums and the installed rTFT against the version that wrote the bundle; a task
#' stopped by the second check is recorded as failed, with the reason. A task whose
#' result is already in `results/` is not run again, so the tasks of a stage can be
#' submitted again to run only those that failed or did not finish.
#'
#' Each task writes into `results/`: a record of when and where it started
#' (`fit-001.started`, for the first fit), and when it ends, a record of how it ended
#' (`fit-001.csv`) and, if it ended well, its result (`fit-001.rds`). The record of
#' how an earlier attempt ended is moved aside when a task starts again
#' (`fit-001.previous.csv`), so that a task running again is seen as started. A
#' result is written under another name and renamed when complete, so a task stopped
#' at the end of its wall time leaves no result that looks complete. Every result is
#' stamped with the checksums of the bundle, so that [collect_hpc()] never assembles
#' it with another bundle.
#'
#' @param dir A character string, the folder of the bundle. Defaults to the working
#'   folder, which is where `task.sh` runs it.
#' @param stage A character string, `"fit"` (the default) for a fit, or `"kfold"` for a
#'   refit for cross-validation.
#' @param task A single whole number, or a character string holding one, the number of
#'   the task among those of its stage, as the tables `bundle/tasks-fit.csv` and
#'   `bundle/tasks-kfold.csv` number them. Defaults to the number Slurm gives the
#'   array task, in the environment variable `SLURM_ARRAY_TASK_ID`.
#' @return The path of the task's result, invisibly. A message says how the task
#'   ended. A task that fails is an error, raised after its reason has been written
#'   into `results/`.
#' @seealso [write_hpc_bundle()], which writes the bundle, and [hpc_status()], which
#'   says what each task has done.
#' @examples
#' \dontrun{
#' # on the cluster, from the folder of the bundle; task.sh runs the same, with the
#' # number of the task given by Slurm
#' run_hpc_task(stage = "fit")
#'
#' # the second fit of the bundle in the folder my-run, run by hand
#' run_hpc_task("my-run", stage = "fit", task = 2)
#' }
#' @export
run_hpc_task <- function(dir = ".", stage = "fit",
                         task = Sys.getenv("SLURM_ARRAY_TASK_ID")) {
  stage <- .check_stage(stage)
  .verify_bundle(dir)
  settings  <- readRDS(file.path(dir, "bundle", "settings.rds"))
  tasks     <- .read_tasks(dir, stage)
  row       <- tasks[.check_task(task, tasks, stage), ]
  id     <- .bundle_id(dir)
  label  <- .task_label(stage, row)
  result <- .task_file(dir, stage, row$task, "rds")

  # a task that has ended well is not run again
  ended <- .read_record(dir, stage, row$task)
  if (!is.null(ended) && identical(ended$status, "done") &&
      identical(ended$bundle, id) && file.exists(result)) {
    cli::cli_inform(c("v" = "Task {label} has its result already, and is not run again."),
                    class = "rTFT_task_skipped")
    return(invisible(result))
  }

  host <- Sys.info()[["nodename"]]
  job  <- .slurm_job()
  now  <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  dir.create(file.path(dir, "results"), showWarnings = FALSE)
  # The record of how an earlier attempt ended is moved aside, so that while this one
  # runs the task is seen as started, and not as failed and to be submitted again;
  # the reason of the earlier attempt is kept beside it.
  if (!is.null(ended)) {
    previous <- .task_file(dir, stage, row$task, "previous.csv")
    if (!suppressWarnings(file.rename(.task_file(dir, stage, row$task, "csv"),
                                      previous)))
      cli::cli_abort("{.path {previous}} could not be written.",
                     class = "rTFT_not_written")
  }
  .write_whole(data.frame(stage = stage, task = row$task, job = job, host = host,
                          started = now()),
               .task_file(dir, stage, row$task, "started"), .write_csv_lf)
  cli::cli_inform("Running task {label}, on {host}.", class = "rTFT_task_progress")

  # Each model is compiled into the folder stan/ of the bundle, where its refits find
  # it: cmdstanr reuses a compiled program whose code is the same. The chains are
  # given their cores directly, and mc.cores is left alone: set, it would make
  # other functions copy this process with forks, which can hang for good once the
  # threads of a loaded sampler are running (seen on Nibi, 2026-10-19).
  old <- options(cmdstanr_write_stan_file_dir = normalizePath(file.path(dir, "stan")))
  on.exit(options(old), add = TRUE)
  begun <- proc.time()[["elapsed"]]
  # The version is checked here, among what the task records, so that a task stopped
  # by it says why in results/, where hpc_status() reads it.
  res <- tryCatch({
    .check_rtft_version(settings)
    if (stage == "fit") .fit_task(dir, settings, row)
    else .kfold_task(dir, settings, row, id)
  }, error = function(e) e)
  seconds <- proc.time()[["elapsed"]] - begun
  failed  <- inherits(res, "error")

  record <- data.frame(
    stage = stage, task = row$task, model = row$model,
    fold = if (stage == "kfold") row$fold else NA_integer_,
    status = if (failed) "failed" else "done",
    message = if (failed) .one_line(res) else NA_character_,
    seconds = round(seconds), bundle = id,
    rtft_version = .installed_rtft_version(), job = job,
    host = host, ended = now())
  # the result first and the record of its end last, so that a record that says done
  # always has its result beside it
  if (!failed)
    .write_whole(c(res, list(task = row$task, model = row$model, fold = record$fold,
                             seconds = seconds, bundle = id)),
                 result, saveRDS)
  .write_whole(record, .task_file(dir, stage, row$task, "csv"), .write_csv_lf)
  if (failed)
    cli::cli_abort(c(
      "Task {label} failed.",
      "x" = "{record$message}",
      "i" = "The reason is kept in {.file results/}, and {.code bash status.sh} lists every task that failed."
    ), class = "rTFT_task_failed")
  cli::cli_inform(c("v" = "Task {label} ended well, in {(.duration(seconds))}."),
                  class = "rTFT_task_done")
  invisible(result)
}


# ===========================================================================
# 5. Following the run
# ===========================================================================

# .array_ranges() / .expand_ranges() / .read_submissions() / .log_reason() ------
# .task_states()

# The whole numbers `x` as Slurm's --array takes them, each run of consecutive
# numbers shortened to its ends: "1-3,5,7-8".
.array_ranges <- function(x) {
  x   <- sort(unique(as.integer(x)))
  run <- cumsum(c(1L, diff(x) != 1L))
  ends <- vapply(split(x, run), \(r) {
    if (length(r) == 1L) as.character(r) else paste0(r[1], "-", r[length(r)])
  }, character(1))
  paste(ends, collapse = ",")
}

# The whole numbers that `x`, a list as Slurm's --array takes it ("1-3,5"), stands
# for: 1, 2, 3 and 5.
.expand_ranges <- function(x) {
  ends <- strsplit(strsplit(x, ",", fixed = TRUE)[[1]], "-", fixed = TRUE)
  unlist(lapply(ends, \(e) seq(as.integer(e[1]), as.integer(e[length(e)]))))
}

# The submissions submit.sh recorded in the folder results/ of `dir`, as a data frame
# with one row per task submitted: its `stage`, `task` number, the `time` it was
# submitted (in seconds since 1970, by the clock of the login node) and the Slurm
# `job`; without rows when none was recorded.
.read_submissions <- function(dir) {
  none <- data.frame(stage = character(0), task = integer(0), time = numeric(0),
                     job = character(0))
  file <- file.path(dir, "results", "submitted.csv")
  if (!file.exists(file)) return(none)
  sent <- utils::read.csv(file, colClasses = "character")
  rows <- purrr::map(seq_len(nrow(sent)), \(i) data.frame(
    stage = sent$stage[i], task = .expand_ranges(sent$tasks[i]),
    time = as.numeric(sent$time[i]), job = sent$job[i]))
  rbind(none, do.call(rbind, rows))
}

# Why the task of `stage` that started as the Slurm job and array task `job` has no
# record of how it ended, as far as the log Slurm wrote for it in the folder logs/ of
# `dir` can tell: it ran out of time, ran out of memory, was stopped by Slurm for a
# reason Slurm gives (a node that failed, say) or was cancelled; or else it may
# still be running. Slurm writes a line of its own into the log when it stops a
# task.
.log_reason <- function(dir, stage, job) {
  log <- file.path("logs", sprintf("%s-%s.out", stage, job))
  if (is.na(job) || !file.exists(file.path(dir, log)))
    return(paste("no result yet: it may still be running, or it stopped before it",
                 "could write one"))
  text <- readLines(file.path(dir, log), warn = FALSE)
  if (any(grepl("DUE TO TIME LIMIT", text, fixed = TRUE)))
    return(sprintf("it ran out of time (raise TIME in settings.sh); see %s", log))
  if (any(grepl("oom[-_]kill|out of memory", text, ignore.case = TRUE)))
    return(sprintf("it ran out of memory (raise MEM_PER_CPU in settings.sh); see %s",
                   log))
  # "CANCELLED AT <time> DUE TO NODE FAILURE": Slurm's reason, in its own words
  stopped <- regmatches(text, regexpr("CANCELLED AT [^ ]+ DUE TO [A-Z ]+", text))
  if (length(stopped))
    return(sprintf("it was stopped by Slurm, due to %s; see %s",
                   tolower(trimws(sub(".* DUE TO ", "", stopped[1]))), log))
  if (any(grepl("CANCELLED", text, fixed = TRUE)))
    return(sprintf("it was cancelled; see %s", log))
  sprintf("no result yet: it may still be running, or it stopped without one; see %s",
          log)
}

# The state of each task of `stage` in the bundle in `dir`, whose identity is `id`,
# as a data frame with one row per task: its `stage`, `task` number, `model` and
# `fold` (NA for a fit); its `status`, one of "done", "failed", "started" (with no
# record of how it ended), "queued" (submitted since its files were last written,
# and waiting to start) and "not started"; the `reason` a task is not done; and the
# `seconds` it took. NULL for the refits of a bundle written without them.
.task_states <- function(dir, stage, id) {
  tasks <- .read_tasks(dir, stage)
  if (is.null(tasks)) return(NULL)
  sent <- .read_submissions(dir)
  sent <- sent[sent$stage == stage, ]
  one <- function(i) {
    task    <- tasks$task[i]
    record  <- .read_record(dir, stage, task)
    started <- .task_file(dir, stage, task, "started")
    status  <- "not started"
    reason  <- NA_character_
    seconds <- NA_real_
    if (!is.null(record)) {
      seconds <- as.numeric(record$seconds)
      if (!identical(record$bundle, id)) {
        status <- "failed"
        reason <- "its result was written for another bundle"
      } else if (record$status == "done" &&
                 !file.exists(.task_file(dir, stage, task, "rds"))) {
        status <- "failed"
        reason <- "its result is missing from results/"
      } else {
        status <- record$status
        if (status == "failed") reason <- record$message
      }
    } else if (file.exists(started)) {
      status <- "started"
      reason <- .log_reason(dir, stage,
                            utils::read.csv(started, colClasses = "character")$job)
    }
    # A task not done that was submitted after its files were last written waits in
    # the queue: its files still tell of an earlier attempt, if any.
    mine <- sent[sent$task == task, ]
    if (status != "done" && nrow(mine)) {
      files <- .task_file(dir, stage, task, c("csv", "started"))
      written <- suppressWarnings(max(as.numeric(file.info(files)$mtime), na.rm = TRUE))
      last <- mine[which.max(mine$time), ]
      if (last$time > written) {
        status <- "queued"
        reason <- sprintf("submitted at %s as job %s, waiting to start",
                          format(as.POSIXct(last$time), "%Y-%m-%d %H:%M"), last$job)
      }
    }
    tibble::tibble(stage = stage, task = task, model = tasks$model[i],
                   fold = if (stage == "kfold") tasks$fold[i] else NA_integer_,
                   status = status, reason = reason, seconds = seconds)
  }
  purrr::list_rbind(purrr::map(seq_len(nrow(tasks)), one))
}


# hpc_status() -----------------------------------------------------------------

#' What each task of a bundle has done, and what to do next
#'
#' Reads the records that the tasks of a bundle written by [write_hpc_bundle()] have
#' left in its folder `results/`, and says, for each stage of the run, how many tasks
#' are done, which are not and why, and what to do next. It is run on the cluster by
#' the bundle's script `status.sh`, and works as well on a copy of the bundle brought
#' back.
#'
#' A task is `"done"` when it ended well, `"failed"` when it ended with an error, whose
#' message is its reason, `"started"` when it began and has no record of how it ended,
#' `"queued"` when `submit.sh` submitted it after its files were last written, so that
#' it waits in the queue, and `"not started"` otherwise. A task that started and has no
#' record is either still running or was stopped by Slurm before it could write one.
#' Slurm says which in the log it wrote for the task, in the folder `logs/`, and the
#' reason given is read from it: out of time, out of memory, stopped by Slurm for a
#' reason it gives (a node that failed, say), or cancelled. Slurm's command
#' `squeue -u $USER` (`sq` on the clusters of the Alliance) lists the tasks still
#' waiting or running.
#'
#' When tasks remain to be run, the message gives the command that submits them
#' again, and only them: it leaves out a task that may still be running or that
#' waits in the queue. A task already done is never run twice.
#'
#' @param dir A character string, the folder of the bundle. Defaults to the working
#'   folder, which is where `status.sh` runs it.
#' @return A data frame with one row per task, invisibly: its `stage` (`"fit"` or
#'   `"kfold"`), `task` number, `model`, `fold` (`NA` for a fit), `status`, the
#'   `reason` a task is not done, and the `seconds` it took. A message summarises it.
#' @seealso [run_hpc_task()], which runs a task, and [collect_hpc()], which assembles
#'   the results.
#' @examples
#' \dontrun{
#' # on the cluster, from the folder of the bundle; status.sh runs the same
#' hpc_status()
#' }
#' @export
hpc_status <- function(dir = ".") {
  .verify_bundle(dir)
  id <- .bundle_id(dir)
  states <- purrr::list_rbind(list(.task_states(dir, "fit", id),
                                   .task_states(dir, "kfold", id)))

  # Text from the records, the names of the models and the reasons, is shown as it
  # is: braces are doubled so that cli does not read them as its own.
  as_is <- function(x) gsub("([{}])", "\\1\\1", x)
  stage_name <- c(fit = "Fits", kfold = "Refits for cross-validation")
  counts <- character(0)
  for (s in unique(states$stage)) {
    of <- states$status[states$stage == s]
    n  <- table(factor(of, levels = c("done", "failed", "started", "queued",
                                      "not started")))
    others <- n[-1][n[-1] > 0]
    counts <- c(counts, paste0(
      stage_name[[s]], ": ", n[["done"]], " of ", length(of), " done",
      if (length(others)) paste0("; ", paste(others, names(others), collapse = ", ")),
      "."))
  }
  not_done <- states[states$status != "done", ]
  shown    <- utils::head(not_done, 10)
  # Built only when there is a task to list: paste0() recycles empty columns
  # against its constant pieces, and would give one empty line.
  listed <- character(0)
  if (nrow(shown)) {
    reason <- ifelse(is.na(shown$reason), "", paste0(", ", shown$reason))
    listed <- stats::setNames(
      as_is(paste0(.task_label(shown$stage, shown), ": ", shown$status, reason)),
      # a cross for a task that failed, a warning for one that started without
      # ending, and a plain bullet for one that waits or has not started
      c(failed = "x", started = "!", queued = "*", `not started` = "*")[shown$status])
  }
  if (nrow(not_done) > nrow(shown))
    listed <- c(listed, " " = paste("and", nrow(not_done) - nrow(shown), "more."))

  # what to run again: the tasks that failed or never started, and those whose log
  # says Slurm stopped them; never one that may be running or waits in the queue
  stopped <- not_done$status == "started" &
    grepl("^it (ran out|was stopped|was cancelled)", not_done$reason)
  again   <- not_done[not_done$status %in% c("failed", "not started") | stopped, ]
  queued  <- not_done[not_done$status == "queued", ]
  submit  <- function(of) paste(c(
    "bash submit.sh",
    if (any(of$stage == "fit"))
      paste("--fit", .array_ranges(of$task[of$stage == "fit"])),
    if (any(of$stage == "kfold"))
      paste("--kfold", .array_ranges(of$task[of$stage == "kfold"]))),
    collapse = " ")
  command  <- submit(again)
  if_gone  <- submit(queued)
  running  <- any(not_done$status == "started" & !stopped)

  # Not wrapped to the width of the console, so that each command to paste is on one
  # line: wrapped, its second line would be pasted as a command of its own. The
  # commands are plain text and not code, which cli puts between backticks where the
  # console has no colours, and bash would run what lies between them.
  old <- options(cli.width = Inf)
  on.exit(options(old), add = TRUE)
  cli::cli_inform(c(
    stats::setNames(as_is(counts), rep("*", length(counts))),
    listed,
    "v" = if (!nrow(not_done)) "Every task is done. Copy {.file results/} back to your own computer as step 8 of {.file README.txt} says, and assemble the results there with {.fn collect_hpc}.",
    "i" = if (nrow(again)) "To submit again the tasks not done, and only those:",
    " " = if (nrow(again)) "{command}",
    "i" = if (running) "A task started without a record of its end may still be running; {.code sq} lists the tasks waiting or running.",
    "i" = if (nrow(queued)) "A task submitted and waiting to start is not submitted again. One that {.code sq} no longer lists was cancelled before it started; to submit those:",
    " " = if (nrow(queued)) "{if_gone}"
  ), class = "rTFT_hpc_status")
  invisible(states)
}


# ===========================================================================
# 6. Collecting the results
# ===========================================================================

# .read_result() / .collect_kfold() / .say_other_version() ---------------------

# The result of task `task` of `stage` in the bundle in `dir`, whose identity is
# `id`: what run_hpc_task() wrote when the task ended well; an error holding its
# reason, with the `seconds` it took, when it failed; or NULL when it has no record
# of how it ended. A result stamped with another bundle is an error, since it was
# computed from other data or other models, and so is a record or a result that
# cannot be read whole, as a copy cut short leaves it.
.read_result <- function(dir, stage, task, id) {
  csv  <- .task_file(dir, stage, task, "csv")
  file <- .task_file(dir, stage, task, "rds")
  unreadable <- function(what) cli::cli_abort(c(
    "{.file {what}} cannot be read whole.",
    "i" = "The copy back may have been cut short: copy {.file results/} back again."
  ), class = "rTFT_result_unreadable")
  record <- tryCatch(.read_record(dir, stage, task), error = function(e) list())
  if (!file.exists(csv)) return(NULL)
  # one value of each, present: a record cut short can have fewer columns, or none
  # of its one row
  one <- function(x) length(x) == 1L && !is.na(x) && nzchar(x)
  if (!one(record$status) || !one(record$bundle)) unreadable(basename(csv))
  foreign <- function(what) cli::cli_abort(c(
    "{.file {what}} was written for another bundle than the one in {.path {dir}}.",
    "i" = "Results are assembled only with the bundle that computed them: copy the {.file results/} of this bundle into its own folder."
  ), class = "rTFT_result_not_of_bundle")
  if (!identical(record$bundle, id)) foreign(basename(csv))
  if (record$status == "failed") {
    err <- simpleError(record$message)
    err$seconds <- as.numeric(record$seconds)
    return(err)
  }
  if (!file.exists(file)) return(NULL)
  res <- tryCatch(readRDS(file), error = function(e) unreadable(basename(file)))
  if (!identical(res$bundle, id)) foreign(basename(file))
  res
}

# Assemble the results by refitting of each model named in `fitted`, from the refits
# of the bundle in `dir` (whose identity is `id`, and whose `settings` give the unit
# and the folds), as kfold_hazard() assembles them, and return a list of two
# elements: `kfold`, the results by model, and `errors`, the reason each model left
# out has none: a fold with no result, a fold that failed, or a log predictive
# density that is undefined (NaN).
.collect_kfold <- function(dir, settings, fitted, id) {
  tasks  <- .read_tasks(dir, "kfold")
  folds  <- settings$folds[settings$units]
  kfold  <- list()
  errors <- character(0)
  for (m in fitted) {
    mine   <- tasks[tasks$model == m, ]
    elpd   <- stats::setNames(rep(NA_real_, length(folds)), names(folds))
    scored <- stats::setNames(rep(FALSE, length(folds)), names(folds))
    problems <- character(0)
    for (j in seq_len(nrow(mine))) {
      k   <- mine$fold[j]
      res <- .read_result(dir, "kfold", mine$task[j], id)
      if (is.null(res)) {
        problems <- c(problems, sprintf("fold %d has no result", k))
      } else if (inherits(res, "error")) {
        problems <- c(problems,
                      sprintf("fold %d failed (%s)", k, conditionMessage(res)))
      } else if (!setequal(names(res$elpd), names(folds)[folds == k])) {
        problems <- c(problems, sprintf("fold %d scored other units than its own", k))
      } else {
        elpd[names(res$elpd)]   <- res$elpd
        scored[names(res$elpd)] <- TRUE
      }
    }
    # An undefined score (NaN) is reported, never summed into an estimate. A score
    # of -Inf is kept, as kfold_hazard() keeps it: the refit gave the unit no
    # chance at all, which is a result and not an undefined one.
    bad <- names(elpd)[scored & is.na(elpd)]
    if (length(bad))
      problems <- c(problems, paste("the log predictive density is undefined for",
                                    paste(bad, collapse = ", ")))
    if (length(problems)) {
      errors[[m]] <- paste0("Not cross-validated by refitting: ",
                            paste(problems, collapse = "; "), ".")
    } else {
      kfold[[m]] <- .kfold_result(elpd, settings$unit, settings$K, folds)
    }
  }
  list(kfold = kfold, errors = errors)
}


# Say, when the rTFT installed here is not the one that wrote the bundle whose
# `settings` are given (another version, or another commit where both are known),
# which made the fits and which reads them, and return nothing. The fits carry the
# settings that the functions deriving and plotting results read, so a difference
# is worth knowing as the fits are collected, rather than when one of those
# functions fails far from its cause.
.say_other_version <- function(settings) {
  installed <- .installed_rtft_version()
  have  <- .installed_rtft_commit()
  wrote <- if (is.null(settings$rtft_commit)) NA_character_ else settings$rtft_commit
  other_commit <- !is.na(have) && !is.na(wrote) && !identical(have, wrote)
  if (identical(installed, settings$rtft_version) && !other_commit)
    return(invisible(NULL))
  of <- function(version, commit)
    paste0(version, if (!is.na(commit)) paste0(" (commit ", substr(commit, 1, 7), ")"))
  made  <- of(settings$rtft_version, wrote)
  reads <- of(installed, have)
  cli::cli_inform(c(
    "i" = "These fits were made by rTFT {made}; this is rTFT {reads}.",
    " " = "A newer rTFT reads fits of an older one; should a function fail on them, install the version that wrote the bundle to work with them."
  ), class = "rTFT_collected_by_other_version")
  invisible(NULL)
}


# collect_hpc() ----------------------------------------------------------------

#' Assemble the results of a bundle run on a cluster
#'
#' Assembles the results of the tasks of a bundle written by [write_hpc_bundle()],
#' once its folder `results/` has been copied back from the cluster, into what
#' [fit_hazard_models()] returns: the fits, their cross-validation results and a
#' summary of one row per model, ready for [hazard_weights()] and the functions that
#' average models.
#'
#' Nothing is assembled from a bundle whose files have changed since it was written,
#' or with a result computed from another bundle: the files are checked against their
#' checksums, and each result against the stamp of the bundle that computed it.
#'
#' The convergence rule of the bundle is applied to each fit and reported as
#' [fit_hazard_models()] reports it. A model whose fit failed is reported with its
#' reason, and one whose fit has no result yet (its task did not run, or has not
#' finished) is marked `"missing"` in the summary and named in a warning, so the
#' results can be looked at before every task is done. A model is cross-validated by
#' refitting only when every one of its folds has a result with a log predictive
#' density defined for each of its units; otherwise the reason is given in the
#' summary. A density of zero, a log predictive density of `-Inf`, is kept, as
#' [kfold_hazard()] keeps it.
#'
#' @param dir A character string, the folder of the bundle, on your own computer,
#'   with the folder `results/` copied back into it. Defaults to the working folder.
#' @return A list with the four elements that [fit_hazard_models()] returns: `fits`,
#'   `loo`, `kfold` and `summary`. In `summary`, `status` is `"fitted"`, `"failed"` or
#'   `"missing"`, and `seconds` is the time each fit and its cross-validation took on
#'   the cluster.
#' @seealso [write_hpc_bundle()], which writes the bundle, and [hpc_status()], which
#'   says what each task has done.
#' @examples
#' \dontrun{
#' # on your own computer, once results/ has been copied back into the bundle
#' out <- collect_hpc("my-run")
#' out$summary
#' weights <- hazard_weights(out$kfold)
#' }
#' @export
collect_hpc <- function(dir = ".") {
  .verify_bundle(dir)
  id       <- .bundle_id(dir)
  settings <- readRDS(file.path(dir, "bundle", "settings.rds"))
  models   <- readRDS(file.path(dir, "bundle", "models.rds"))
  tasks    <- .read_tasks(dir, "fit")
  rule     <- settings$rule
  kfold    <- isTRUE(settings$kfold)
  .say_other_version(settings)

  n <- nrow(models)
  results <- vector("list", n)
  rows    <- vector("list", n)
  for (i in seq_len(n)) {
    res <- .read_result(dir, "fit", tasks$task[tasks$model == models$model[i]], id)
    no_result <- is.null(res)
    if (no_result)
      res <- simpleError("No result: its task did not run, or has not finished.")
    seconds <- if (is.null(res$seconds)) NA_real_ else res$seconds
    results[[i]] <- res
    rows[[i]] <- .model_row(models[i, ], res, seconds, rule$rhat, rule$ess,
                            rule$divergent)
    if (no_result) rows[[i]]$status <- "missing"
  }
  model_summary <- purrr::list_rbind(rows)
  fitted <- model_summary$status == "fitted"
  named  <- function(what) stats::setNames(purrr::map(results[fitted], what),
                                           models$model[fitted])
  out <- list(fits = named("fit"), loo = named("loo"), kfold = NULL,
              summary = model_summary)

  if (kfold) {
    cv <- .collect_kfold(dir, settings, names(out$fits), id)
    out$kfold   <- cv$kfold
    out$summary <- .kfold_columns(out$summary, cv$kfold, cv$errors)
  }
  absent <- model_summary$model[model_summary$status == "missing"]
  if (length(absent))
    cli::cli_warn(c(
      "!" = "{length(absent)} of the {n} model{?s} {cli::qty(length(absent))}{?has/have} no result yet: {.val {absent}}.",
      "i" = "{.code bash status.sh}, on the cluster, says why; copy {.file results/} back again when {cli::qty(length(absent))}{?it is/they are} done."
    ), class = "rTFT_tasks_missing")
  .report_models(out$summary, rule$rhat, rule$ess, rule$divergent, kfold)
  refitted <- if (kfold)
    paste0(", ", length(out$kfold), " cross-validated by refitting") else ""
  cli::cli_inform(c(
    "v" = "Collected {sum(fitted)} of the {n} model{?s}{refitted}."
  ), class = "rTFT_hpc_collected")
  out
}
