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
      "i" = "{.arg account} is the name of your allocation (such as {.val def-someone}), {.arg time} is written as {.val 12:00:00} or {.val 1-00:00:00}, and {.arg mem_per_cpu} as {.val 4G} or {.val 4000M}.",
      "i" = "{.arg remote} is the place on the cluster to copy the bundle to, as {.val user@nibi.alliancecan.ca:scratch/my-run}."
    ), class = "rTFT_bad_hpc_setting")
  invisible(TRUE)
}

# The files of a template, which write_hpc_bundle() fills in and writes into the
# bundle.
.script_files <- c("settings.sh", "install.sh", "install.R", "check.sh", "README.txt")

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


# .copy_command() ---------------------------------------------------------------

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


# .bundle_files() / .bundle_checksums() / .verify_bundle() ---------------------

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
  missing <- needed[!file.exists(file.path(dir, needed))]
  if (length(missing))
    cli::cli_abort(c(
      "{.path {dir}} is not a complete bundle: {.file {missing}} {?is/are} missing.",
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
#'   Slurm's `--mem-per-cpu` takes it. Defaults to `"4G"`.
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
#'   holding `settings.sh`, `install.sh`, `install.R`, `check.sh` and `README.txt`
#'   with the same placeholders.
#' @param ... Further arguments for every fit, as in [fit_hazard_models()], such as
#'   `chains`, `iter`, `warmup`, `seed` or `control`, or the priors that
#'   [fit_hazard()] takes.
#' @return The path of the bundle, invisibly. A message says what was written and
#'   what to do next.
#' @seealso [check_hpc_bundle()], which checks the bundle on the cluster.
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
                             time = "12:00:00", mem_per_cpu = "4G", email = NULL,
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
#' compare, so for it only the number is checked.
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
    "v" = "{nrow(models)} model{?s} assembled{refits}."
  ), class = "rTFT_bundle_checked")
  invisible(TRUE)
}
