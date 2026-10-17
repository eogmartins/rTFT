# install.R - install rTFT at the version that wrote the bundle, cmdstanr and CmdStan
#
# Run by install.sh, after it has loaded R as settings.sh says. It is plain R and not
# a function of rTFT, because rTFT may not be installed yet when it runs. Each step is
# skipped when it is already done.

# The library to install into: the one settings.sh names, or else R's own personal
# library, created if it does not exist yet. R adds a personal library to its search
# path only if it exists when R starts, so it is added here by hand the first time.
library_dir <- strsplit(Sys.getenv("R_LIBS"), ":", fixed = TRUE)[[1]][1]
if (is.na(library_dir) || !nzchar(library_dir))
  library_dir <- Sys.getenv("R_LIBS_USER")
library_dir <- path.expand(library_dir)
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))
cat("Installing into", library_dir, "\n\n")

repos <- c(stan = "https://stan-dev.r-universe.dev",
           CRAN = "https://cloud.r-project.org")
not_installed <- function(package) !requireNamespace(package, quietly = TRUE)

# rTFT, at the version that wrote the bundle, with brms and the other packages it
# needs.
wanted <- "{{rtft_version}}"
have   <- if (not_installed("rTFT")) "none" else
  as.character(utils::packageVersion("rTFT"))
if (!identical(have, wanted)) {
  cat("rTFT: installing version", wanted, "from {{rtft_source}}",
      "(installed:", have, ")\n")
  if (not_installed("remotes"))
    install.packages("remotes", lib = library_dir, repos = repos)
  remotes::install_github("{{rtft_source}}", lib = library_dir, repos = repos,
                          upgrade = "never", force = TRUE)
  # The source may hold another version than the one that wrote the bundle (e.g., 
  # a branch that has moved on); then the bundle cannot be run with it.
  got <- as.character(utils::packageVersion("rTFT", lib.loc = library_dir))
  if (!identical(got, wanted))
    stop("{{rtft_source}} holds rTFT ", got, ", and the bundle was written by ",
         wanted, ". Write the bundle again, with rtft_source naming the commit of ",
         wanted, ".", call. = FALSE)
} else {
  cat("rTFT", wanted, "is installed\n")
}

# cmdstanr, the interface between R and CmdStan.
if (not_installed("cmdstanr")) {
  cat("cmdstanr: installing\n")
  install.packages("cmdstanr", lib = library_dir, repos = repos)
} else {
  cat("cmdstanr is installed\n")
}

# CmdStan, which cmdstanr downloads and compiles, into a folder of the home directory.
found <- tryCatch(suppressWarnings(suppressMessages(cmdstanr::cmdstan_path())),
                  error = function(e) NULL)
if (is.null(found) || !dir.exists(found)) {
  # the compiler and make that CmdStan needs, checked first so that a missing one is
  # named rather than found partway through the compile
  cmdstanr::check_cmdstan_toolchain()
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))
  cat("CmdStan: downloading and compiling it on", cores, "cores\n")
  cmdstanr::install_cmdstan(cores = cores)
} else {
  cat("CmdStan is installed at", found, "\n")
}
