# check-expectation-counts.R
#
# This reads a manifest of how many expectations each test file is
# meant to execute and refuses any departure in either direction: a drop means a
# check stopped checking; a rise is usually welcome but still has to be written
# down. It also refuses a test file that is on disk but absent from the
# manifest, since such a file runs green and invisible to the count.
#
# Run from the repository root, with skips disabled so the count is comparable:
#   NOT_CRAN=true Rscript tools/check-expectation-counts.R

suppressPackageStartupMessages(library(testthat))
options(rlang_backtrace_on_error = "none")

if (!requireNamespace("pkgload", quietly = TRUE)) {
  cli::cli_abort(c(
    "{.pkg pkgload} is needed to load the package from source.",
    "i" = "Install it rather than skipping the check."
  ))
}

root <- getwd()
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  cli::cli_abort(c(
    "No {.file DESCRIPTION} found in {.path {root}}.",
    "i" = "Run this from the repository root."
  ))
}

pkgload::load_all(root, quiet = TRUE)

manifest <- utils::read.csv(file.path(root, "tests", "expectation-manifest.csv"),
                            stringsAsFactors = FALSE)

test_dir <- file.path(root, "tests", "testthat")
on_disk  <- sort(basename(list.files(test_dir, pattern = "^test-.*\\.R$")))
unlisted <- setdiff(on_disk, manifest$file)
absent   <- setdiff(manifest$file, on_disk)
if (length(unlisted)) {
  cli::cli_abort(c(
    "Test file{?s} on disk but not in the manifest: {.file {unlisted}}.",
    "i" = "Add {?it/them} to {.file tests/expectation-manifest.csv} with {?its/their} count{?s}."
  ))
}
if (length(absent)) {
  cli::cli_abort("The manifest lists file{?s} that {?is/are} not on disk: {.file {absent}}.")
}

observed <- vapply(manifest$file, function(f) {
  r <- as.data.frame(testthat::test_file(file.path(test_dir, f), reporter = "silent"))
  c(passed = sum(r$passed), failed = sum(r$failed), skipped = sum(r$skipped))
}, numeric(3)) |> t() |> as.data.frame()
observed$file <- manifest$file
res <- merge(manifest, observed, by = "file")
res$delta <- res$passed - res$expectations

cat("\n")
print(res[, c("file", "expectations", "passed", "delta", "failed", "skipped")],
      row.names = FALSE)
cat("\n")

dropped <- res[res$delta < 0, , drop = FALSE]
grew    <- res[res$delta > 0, , drop = FALSE]

if (nrow(dropped)) {
  detail <- stats::setNames(
    sprintf("{.file %s}: %d expected, %d run", dropped$file,
            dropped$expectations, dropped$passed),
    rep("*", nrow(dropped)))
  cli::cli_abort(c(
    "Expectation count fell: a check has stopped checking.",
    detail,
    "i" = "Find out which before touching {.file tests/expectation-manifest.csv}."
  ))
}
if (nrow(grew)) {
  detail <- stats::setNames(
    sprintf("{.file %s}: %d expected, %d run", grew$file,
            grew$expectations, grew$passed),
    rep("*", nrow(grew)))
  cli::cli_abort(c(
    "Expectation count rose.",
    detail,
    "i" = "Update {.file tests/expectation-manifest.csv} deliberately."
  ))
}
if (sum(res$failed) > 0) {
  cli::cli_abort("There {?is/are} {sum(res$failed)} failing expectation{?s}.")
}
if (sum(res$skipped) > 0) {
  cli::cli_abort(c(
    "{sum(res$skipped)} expectation{?s} {?was/were} skipped.",
    "i" = "Run with {.envvar NOT_CRAN=true} so the count is comparable."
  ))
}

cli::cli_alert_success(
  "{sum(res$passed)} expectation{?s} across {nrow(res)} file{?s}: none failed, none skipped."
)
