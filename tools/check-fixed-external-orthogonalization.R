# Run against an installed package, not source-overlaid functions. From the
# repository root: Rscript tools/check-fixed-external-orthogonalization.R [R_LIBRARY]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply at most one installed R library directory")
if (length(args)) .libPaths(c(normalizePath(args[1L], mustWork = TRUE), .libPaths()))
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
suppressPackageStartupMessages(library(stochtree))
suppressPackageStartupMessages(library(testthat))
description <- utils::packageDescription("stochtree")
if (!identical(description[["Config/Orthogonalization"]],
               "fixed-external-propensity-20260921")) {
  stop("Install the fixed-external-propensity package before running these checks")
}
cat("Installed package:", find.package("stochtree"), "\n")
cat("Version:", as.character(utils::packageVersion("stochtree")), "\n")
Sys.setenv(NOT_CRAN = "true")
results <- testthat::test_dir(
  file.path(root, "test", "R", "testthat"),
  filter = "^(predict|linear-orthogonalization)$",
  env = new.env(parent = asNamespace("stochtree")),
  reporter = "summary", stop_on_failure = TRUE)
summary <- as.data.frame(results)
cat("Expectations passed:", sum(summary$passed), "\n")
cat("Failures:", sum(summary$failed), "; errors:", sum(summary$error), "\n")
cat("These are integration/regression checks, not a convergence or coverage study.\n")
