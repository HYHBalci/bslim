# Rscript test/R/run-coexistence-test.R <spbart library> [<stochtree library>]
# Checks the two installed packages and their compiled code in one R session.
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args, .libPaths()))
invisible(loadNamespace("stochtree"))
original_predict <- getS3method("predict", "bartmodel")
invisible(loadNamespace("spbart"))
stopifnot(identical(as.character(getNamespaceName(asNamespace("spbart"))), "spbart"))
stopifnot(identical(getS3method("predict", "bartmodel"), original_predict))
stopifnot(identical(environment(getS3method("predict", "spbart_bartmodel")), asNamespace("spbart")))
dlls <- getLoadedDLLs()
stopifnot(all(c("stochtree", "spbart") %in% names(dlls)))
stopifnot(!identical(dlls[["stochtree"]][["path"]], dlls[["spbart"]][["path"]]))
X <- matrix(as.double(seq_len(24)), 8, 3)
for (pkg in c("stochtree", "spbart")) {
  ns <- asNamespace(pkg)
  dataset <- get("createForestDataset", ns)(X)
  rows <- get("dataset_num_rows_cpp", ns)(dataset$data_ptr)
  cols <- get("dataset_num_covariates_cpp", ns)(dataset$data_ptr)
  stopifnot(rows == nrow(X), cols == ncol(X))
  routine <- get(paste0("_", pkg, "_create_forest_dataset_cpp"), ns)
  stopifnot(identical(routine$dll[["name"]], pkg))
}
stopifnot(is.function(spbart::bcf_linear_probit), is.function(stochtree::bcf))
set.seed(835)
X <- matrix(rnorm(120), 40, 3)
y <- X[, 1] + rnorm(40)
for (pkg in c("stochtree", "spbart")) {
  fit <- getExportedValue(pkg, "bart")(
    X_train = X, y_train = y, num_gfr = 2, num_burnin = 2, num_mcmc = 3,
    general_params = list(random_seed = 835, verbose = FALSE),
    mean_forest_params = list(num_trees = 5))
  expected_class <- if (pkg == "spbart") "spbart_bartmodel" else "bartmodel"
  prediction <- predict(fit, X)$y_hat
  stopifnot(inherits(fit, expected_class), nrow(prediction) == nrow(X),
            ncol(prediction) > 0, all(is.finite(prediction)))
}
cat("PASS: stochtree and spbart namespaces, native routines and predictions coexist.\n")
