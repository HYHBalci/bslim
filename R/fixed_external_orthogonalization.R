# Fixed external-propensity treatment coding for the linear BCF model.
# These helpers deliberately do not fit, sample, or clip propensity scores.
.validate_external_propensity <- function(propensity, n, name = "propensity") {
  shape_ok <- is.null(dim(propensity)) ||
    (is.matrix(propensity) && ncol(propensity) == 1L)
  if (!is.numeric(propensity) || !shape_ok || length(propensity) != n ||
      any(!is.finite(propensity)) || any(propensity < 0 | propensity > 1)) {
    stop(paste0(name, " must be a numeric vector or one-column matrix with ",
                n, " finite externally supplied probabilities in [0, 1]."),
         call. = FALSE)
  }
  as.numeric(propensity)
}

.validate_binary_treatment <- function(Z, n, name = "Z") {
  shape_ok <- is.null(dim(Z)) || (is.matrix(Z) && ncol(Z) == 1L)
  if (!is.numeric(Z) || !shape_ok || length(Z) != n ||
      any(!is.finite(Z)) || any(!Z %in% c(0, 1))) {
    stop(paste0(name, " must be a numeric vector or one-column matrix with ",
                n, " finite binary 0/1 treatments; do not pre-centre treatment ",
                "when orthogonalize = TRUE."), call. = FALSE)
  }
  as.numeric(Z)
}

# Each row of coefficients is one retained forest-aligned draw, in the order
# alpha, main-effect coefficients, interaction coefficients. The CATE is not
# multiplied by Z - propensity: only the outcome prediction uses that coding.
.fixed_external_linear_draws <- function(X, coefficients, pairs, outcome_scale) {
  X <- as.matrix(X)
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("Prediction covariates must be finite numeric values.", call. = FALSE)
  }
  if (!is.matrix(coefficients) || !is.numeric(coefficients) ||
      nrow(coefficients) < 1L || any(!is.finite(coefficients))) {
    stop("The model must contain finite retained orthogonal_coefficients draws.",
         call. = FALSE)
  }
  if (is.null(pairs)) pairs <- matrix(integer(), nrow = 2L, ncol = 0L)
  if (!is.matrix(pairs) || nrow(pairs) != 2L ||
      (length(pairs) && (!is.numeric(pairs) || any(!is.finite(pairs)) ||
                         any(pairs != as.integer(pairs)) ||
                         any(pairs < 1L | pairs > ncol(X))))) {
    stop("The model contains invalid interaction_pairs.", call. = FALSE)
  }
  if (!is.numeric(outcome_scale) || length(outcome_scale) != 1L ||
      !is.finite(outcome_scale) || outcome_scale <= 0) {
    stop("The model contains an invalid outcome_scale.", call. = FALSE)
  }
  if (ncol(coefficients) != 1L + ncol(X) + ncol(pairs)) {
    stop("Retained coefficient dimensions do not match the prediction covariates and interactions.",
         call. = FALSE)
  }
  design <- cbind(1, X)
  if (ncol(pairs)) {
    interactions <- matrix(0, nrow(X), ncol(pairs))
    for (k in seq_len(ncol(pairs))) {
      interactions[, k] <- X[, pairs[1L, k]] * X[, pairs[2L, k]]
    }
    design <- cbind(design, interactions)
  }
  (design %*% t(coefficients)) * outcome_scale
}

.predict_linear_bcf_fixed_external <- function(object, X, Z, propensity = NULL,
                                                rfx_group_ids = NULL,
                                                rfx_basis = NULL, ...) {
  params <- object$model_params
  if (!identical(params$orthogonalization, "fixed-external") ||
      !isTRUE(params$orthogonalize) ||
      !identical(params$propensity_covariate, "none") ||
      !identical(params$propensity_seperate, "none") ||
      isTRUE(params$adaptive_coding) || isTRUE(params$internal_propensity_model)) {
    stop("The model does not contain compatible fixed-external orthogonalisation metadata.",
         call. = FALSE)
  }
  if ((!is.data.frame(X)) && (!is.matrix(X))) {
    stop("X must be a matrix or dataframe.", call. = FALSE)
  }
  if (is.null(object$train_set_metadata)) {
    stop("The model is missing train_set_metadata.", call. = FALSE)
  }
  Z <- .validate_binary_treatment(Z, nrow(X), "Z")
  propensity <- .validate_external_propensity(propensity, nrow(X), "propensity")
  X_processed <- preprocessPredictionData(X, object$train_set_metadata)
  if (!identical(as.integer(ncol(X_processed)), as.integer(params$num_covariates))) {
    stop("Prediction covariates must match those used to train the model.", call. = FALSE)
  }
  treatment_coding <- Z - propensity
  # Propensity is exclusively in the treatment basis, never an extra split
  # variable in the prognostic forest or an extra linear CATE covariate.
  forest_dataset <- createForestDataset(X_processed, as.matrix(treatment_coding))
  mu_hat <- object$forests_mu$predict(forest_dataset) * params$outcome_scale +
    params$outcome_mean
  tau_hat <- .fixed_external_linear_draws(
    X_processed, object$orthogonal_coefficients, object$interaction_pairs,
    params$outcome_scale)
  if (!identical(dim(mu_hat), dim(tau_hat))) {
    stop("Prognostic forest draws and retained orthogonal coefficients are not aligned.",
         call. = FALSE)
  }
  y_hat <- mu_hat + tau_hat * treatment_coding
  result <- list(mu_hat = mu_hat, tau_hat = tau_hat, y_hat = y_hat)
  if (isTRUE(params$has_rfx)) {
    if (is.null(rfx_group_ids) || length(rfx_group_ids) != nrow(X)) {
      stop("One random-effect group ID is required per prediction row.", call. = FALSE)
    }
    groups <- factor(rfx_group_ids, levels = object$rfx_unique_group_ids)
    if (anyNA(groups)) {
      stop("Prediction random-effect groups must be present in the training data.",
           call. = FALSE)
    }
    if (is.null(rfx_basis)) {
      if (isTRUE(params$has_rfx_basis)) {
        stop("rfx_basis is required for this model.", call. = FALSE)
      }
      rfx_basis <- matrix(1, nrow(X), 1L)
    } else {
      rfx_basis <- as.matrix(rfx_basis)
    }
    if (!is.numeric(rfx_basis) || nrow(rfx_basis) != nrow(X) ||
        any(!is.finite(rfx_basis)) || ncol(rfx_basis) != params$num_rfx_basis) {
      stop("rfx_basis must match the training basis dimension and prediction rows.",
           call. = FALSE)
    }
    rfx_predictions <- object$rfx_samples$predict(as.integer(groups), rfx_basis) *
      params$outcome_scale
    if (!identical(dim(rfx_predictions), dim(y_hat))) {
      stop("Random-effect draws are not aligned with outcome draws.", call. = FALSE)
    }
    result$y_hat <- result$y_hat + rfx_predictions
    result$rfx_predictions <- rfx_predictions
  }
  result
}
