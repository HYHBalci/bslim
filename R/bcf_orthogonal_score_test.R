#' Experimental cross-fitted orthogonal test of mean treatment heterogeneity
#'
#' Tests fixed, prespecified directions of departure from a constant CATE using
#' out-of-fold marginal outcome and propensity regressions. Bayesian nuisance
#' estimation does not make the resulting frequentist test fully Bayesian.
#' @param X_train Numeric covariate matrix used for nuisance learning.
#' @param Z_train Binary numeric treatment vector.
#' @param y_train Numeric continuous outcome vector.
#' @param score_basis Prespecified numeric modifier matrix; defaults to X_train.
#' @param nuisance_method One of "bart", "glm", "supplied", or "custom".
#' @param m_hat,e_hat Supplied marginal-outcome and propensity predictions.
#' @param nuisance_source For supplied predictions, "cross_fitted" or "oracle".
#'   The latter explicitly asserts known nuisance functions, not in-sample fits.
#' @param known_propensity Known design propensity, scalar or length n. Never
#'   use an estimated full-sample propensity under this name.
#' @param fold_id Vector identifying held-out folds. Generated if absent.
#' @param num_folds Number of generated folds, assigned independently of the observations.
#' @param training_indices For supplied cross-fitted predictions, a list of
#'   training row indices in order of unique(fold_id). Each must exclude its
#'   held-out fold. This audits declared provenance, not the prediction algorithm.
#' @param nuisance_fitter Custom function with arguments X_train, Z_train,
#'   y_train, X_test, returning held-out m_hat and e_hat. Labels from held-out
#'   rows are never passed. e_hat may be omitted if known_propensity is supplied.
#' @param variance_correction "HC0" (default) or "HC3". HC3 divides score contributions by
#'   one minus leverage of the residualized treatment design, then recentres.
#'   It is a conservative sensitivity adjustment, not a validated small-sample improvement.
#' @param num_null_draws Gaussian multiplier draws, at least 99, or zero to skip.
#' @param seed Nonnegative integer seed. The caller's R RNG state is restored.
#' @param propensity_clip Optional explicit bound in (0,.5); no clipping by
#'   default. Clipping changes the nuisance and can introduce bias; rate claims require true propensities inside the clipping interval.
#' @param rank_action "error" or "drop" for redundant score-basis columns.
#' @param num_gfr,num_burnin,num_mcmc BART nuisance sampling budgets.
#' @param bart_params General BART nuisance parameters, overriding defaults.
#' @param mean_forest_params BART mean-forest nuisance parameters.
#' @return List of class bcforthogonalscore with p_value_quad (chi-square),
#'   p_value_quad_multiplier and p_value_max (plus-one Gaussian multiplier),
#'   statistics, alpha, scores, covariance, fold provenance and diagnostics.
#' @details Under the constant-CATE null, the population score is
#'   (H-c)(Z-e){Y-m-alpha(Z-e)}, with c=E[H(Z-e)^2]/E[(Z-e)^2]. Its derivatives
#'   in m, e and alpha vanish at their true values. Fixed-dimensional inference
#'   requires independent observations, overlap, nonsingular score covariance,
#'   moments and sufficiently fast nuisance convergence, including a propensity
#'   error squared term and the product of outcome and propensity errors being
#'   o(n^-1/2). This is not unrestricted double robustness. Fixed-score multiplier
#'   draws do not re-estimate nuisance functions or include their finite-sample
#'   uncertainty. Constant CATE implies the tested moments vanish; their converse
#'   need not hold. Rank and basis selection using outcomes is not supported.
#'   BART uses held-out predictions and a probit propensity model; reported
#'   propensity predictions average transformed latent draws.
#' @export
bcf_orthogonal_score_test <- function(
    X_train, Z_train, y_train, score_basis = NULL,
    nuisance_method = c("bart", "glm", "supplied", "custom"),
    m_hat = NULL, e_hat = NULL, nuisance_source = c("cross_fitted", "oracle"),
    known_propensity = NULL, fold_id = NULL, num_folds = 5L,
    training_indices = NULL, nuisance_fitter = NULL,
    variance_correction = c("HC0", "HC3"), num_null_draws = 1999L,
    seed = 1L, propensity_clip = NULL, rank_action = c("error", "drop"),
    num_gfr = 5L, num_burnin = 50L, num_mcmc = 100L,
    bart_params = list(), mean_forest_params = list()) {
  nuisance_method <- match.arg(nuisance_method)
  nuisance_source <- match.arg(nuisance_source)
  variance_correction <- match.arg(variance_correction)
  rank_action <- match.arg(rank_action)
  if (!is.list(bart_params) || !is.list(mean_forest_params)) stop("bart_params and mean_forest_params must be lists.")
  if (nuisance_method == "bart") {
    budgets <- c(num_gfr, num_burnin, num_mcmc)
    if (length(budgets) != 3L || any(!is.finite(budgets)) || any(budgets != floor(budgets)) || any(budgets < 0) || num_mcmc < 1L)
      stop("BART budgets must be nonnegative integers with num_mcmc >= 1 retained draw.")
  }
  numeric_matrix <- function(x, name) {
    x <- as.matrix(x)
    if (!is.numeric(x) || length(dim(x)) != 2L || any(!is.finite(x)))
      stop(name, " must be a finite numeric matrix (encode categorical variables explicitly).")
    x
  }
  X <- numeric_matrix(X_train, "X_train"); n <- nrow(X)
  if (n < 10L || ncol(X) < 1L) stop("At least ten observations and one covariate are required.")
  if (!is.numeric(y_train) || length(y_train) != n || any(!is.finite(y_train)))
    stop("y_train must be a finite numeric vector of length n.")
  if (!is.numeric(Z_train) || length(Z_train) != n || any(!Z_train %in% c(0, 1)) || length(unique(Z_train)) != 2L)
    stop("Z_train must contain both binary numeric treatment levels 0 and 1.")
  Y <- as.numeric(y_train); Z <- as.numeric(Z_train)
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 || seed != floor(seed) || seed > .Machine$integer.max - 10000L)
    stop("seed must be a nonnegative integer below the integer limit.")
  if (length(num_null_draws) != 1L || !is.finite(num_null_draws) || num_null_draws != floor(num_null_draws) ||
      !(num_null_draws == 0 || num_null_draws >= 99)) stop("num_null_draws must be zero or an integer at least 99.")
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(as.integer(seed))
  H <- numeric_matrix(if (is.null(score_basis)) X else score_basis, "score_basis")
  if (nrow(H) != n || ncol(H) < 1L) stop("score_basis must have n rows and at least one column.")
  if (is.null(colnames(H))) colnames(H) <- paste0("H", seq_len(ncol(H)))
  original_names <- colnames(H)
  H <- sweep(H, 2L, colMeans(H), "-")
  scale_h <- sqrt(colSums(H^2) / n)
  active <- which(scale_h > 0)
  if (!length(active)) stop("score_basis has no varying columns.")
  H1 <- sweep(H[, active, drop = FALSE], 2L, scale_h[active], "/")
  qr_h <- qr(H1, tol = 1e-9)
  keep <- active[sort(qr_h$pivot[seq_len(qr_h$rank)])]
  dropped <- setdiff(seq_len(ncol(H)), keep)
  if (length(dropped) && rank_action == "error") stop("score_basis is rank deficient after centering; remove constant/redundant columns or set rank_action='drop'.")
  H <- sweep(H[, keep, drop = FALSE], 2L, scale_h[keep], "/")
  q <- ncol(H)
  if (n <= q + 2L) stop("Insufficient observations for the score dimension.")
  if (q > n / 10) warning("Score dimension is large relative to n; asymptotic calibration can be inaccurate.")
  known_e <- NULL
  if (!is.null(known_propensity)) {
    known_e <- if (length(known_propensity) == 1L) rep(known_propensity, n) else known_propensity
    if (!is.numeric(known_e) || length(known_e) != n || any(!is.finite(known_e)) || any(known_e <= 0 | known_e >= 1))
      stop("known_propensity must lie strictly inside (0,1) and be scalar or length n.")
  }
  if (nuisance_method != "supplied" && (!is.null(m_hat) || !is.null(e_hat))) stop("Use nuisance_method='supplied' for external predictions.")
  need_folds <- nuisance_method != "supplied" || nuisance_source == "cross_fitted"
  if (need_folds) {
    if (is.null(fold_id)) {
      if (nuisance_method == "supplied") stop("Supplied cross-fitted predictions require fold_id and training_indices.")
      if (length(num_folds) != 1L || !is.finite(num_folds) || num_folds != floor(num_folds) || num_folds < 2L || num_folds > min(table(Z)))
        stop("num_folds must be between 2 and the smaller treatment-arm size.")
      fold_id <- sample(rep(seq_len(num_folds), length.out = n))

    }
    if (length(fold_id) != n || anyNA(fold_id) || length(unique(fold_id)) < 2L) stop("fold_id must identify at least two nonmissing held-out folds for every row.")
    levels_fold <- unique(fold_id)
    held_out <- lapply(levels_fold, function(k) which(fold_id == k))
    if (nuisance_method == "supplied") {
      if (!is.list(training_indices) || length(training_indices) != length(held_out)) stop("training_indices must contain one training-index vector per unique(fold_id), in that order.")
      for (k in seq_along(held_out)) {
        tr <- training_indices[[k]]
        if (!is.numeric(tr) || !length(tr) || anyNA(tr) || any(tr != floor(tr)) || any(tr < 1 | tr > n) || anyDuplicated(tr) || any(tr %in% held_out[[k]]))
          stop("training_indices contains invalid indices or held-out fold leakage.")
      }
    } else training_indices <- lapply(held_out, function(te) setdiff(seq_len(n), te))
  } else { held_out <- list(); training_indices <- list(); fold_id <- rep(NA_integer_, n) }
  if (nuisance_method == "supplied") {
    if (!is.null(known_e) && !is.null(e_hat)) stop("Supply either e_hat or known_propensity, not both.")
    if (!is.null(known_e)) e_hat <- known_e
  } else {
    m_hat <- e_hat <- rep(NA_real_, n)
    if (nuisance_method == "custom" && !is.function(nuisance_fitter)) stop("custom nuisances require nuisance_fitter.")
    for (k in seq_along(held_out)) {
      te <- held_out[[k]]; tr <- training_indices[[k]]
      if (length(unique(Z[tr])) != 2L) stop("Every nuisance training fold must contain both treatment levels.")
      if (nuisance_method == "glm") {
        design_tr <- cbind(1, X[tr, , drop = FALSE]); design_te <- cbind(1, X[te, , drop = FALSE])
        outcome_fit <- stats::lm.fit(design_tr, Y[tr])
        if (outcome_fit$rank != ncol(design_tr)) stop("GLM outcome training design is rank deficient; simplify nuisance covariates.")
        m_hat[te] <- drop(design_te %*% outcome_fit$coefficients)
        if (is.null(known_e)) {
          propensity_fit <- stats::glm.fit(design_tr, Z[tr], family = stats::binomial())
          if (!propensity_fit$converged || any(!is.finite(propensity_fit$coefficients))) stop("Propensity GLM did not converge with finite coefficients.")
          e_hat[te] <- stats::plogis(drop(design_te %*% propensity_fit$coefficients))
        }
      } else if (nuisance_method == "custom") {
        pred <- nuisance_fitter(X_train = X[tr, , drop = FALSE], Z_train = Z[tr], y_train = Y[tr], X_test = X[te, , drop = FALSE])
        if (!is.list(pred) || length(pred$m_hat) != length(te) || (is.null(known_e) && length(pred$e_hat) != length(te))) stop("nuisance_fitter must return held-out m_hat and e_hat vectors.")
        m_hat[te] <- pred$m_hat
        if (is.null(known_e)) e_hat[te] <- pred$e_hat
      } else {
        gp <- utils::modifyList(list(verbose = FALSE, random_seed = as.integer(seed + k), keep_gfr = FALSE), bart_params)
        gp$random_seed <- as.integer(seed + k)
        mp <- utils::modifyList(list(num_trees = 50L, probit_outcome_model = FALSE), mean_forest_params)
        mp$probit_outcome_model <- FALSE
        fo <- bart(X_train = X[tr, , drop = FALSE], y_train = Y[tr], X_test = X[te, , drop = FALSE], num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc, general_params = gp, mean_forest_params = mp)
        if (is.null(fo$y_hat_test) || ncol(as.matrix(fo$y_hat_test)) < 1L) stop("BART returned no retained held-out outcome draws.")
        m_hat[te] <- rowMeans(as.matrix(fo$y_hat_test))
        if (is.null(known_e)) {
          gp$random_seed <- as.integer(seed + 1000L + k); gp$sample_sigma2_global <- FALSE
          mp$probit_outcome_model <- TRUE
          fe <- bart(X_train = X[tr, , drop = FALSE], y_train = Z[tr], X_test = X[te, , drop = FALSE], num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc, general_params = gp, mean_forest_params = mp)
          if (is.null(fe$y_hat_test) || ncol(as.matrix(fe$y_hat_test)) < 1L) stop("BART returned no retained held-out propensity draws.")
          e_hat[te] <- rowMeans(stats::pnorm(as.matrix(fe$y_hat_test)))
        }
      }
    }
    if (!is.null(known_e)) e_hat <- known_e
  }
  if (!is.numeric(m_hat) || length(m_hat) != n || any(!is.finite(m_hat)) || !is.numeric(e_hat) || length(e_hat) != n || any(!is.finite(e_hat)))
    stop("Nuisance predictions must be finite numeric vectors of length n.")
  e_raw <- as.numeric(e_hat)
  if (!is.null(propensity_clip)) {
    if (length(propensity_clip) != 1L || !is.finite(propensity_clip) || propensity_clip <= 0 || propensity_clip >= .5) stop("propensity_clip must lie in (0,.5).")
    e_hat <- pmin(pmax(e_hat, propensity_clip), 1 - propensity_clip)
  }
  if (any(e_hat <= 0 | e_hat >= 1)) stop("Propensity predictions must lie strictly inside (0,1); inspect overlap or explicitly set propensity_clip.")
  clipped <- sum(e_raw != e_hat)
  if (clipped) warning(clipped, " propensity predictions were clipped; this can introduce bias.")
  W <- Z - e_hat; denominator <- sum(W^2)
  if (!is.finite(denominator) || denominator <= n * .Machine$double.eps) stop("Residual treatment variation is numerically zero.")
  alpha <- sum(W * (Y - m_hat)) / denominator
  projection <- colSums(H * W^2) / denominator
  H_perp <- sweep(H, 2L, projection, "-")
  residual <- Y - m_hat - alpha * W
  psi <- H_perp * (W * residual)
  U <- colSums(psi)
  design <- cbind(W, H_perp * W)
  qr_design <- qr(design, tol = 1e-9)
  if (qr_design$rank != q + 1L) stop("Residualized treatment-score design is rank deficient.")
  leverage <- rowSums(qr.Q(qr_design)^2)
  if (any(leverage >= 1 - 1e-8)) stop("Near-unit leverage prevents stable score covariance; simplify the basis or increase sample size.")
  adjusted <- if (variance_correction == "HC3") psi / (1 - leverage) else psi
  adjusted <- sweep(adjusted, 2L, colMeans(adjusted), "-")
  covariance <- crossprod(adjusted)
  eig <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eig)) || min(eig) <= max(eig) * 1e-10 || max(eig) <= 0) stop("Score covariance is singular or ill-conditioned; no ridge-calibrated p-value is returned.")
  inverse <- chol2inv(chol(covariance))
  statistic_quad <- drop(crossprod(U, inverse %*% U))
  standard_error <- sqrt(diag(covariance))
  statistic_max <- max(abs(U / standard_error))
  p_quad_multiplier <- p_max <- NA_real_
  if (num_null_draws > 0) {
    exceed_quad <- exceed_max <- 0L; remaining <- num_null_draws
    while (remaining > 0L) {
      b <- min(250L, remaining)
      star <- crossprod(adjusted, matrix(stats::rnorm(n * b), n, b))
      q_star <- colSums(star * (inverse %*% star))
      m_star <- apply(abs(star / standard_error), 2L, max)
      exceed_quad <- exceed_quad + sum(q_star >= statistic_quad)
      exceed_max <- exceed_max + sum(m_star >= statistic_max)
      remaining <- remaining - b
    }
    p_quad_multiplier <- (1 + exceed_quad) / (num_null_draws + 1)
    p_max <- (1 + exceed_max) / (num_null_draws + 1)
  }
  structure(list(p_value_quad = stats::pchisq(statistic_quad, df = q, lower.tail = FALSE),
    p_value_quad_multiplier = p_quad_multiplier, p_value_max = p_max,
    statistic_quad = statistic_quad, statistic_max = statistic_max, df = q,
    alpha = alpha, score_sum = U, covariance = covariance, scores = psi,
    multiplier_scores = adjusted, basis = H, projected_basis = H_perp,
    alpha_projection = projection, fold_id = fold_id, training_indices = training_indices,
    nuisance_predictions = list(m_hat = as.numeric(m_hat), e_hat = as.numeric(e_hat)),
    diagnostics = list(n = n, q = q, nuisance_method = nuisance_method,
      nuisance_source = if (nuisance_method == "supplied") nuisance_source else "cross_fitted",
      known_propensity = !is.null(known_e), variance_correction = variance_correction,
      propensity_range = range(e_hat), raw_propensity_range = range(e_raw), clipped = clipped, max_leverage = max(leverage),
      residual_treatment_effective_n = denominator^2 / sum(W^4),
      covariance_condition = max(eig) / min(eig), dropped_columns = original_names[dropped],
      retained_columns = original_names[keep], num_null_draws = num_null_draws,
      nuisance_uncertainty_reestimated = FALSE), call = match.call()), class = "bcforthogonalscore")
}
