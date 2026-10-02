#' Experimental cross-fitted weighted AIPW test of mean treatment heterogeneity
#'
#' Tests fixed, prespecified CATE directions using arm-specific outcome means
#' and propensity predictions. This combines nuisance learning with frequentist
#' score inference; it is not fully Bayesian inference.
#' @param X_train Finite numeric matrix for nuisance learning.
#' @param Z_train Binary numeric treatment vector containing both arms.
#' @param y_train Finite numeric continuous outcome vector.
#' @param score_basis Prespecified finite numeric modifier matrix; defaults to X_train.
#' @param nuisance_method One of "bart", "glm", "supplied", or "custom".
#' @param m0_hat,m1_hat,e_hat Supplied control/treated outcome means and propensity.
#' @param nuisance_source For supplied inputs, "cross_fitted" or "oracle", either
#'   common to all nuisances or a named vector/list with keys m0, m1, e.
#'   Oracle means the true function, not a full-sample fitted prediction.
#' @param known_propensity Known design propensity, scalar or length n.
#' @param fold_id Common held-out fold vector, or (for supplied inputs) a named
#'   list with a separate vector for each cross-fitted nuisance.
#' @param num_folds Number of generated, observation-independent balanced folds.
#' @param training_indices Common list in order unique(fold_id), or named lists
#'   for each supplied cross-fitted nuisance. Training rows must exclude that
#'   nuisance's held-out fold; arm-specific subsets are allowed.
#' @param nuisance_fitter Custom function accepting X_train, Z_train, y_train,
#'   X_test and returning held-out m0_hat, m1_hat, and e_hat (unless known).
#'   It never receives held-out outcomes or treatments.
#' @param weighting "overlap" (default) uses e(1-e); "uniform" uses unit weights.
#' @param variance_correction "HC0" (default) or "HC3", a leverage heuristic.
#' @param num_null_draws At least 99 Gaussian multiplier draws, or zero to skip.
#' @param seed Nonnegative integer seed. The caller's R RNG state is restored.
#' @param propensity_clip Optional explicit clipping bound in (0,.5).
#' @param rank_action Error on redundant basis columns or explicitly drop them.
#' @param num_gfr,num_burnin,num_mcmc BART budgets with at least one retained draw.
#' @param bart_params General BART parameters. Seeds and warm-up retention are controlled here.
#' @param mean_forest_params BART mean-forest parameters; outcome/probit type is controlled here.
#' @return A bcfaipwscore list with quadratic and maximum statistics, asymptotic
#'   and multiplier p-values, weighted common effect, scores, covariance,
#'   nuisance predictions and per-nuisance provenance.
#' @details The uniform AIPW pseudo-outcome is m1-m0+Z(Y-m1)/e-(1-Z)(Y-m0)/(1-e).
#'   With weight omega, let t=omega times this pseudo-outcome. Estimate alpha by
#'   sum(t)/sum(omega), c by sum(omega H)/sum(omega), and use contributions
#'   (H-c)(t-alpha omega). HC0 centres the contributions for covariance. HC3
#'   first divides by one minus leverage of sqrt(omega)[1,H-c], then centres;
#'   it is an unvalidated sensitivity heuristic, not regression HC3 inference.
#'   Constant CATE implies zero tested moments, but a finite basis can miss
#'   heterogeneous alternatives. AIPW conditional means have a double-robust
#'   identification property; valid inference still requires overlap, fixed
#'   dimension, moments, nonsingular covariance and suitable nuisance-rate and
#'   product-error control. No BART nuisance convergence rate is established
#'   by this implementation. Clipping can alter these conditions. Supplied
#'   provenance is checked as a declaration, not proof of actual training honesty. Separate nuisance partitions require suitable common-refinement and rate assumptions.
#'   Multiplier draws keep nuisances fixed and omit finite-sample refitting uncertainty.
#' @export
bcf_aipw_score_test <- function(
    X_train, Z_train, y_train, score_basis = NULL,
    nuisance_method = c("bart", "glm", "supplied", "custom"),
    m0_hat = NULL, m1_hat = NULL, e_hat = NULL,
    nuisance_source = "cross_fitted", known_propensity = NULL,
    fold_id = NULL, num_folds = 5L, training_indices = NULL,
    nuisance_fitter = NULL, weighting = c("overlap", "uniform"),
    variance_correction = c("HC0", "HC3"), num_null_draws = 1999L,
    seed = 1L, propensity_clip = NULL, rank_action = c("error", "drop"),
    num_gfr = 5L, num_burnin = 50L, num_mcmc = 100L,
    bart_params = list(), mean_forest_params = list()) {
  nuisance_method <- match.arg(nuisance_method)
  weighting <- match.arg(weighting)
  variance_correction <- match.arg(variance_correction)
  rank_action <- match.arg(rank_action)
  if (!is.list(bart_params) || !is.list(mean_forest_params)) stop("BART parameter arguments must be lists.")
  scalar_integer <- function(x, minimum = 0) is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x) && x >= minimum && x <= .Machine$integer.max
  if (!scalar_integer(seed)) stop("seed must be a nonnegative integer within the integer range.")
  if (!scalar_integer(num_null_draws) || !(num_null_draws == 0 || num_null_draws >= 99)) stop("num_null_draws must be zero or an integer at least 99.")
  if (nuisance_method == "bart" && !(scalar_integer(num_gfr) && scalar_integer(num_burnin) && scalar_integer(num_mcmc, 1))) stop("BART budgets require nonnegative integers and at least one retained draw.")
  matrix_input <- function(x, label) {
    x <- as.matrix(x)
    if (!is.numeric(x) || length(dim(x)) != 2L || any(!is.finite(x))) stop(label, " must be a finite numeric matrix.")
    x
  }
  vector_input <- function(x, size, label) {
    if (!is.numeric(x) || !is.null(dim(x)) || length(x) != size || any(!is.finite(x))) stop(label, " must be a finite numeric vector of the required length.")
    as.numeric(x)
  }
  X <- matrix_input(X_train, "X_train"); n <- nrow(X)
  if (n < 10L || ncol(X) < 1L) stop("At least ten observations and one covariate are required.")
  Y <- vector_input(y_train, n, "y_train"); Z <- vector_input(Z_train, n, "Z_train")
  if (!all(Z %in% c(0, 1)) || length(unique(Z)) != 2L) stop("Z_train must contain both binary treatment levels.")
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(as.integer(seed))
  learner_seed <- function(k, component) as.integer((as.double(seed) + 100 * k + component) %% .Machine$integer.max)
  H <- matrix_input(if (is.null(score_basis)) X else score_basis, "score_basis")
  if (nrow(H) != n || ncol(H) < 1L) stop("score_basis must contain n rows and at least one column.")
  if (is.null(colnames(H))) colnames(H) <- paste0("H", seq_len(ncol(H)))
  original_names <- colnames(H)
  centered <- sweep(H, 2L, colMeans(H), "-")
  scales <- sqrt(colMeans(centered^2)); active <- which(scales > 0)
  if (!length(active)) stop("score_basis has no varying columns.")
  rank_qr <- qr(sweep(centered[, active, drop = FALSE], 2L, scales[active], "/"), tol = 1e-9)
  keep <- active[sort(rank_qr$pivot[seq_len(rank_qr$rank)])]
  dropped <- setdiff(seq_len(ncol(H)), keep)
  if (length(dropped) && rank_action == "error") stop("score_basis is rank deficient after centering.")
  H <- H[, keep, drop = FALSE]; q <- ncol(H)
  if (n <= q + 2L) stop("Insufficient observations for the score dimension.")
  if (q > n / 10) warning("Score dimension is large relative to n; asymptotic calibration can be inaccurate.")
  known_e <- NULL
  if (!is.null(known_propensity)) {
    known_e <- vector_input(if (length(known_propensity) == 1L) rep(known_propensity, n) else known_propensity, n, "known_propensity")
    if (any(known_e <= 0 | known_e >= 1)) stop("known_propensity must lie strictly inside (0,1).")
    if (!is.null(e_hat)) stop("Supply either e_hat or known_propensity, not both.")
  }
  keys <- c("m0", "m1", "e")
  sources <- stats::setNames(rep("cross_fitted", 3L), keys)
  if (is.list(nuisance_source)) nuisance_source <- unlist(nuisance_source, use.names = TRUE)
  if (!is.character(nuisance_source) || anyNA(nuisance_source) || !length(nuisance_source) || any(!nuisance_source %in% c("cross_fitted", "oracle"))) stop("nuisance_source entries must be cross_fitted or oracle.")
  if (length(nuisance_source) == 1L && is.null(names(nuisance_source))) sources[] <- nuisance_source else {
    required <- if (is.null(known_e)) keys else keys[1:2]
    if (is.null(names(nuisance_source)) || anyDuplicated(names(nuisance_source)) || any(!names(nuisance_source) %in% keys) || !all(required %in% names(nuisance_source))) stop("Named nuisance_source must declare m0, m1 and, unless known, e.")
    sources[names(nuisance_source)] <- nuisance_source
  }
  if (!is.null(known_e)) sources["e"] <- "known"
  provenance <- stats::setNames(vector("list", 3L), keys)
  validate_folds <- function(folds) {
    if (!is.atomic(folds) || !is.null(dim(folds)) || length(folds) != n || anyNA(folds) || length(unique(folds)) < 2L) stop("Each cross-fitted nuisance requires a fold_id vector with at least two nonmissing folds.")
    folds
  }
  validate_training <- function(folds, indices, key) {
    held <- lapply(unique(folds), function(k) which(folds == k))
    if (!is.list(indices) || length(indices) != length(held)) stop(key, " training_indices must contain one vector per unique(fold_id), in that order.")
    for (k in seq_along(held)) {
      tr <- indices[[k]]
      if (!is.numeric(tr) || !is.null(dim(tr)) || !length(tr) || anyNA(tr) || any(!is.finite(tr)) || any(tr != floor(tr)) || any(tr < 1 | tr > n) || anyDuplicated(tr) || any(tr %in% held[[k]])) stop(key, " training_indices contain invalid rows or held-out leakage.")
    }
    indices
  }
  if (nuisance_method == "supplied") {
    nested_training <- is.list(training_indices) && any(vapply(training_indices, is.list, logical(1)))
    for (key in keys) {
      if (sources[key] == "cross_fitted") {
        folds <- validate_folds(if (is.list(fold_id)) fold_id[[key]] else fold_id)
        tr <- validate_training(folds, if (nested_training) training_indices[[key]] else training_indices, key)
        provenance[[key]] <- list(source = "cross_fitted", learner = "supplied", fold_id = folds, training_indices = tr)
      } else provenance[[key]] <- list(source = unname(sources[key]), learner = "supplied", fold_id = NULL, training_indices = list())
    }
    if (!is.null(known_e)) e_hat <- known_e
  } else {
    if (!is.null(m0_hat) || !is.null(m1_hat) || !is.null(e_hat)) stop("Use nuisance_method='supplied' for external predictions.")
    if (any(sources[sources != "known"] != "cross_fitted")) stop("Oracle declarations apply only to supplied nuisance predictions.")
    if (!is.null(training_indices)) stop("training_indices are generated internally; supply them only with supplied predictions.")
    if (is.list(fold_id)) stop("Internally fitted nuisances require one common fold_id vector.")
    if (is.null(fold_id)) {
      if (!scalar_integer(num_folds, 2) || num_folds > n) stop("num_folds must be an integer between 2 and n.")
      fold_id <- sample(rep(seq_len(num_folds), length.out = n))
    }
    fold_id <- validate_folds(fold_id)
    held <- lapply(unique(fold_id), function(k) which(fold_id == k))
    full_training <- lapply(held, function(te) setdiff(seq_len(n), te))
    for (key in keys) provenance[[key]] <- list(source = unname(sources[key]), learner = if (sources[key] == "known") "design" else nuisance_method, fold_id = if (sources[key] == "known") NULL else fold_id, training_indices = list(), seeds = integer())
    if (nuisance_method == "custom" && !is.function(nuisance_fitter)) stop("custom nuisances require nuisance_fitter.")
    m0_hat <- m1_hat <- e_hat <- rep(NA_real_, n)
    bart_prediction <- function(tr, te, response, probit, fit_seed) {
      if (length(tr) < 2L) stop("Each BART armwise training sample must contain at least two observations.")
      set.seed(fit_seed)
      gp <- utils::modifyList(list(verbose = FALSE), bart_params)
      gp$random_seed <- as.numeric(fit_seed); gp$keep_gfr <- FALSE; gp$keep_burnin <- FALSE
      mp <- utils::modifyList(list(num_trees = 50L), mean_forest_params)
      mp$probit_outcome_model <- probit
      if (probit) gp$sample_sigma2_global <- FALSE
      fit <- bart(X_train = X[tr, , drop = FALSE], y_train = response, X_test = X[te, , drop = FALSE], num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc, general_params = gp, mean_forest_params = mp)
      draws <- as.matrix(fit$y_hat_test)
      if (!is.numeric(draws) || nrow(draws) != length(te) || ncol(draws) < 1L || any(!is.finite(draws))) stop("BART returned invalid or missing retained held-out draws.")
      rowMeans(if (probit) stats::pnorm(draws) else draws)
    }
    for (k in seq_along(held)) {
      te <- held[[k]]; tr <- full_training[[k]]; tr0 <- tr[Z[tr] == 0]; tr1 <- tr[Z[tr] == 1]
      if (!length(tr0) || !length(tr1)) stop("Every nuisance training fold must contain both treatment arms; generated folds are not outcome-stratified or repaired.")
      if (nuisance_method == "custom") {
        fit_seed <- learner_seed(k, 3); set.seed(fit_seed)
        pred <- nuisance_fitter(X_train = X[tr, , drop = FALSE], Z_train = Z[tr], y_train = Y[tr], X_test = X[te, , drop = FALSE])
        if (!is.list(pred)) stop("nuisance_fitter must return a list of held-out predictions.")
        m0_hat[te] <- vector_input(pred$m0_hat, length(te), "custom m0_hat")
        m1_hat[te] <- vector_input(pred$m1_hat, length(te), "custom m1_hat")
        if (is.null(known_e)) e_hat[te] <- vector_input(pred$e_hat, length(te), "custom e_hat")
        for (key in keys[sources != "known"]) { provenance[[key]]$training_indices[[k]] <- tr; provenance[[key]]$seeds[k] <- fit_seed }
      } else {
        provenance$m0$training_indices[[k]] <- tr0; provenance$m1$training_indices[[k]] <- tr1
        if (is.null(known_e)) provenance$e$training_indices[[k]] <- tr
        if (nuisance_method == "glm") {
          des_te <- cbind(1, X[te, , drop = FALSE])
          arm_prediction <- function(indices) {
            des <- cbind(1, X[indices, , drop = FALSE]); fit <- stats::lm.fit(des, Y[indices])
            if (fit$rank != ncol(des) || any(!is.finite(fit$coefficients))) stop("GLM armwise outcome training design is rank deficient.")
            drop(des_te %*% fit$coefficients)
          }
          m0_hat[te] <- arm_prediction(tr0); m1_hat[te] <- arm_prediction(tr1)
          if (is.null(known_e)) {
            des <- cbind(1, X[tr, , drop = FALSE]); fit <- stats::glm.fit(des, Z[tr], family = stats::binomial())
            if (!fit$converged || fit$rank != ncol(des) || any(!is.finite(fit$coefficients))) stop("Propensity GLM failed to converge with a full-rank finite fit.")
            e_hat[te] <- stats::plogis(drop(des_te %*% fit$coefficients))
          }
        } else {
          seeds <- vapply(seq_len(3L), function(j) learner_seed(k, j), integer(1))
          m0_hat[te] <- bart_prediction(tr0, te, Y[tr0], FALSE, seeds[1])
          m1_hat[te] <- bart_prediction(tr1, te, Y[tr1], FALSE, seeds[2])
          provenance$m0$seeds[k] <- seeds[1]; provenance$m1$seeds[k] <- seeds[2]
          if (is.null(known_e)) { e_hat[te] <- bart_prediction(tr, te, Z[tr], TRUE, seeds[3]); provenance$e$seeds[k] <- seeds[3] }
        }
      }
    }
    if (!is.null(known_e)) e_hat <- known_e
  }
  m0_hat <- vector_input(m0_hat, n, "m0_hat"); m1_hat <- vector_input(m1_hat, n, "m1_hat"); e_raw <- vector_input(e_hat, n, "e_hat")
  if (any(e_raw < 0 | e_raw > 1)) stop("Propensity predictions must be probabilities in [0,1], even when clipping is requested.")
  e_hat <- e_raw
  if (!is.null(propensity_clip)) {
    if (!is.numeric(propensity_clip) || length(propensity_clip) != 1L || !is.finite(propensity_clip) || propensity_clip <= 0 || propensity_clip >= .5) stop("propensity_clip must lie in (0,.5).")
    e_hat <- pmin(pmax(e_hat, propensity_clip), 1 - propensity_clip)
  }
  if (any(e_hat <= 0 | e_hat >= 1)) stop("Propensities must lie strictly inside (0,1); clipping is never implicit.")
  clipped <- sum(e_raw != e_hat)
  if (clipped) warning(clipped, " propensity predictions were explicitly clipped; this may introduce bias.")
  provenance$e$prediction_transform <- list(propensity_clip = propensity_clip, clipped_rows = which(e_raw != e_hat))
  if (!is.null(known_e)) provenance$e$learner <- "design"
  if (weighting == "overlap") {
    omega <- e_hat * (1 - e_hat)
    t <- omega * (m1_hat - m0_hat) + Z * (1 - e_hat) * (Y - m1_hat) - (1 - Z) * e_hat * (Y - m0_hat)
  } else {
    omega <- rep(1, n)
    t <- m1_hat - m0_hat + Z / e_hat * (Y - m1_hat) - (1 - Z) / (1 - e_hat) * (Y - m0_hat)
  }
  if (any(!is.finite(t)) || sum(omega) <= n * .Machine$double.eps) stop("Invalid weighted pseudo-outcome or negligible total weight.")
  alpha <- sum(t) / sum(omega)
  projection <- colSums(H * omega) / sum(omega); B <- sweep(H, 2L, projection, "-")
  psi <- B * (t - alpha * omega); U <- colSums(psi)
  design <- cbind(1, B) * sqrt(omega); qd <- qr(design, tol = 1e-9)
  if (qd$rank != q + 1L) stop("Weighted leverage design is rank deficient.")
  leverage <- rowSums(qr.Q(qd)^2)
  if (any(leverage >= 1 - 1e-8)) stop("Near-unit leverage prevents stable covariance.")
  adjusted <- if (variance_correction == "HC3") psi / (1 - leverage) else psi
  adjusted <- sweep(adjusted, 2L, colMeans(adjusted), "-"); V <- crossprod(adjusted)
  eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eigenvalues)) || max(eigenvalues) <= 0 || min(eigenvalues) <= max(eigenvalues) * 1e-10) stop("Score covariance is singular or ill-conditioned; no ridge rescue is applied.")
  inverse <- chol2inv(chol(V)); statistic_quad <- drop(crossprod(U, inverse %*% U))
  se <- sqrt(diag(V)); statistic_max <- max(abs(U / se)); p_quad_multiplier <- p_max <- NA_real_
  if (num_null_draws > 0) {
    eq <- em <- 0L; remaining <- num_null_draws
    while (remaining > 0L) {
      count <- min(250L, remaining); star <- crossprod(adjusted, matrix(stats::rnorm(n * count), n, count))
      eq <- eq + sum(colSums(star * (inverse %*% star)) >= statistic_quad)
      em <- em + sum(apply(abs(star / se), 2L, max) >= statistic_max)
      remaining <- remaining - count
    }
    p_quad_multiplier <- (1 + eq) / (num_null_draws + 1); p_max <- (1 + em) / (num_null_draws + 1)
  }
  structure(list(
    p_value_quad = stats::pchisq(statistic_quad, q, lower.tail = FALSE),
    p_value_quad_multiplier = p_quad_multiplier, p_value_max = p_max,
    statistic_quad = statistic_quad, statistic_max = statistic_max, df = q,
    alpha = alpha, projection = projection, alpha_projection = projection,
    score_sum = U, covariance = V, scores = psi, multiplier_scores = adjusted,
    basis = H, projected_basis = B, t = t, omega = omega, fold_id = fold_id,
    training_indices = lapply(provenance, `[[`, "training_indices"),
    nuisance_provenance = provenance,
    nuisance_predictions = list(m0_hat = m0_hat, m1_hat = m1_hat,
      e_hat = e_hat, raw_e_hat = e_raw),
    diagnostics = list(
      n = n, q = q, nuisance_method = nuisance_method, nuisance_source = sources,
      known_propensity = !is.null(known_e),
      known_propensity_modified = !is.null(known_e) && clipped > 0,
      used_propensity_source = if (clipped) paste0("clipped_", sources[["e"]]) else sources[["e"]],
      weighting = weighting, variance_correction = variance_correction,
      propensity_range = range(e_hat), raw_propensity_range = range(e_raw),
      clipped = clipped, max_leverage = max(leverage),
      weight_effective_n = sum(omega)^2 / sum(omega^2),
      covariance_condition = max(eigenvalues) / min(eigenvalues),
      dropped_columns = original_names[dropped], retained_columns = original_names[keep],
      num_null_draws = num_null_draws, nuisance_uncertainty_reestimated = FALSE),
    call = match.call()), class = "bcfaipwscore")
}