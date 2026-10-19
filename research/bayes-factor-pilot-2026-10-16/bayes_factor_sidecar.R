# Opt-in research sidecar. This file is not part of the package namespace.
# For the model, identities, and limitations see README.md in this directory.

log_sum_exp <- function(x) {
  a <- max(x)
  if (!is.finite(a)) return(a)
  a + log(sum(exp(x - a)))
}

log_mean_exp <- function(x) log_sum_exp(x) - log(length(x))

standardize_modifier_basis <- function(H) {
  H <- as.matrix(H)
  if (!is.numeric(H) || any(!is.finite(H))) stop("H must be finite numeric.")
  H <- sweep(H, 2L, colMeans(H), "-")
  scales <- sqrt(colMeans(H^2))
  if (any(scales <= 0)) stop("Constant modifiers are not supported.")
  sweep(H, 2L, scales, "/")
}

# Apply (sigma2 I + alpha_variance W W')^{-1} without making an n by n matrix.
apply_null_precision <- function(A, W, sigma2 = 1, alpha_variance = 4) {
  A <- as.matrix(A)
  if (length(W) != nrow(A) || !is.finite(sigma2) || sigma2 <= 0 ||
      !is.finite(alpha_variance) || alpha_variance < 0) stop("Invalid null covariance.")
  (A - W %o% (alpha_variance * drop(crossprod(W, A)) /
                 (sigma2 + alpha_variance * sum(W^2)))) / sigma2
}

# Each column of residual is y-f from one null-posterior forest draw.
# Alpha is analytically integrated under its proper, shared N(0,va) prior.
# Delta|S ~ N(0, slab_sd^2 I). Returns log(m_S(y|f)/m_0(y|f)).
conditional_log_bf <- function(residual, D, W, slab_sd = .5,
                               sigma2 = 1, alpha_variance = 4) {
  residual <- as.matrix(residual); D <- as.matrix(D)
  if (nrow(D) != nrow(residual) || length(W) != nrow(D) ||
      any(!is.finite(c(residual, D, W))) || length(slab_sd) != 1L ||
      !is.finite(slab_sd) || slab_sd < 0) stop("Invalid BF inputs.")
  if (slab_sd == 0 || ncol(D) == 0L) return(rep(0, ncol(residual)))
  precision_D <- apply_null_precision(D, W, sigma2, alpha_variance)
  A <- diag(ncol(D)) + slab_sd^2 * crossprod(D, precision_D)
  L <- chol(A)
  score <- crossprod(precision_D, residual)
  solved <- forwardsolve(t(L), score)
  -sum(log(diag(L))) + .5 * slab_sd^2 * colSums(solved^2)
}

subset_prior <- function(q, inclusion_probability = 1 / q) {
  if (q < 1 || q > 12 || inclusion_probability <= 0 || inclusion_probability >= 1)
    stop("Use 1 <= q <= 12 and 0 < inclusion_probability < 1.")
  subsets <- lapply(seq_len(2^q - 1L), function(mask)
    which(as.logical(intToBits(mask)[seq_len(q)])))
  size <- lengths(subsets)
  log_weight <- size * log(inclusion_probability) +
    (q - size) * log1p(-inclusion_probability)
  log_weight <- log_weight - log_sum_exp(log_weight)
  list(subsets = subsets, log_weight = log_weight)
}

evaluate_bf_priors <- function(residual, D, W, slab_sd = .5,
                               sigma2 = 1, alpha_variance = 4) {
  D <- as.matrix(D); residual <- as.matrix(residual)
  dense <- conditional_log_bf(residual, D, W, slab_sd, sigma2, alpha_variance)
  p <- subset_prior(ncol(D), inclusion_probability = min(.5, 1/ncol(D)))
  lbf <- vapply(p$subsets, function(s)
    conditional_log_bf(residual, D[, s, drop = FALSE], W,
                       slab_sd, sigma2, alpha_variance), numeric(ncol(residual)))
  lbf <- matrix(lbf, nrow = ncol(residual), ncol = length(p$subsets))
  weighted <- sweep(lbf, 2L, p$log_weight, "+")
  sparse <- apply(weighted, 1L, log_sum_exp)
  list(dense = dense, sparse = sparse)
}

importance_summary <- function(log_ratios, chain_id = rep(1L, length(log_ratios))) {
  if (any(!is.finite(log_ratios))) stop("Non-finite likelihood ratio.")
  w <- exp(log_ratios - max(log_ratios)); wn <- w / sum(w)
  by_chain <- vapply(split(log_ratios, chain_id), log_mean_exp, numeric(1))
  blocks <- unlist(lapply(split(log_ratios, chain_id), function(x) {
    id <- cut(seq_along(x), breaks = min(10L, floor(length(x)/10L)), labels = FALSE)
    vapply(split(x, id), log_mean_exp, numeric(1))
  }), use.names = FALSE)
  # Batch-means relative SE is descriptive and unreliable if important tails are missed.
  bw <- exp(blocks - max(log_ratios))
  relative_se <- if (length(bw) > 1) sd(bw) / sqrt(length(bw)) / mean(w) else NA_real_
  list(log_bf = log_mean_exp(log_ratios),
       importance_ess = 1 / sum(wn^2),
       importance_ess_fraction = 1 / (length(wn) * sum(wn^2)),
       largest_weight_fraction = max(wn),
       top_one_percent_weight_fraction = sum(sort(wn, decreasing = TRUE)[
         seq_len(max(1, ceiling(.01 * length(wn))))]),
       chain_log_bf_range = diff(range(by_chain)),
       batch_relative_se = relative_se,
       chain_log_bf = by_chain,
       warning = if (max(wn) > .1 || diff(range(by_chain)) > log(3) ||
                     1/sum(wn^2) < 100) "importance_concentration_or_chain_disagreement" else "none")
}

restricted_moment_test <- function(y, m_hat, alpha_hat, W, H) {
  H <- standardize_modifier_basis(H)
  # This is the external centered-contribution endpoint described in the draft.
  c <- W * (y - m_hat - alpha_hat * W)
  c <- c - mean(c)
  psi <- H * c
  U <- colSums(psi); V <- crossprod(psi)
  Q <- drop(crossprod(U, chol2inv(chol(V)) %*% U))
  list(Q = Q, df = ncol(H), p = pchisq(Q, ncol(H), lower.tail = FALSE))
}

# Exact in the stipulated Gaussian score experiment only, approximate for a
# data-derived sandwich score. This is not a forest marginal likelihood BF.
score_experiment_log_bf <- function(Q, q, g) {
  if (any(Q < 0) || q < 1 || any(g < 0)) stop("Invalid score experiment inputs.")
  -.5 * q * log1p(g) + .5 * g / (1 + g) * Q
}

wilson_interval <- function(k, n, level = .95) {
  if (!n) return(c(NA_real_, NA_real_))
  z <- qnorm((1 + level)/2); p <- k/n; den <- 1 + z^2/n
  mid <- (p + z^2/(2*n))/den
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2))/den
  c(max(0, mid-half), min(1, mid+half))
}
