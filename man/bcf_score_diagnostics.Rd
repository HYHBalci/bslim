\name{bcf_score_diagnostics}
\alias{bcf_score_diagnostics}
\title{Diagnose pooled restricted-BCF mean scores and Monte Carlo precision}
\usage{
bcf_score_diagnostics(fit, y, z, propensity, H, batch_size = NULL,
  alpha = 0.05, rank_tol = 1e-10, ridge = 1e-8, keep_draws = FALSE)
}
\arguments{
\item{fit}{A returned posterior-averaged mean-test bcfscoretest object.}
\item{y, z}{Original numeric outcomes and binary treatment assignments.}
\item{propensity}{The actual propensity marks used by the fit, scalar or length n. Internally fitted marks must be retained externally.}
\item{H}{A prespecified finite n-row score matrix. Columns are recentered. The exact prepared package score matrix is required to reconstruct per-chain ridge p-values.}
\item{batch_size}{Contiguous draws per batch; defaults to floor(sqrt(retained draws)). At least eight complete batches per chain are required.}
\item{alpha}{Decision threshold for descriptive flags.}
\item{rank_tol}{Relative eigenvalue tolerance for covariance rank.}
\item{ridge}{Fixed nonnegative ridge in sampler units; default 1e-8.}
\item{keep_draws}{Whether to return linearized score/covariance trajectories.}
}
\value{A list with pooled U,V,quadratic statistic and p-value, per-chain tests, unregularized covariance rank and condition, batch Monte Carlo covariance, delta MCSE of statistic and p-value, and optional posterior Rhat/ESS.}
\description{Reconstructs a score after pooling retained chain contributions and provides covariance-rank and score-specific Monte Carlo diagnostics. The fitted sampler and legacy per-chain p-values are unchanged.}
\details{Returned draws are already on original outcome units. Contributions are pooled before covariance is computed; averaging p-values is a different procedure. Per-draw residual standardization is reconstructed when enabled. Only the Gaussian quadratic mean-test reference is supported. A matching H is required for package p-value agreement. Singular or ill-conditioned unregularized covariance produces an NA primary p-value even if a separate legacy ridge number is computable.

Delta Monte Carlo error includes score and covariance changes. Batch estimates assume stationary independent chains and adequate batch length; they exclude burn-in bias and between-chain nonstationarity. The suggested posterior package supplies rank-normalized Rhat and ESS if installed. Numerical diagnostics do not establish frequentist calibration, nuisance convergence, same-sample inferential guarantees or observational validity.}