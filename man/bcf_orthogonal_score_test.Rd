\name{bcf_orthogonal_score_test}
\alias{bcf_orthogonal_score_test}
\title{Experimental Cross-Fitted Orthogonal Test of Mean Treatment Heterogeneity}
\usage{
bcf_orthogonal_score_test(
  X_train, Z_train, y_train, score_basis = NULL,
  nuisance_method = c("bart", "glm", "supplied", "custom"),
  m_hat = NULL, e_hat = NULL,
  nuisance_source = c("cross_fitted", "oracle"),
  known_propensity = NULL, fold_id = NULL, num_folds = 5L,
  training_indices = NULL, nuisance_fitter = NULL,
  variance_correction = c("HC0", "HC3"), num_null_draws = 1999L,
  seed = 1L, propensity_clip = NULL, rank_action = c("error", "drop"),
  num_gfr = 5L, num_burnin = 50L, num_mcmc = 100L,
  bart_params = list(), mean_forest_params = list()
)
}
\arguments{
\item{X_train}{Finite numeric covariate matrix for nuisance learning. Encode categorical variables explicitly.}
\item{Z_train}{Binary numeric treatment vector containing both zero and one.}
\item{y_train}{Finite numeric continuous outcome vector.}
\item{score_basis}{Prespecified finite numeric modifier matrix; defaults to X_train. Nuisance covariates and test directions may differ.}
\item{nuisance_method}{BART, linear/logistic GLM, externally supplied predictions, or a custom learner.}
\item{m_hat, e_hat}{Externally supplied predictions of the marginal outcome mean and propensity.}
\item{nuisance_source}{For supplied predictions, cross-fitted or oracle. Oracle explicitly asserts true nuisance functions, not in-sample estimates.}
\item{known_propensity}{Known design propensity, scalar or length n. Do not pass an estimated full-sample propensity under this name.}
\item{fold_id}{Held-out fold identifier for every observation. Generated if omitted for internally fitted nuisances.}
\item{num_folds}{Number of generated folds, assigned independently of the observations.}
\item{training_indices}{For supplied cross-fitted predictions, a list of training row indices in the order of unique(fold_id). Each must exclude the corresponding held-out fold. The function checks declared provenance, not actual model honesty.}
\item{nuisance_fitter}{Custom function accepting X_train, Z_train, y_train and X_test; returns held-out m_hat and e_hat vectors. Test labels are not passed. e_hat may be omitted when known_propensity is supplied.}
\item{variance_correction}{HC0 (default) is the canonical empirical sandwich. HC3 is a conservative sensitivity adjustment using leverage of the full residualized treatment design, divides score contributions by one minus leverage, and recentres. HC0 uses centred unadjusted contributions. Neither gives a finite-sample guarantee.}
\item{num_null_draws}{At least 99 Gaussian multiplier draws, or zero to omit simulated p-values.}
\item{seed}{Nonnegative integer seed. Caller R RNG state is restored.}
\item{propensity_clip}{Optional explicit bound in (0,.5). No clipping by default. Clipping is reported and can introduce bias; rate arguments require true propensities inside the chosen interval.}
\item{rank_action}{Error on redundant modifier columns, or explicitly drop them.}
\item{num_gfr,num_burnin,num_mcmc}{BART budgets; at least one retained MCMC draw is required.}
\item{bart_params}{General BART nuisance parameters.}
\item{mean_forest_params}{BART mean-forest nuisance parameters.}
}
\description{Tests prespecified finite-dimensional directions against constant conditional average treatment effect, with cross-fitted marginal outcome and propensity regressions and a common-effect nuisance projection. Bayesian nuisance estimation is combined with frequentist inference; this is not fully Bayesian inference or a claim of a novel general testing principle.}
\details{
Writing W=Z-e, the score is (H-c)W(Y-m-alpha W), where c=E[H W^2]/E[W^2]. Its first derivatives in the marginal mean m, propensity e and common effect alpha vanish at their true values under the constant-CATE null. The empirical implementation estimates alpha by residual regression, c by W-squared weighting, and uses a robust covariance. It rejects singular covariance instead of adding a ridge and reporting a calibrated p-value.

Asymptotic justification requires independent observations, overlap, fixed prespecified score dimension, finite moments, nonsingular covariance and sufficiently fast nuisance convergence. In particular, a product of outcome and propensity errors and a squared propensity-error term must be o(n^-1/2), with appropriate feature-weighted control. This is not unrestricted double robustness. Fixed-score Gaussian multiplier draws do not re-estimate nuisances or capture their finite-sample uncertainty. HC3 can be substantially conservative even with true nuisances. It is a leverage adjustment of the restricted score, not unrestricted-regression HC3 Wald inference.

Constant CATE implies zero tested moments; a finite basis may miss nonconstant alternatives. No clustered-data, adaptive outcome-selected-basis, binary-link or distributional extension is claimed. BART learns the marginal outcome using training-fold outcomes and the propensity using a probit model. Propensity predictions average the normal CDF of held-out latent draws, not the CDF of their posterior mean.
}
\value{A bcforthogonalscore list containing p_value_quad (chi-square), p_value_quad_multiplier and p_value_max (Gaussian multiplier plus-one estimates), statistics, degrees of freedom, alpha, scores, covariance, projected basis, nuisance predictions, fold provenance and overlap/leverage/rank diagnostics. Simulated p-values are NA when num_null_draws is zero.}
\references{Chernozhukov et al. (2018), Double/debiased machine learning for treatment and structural parameters, The Econometrics Journal 21:C1--C68, doi:10.1111/ectj.12097.}
\examples{
set.seed(1)
X <- matrix(rnorm(300), 100, 3)
e <- plogis(X[, 1])
z <- rbinom(100, 1, e)
m <- X[, 1] + X[, 2]
y <- m + 0.5 * (z - e) + rnorm(100)
fit <- bcf_orthogonal_score_test(X, z, y,
  nuisance_method = "supplied", nuisance_source = "oracle",
  m_hat = m, e_hat = e, num_null_draws = 99)
fit$p_value_quad
}
