\name{bcf_aipw_score_test}
\alias{bcf_aipw_score_test}
\title{Experimental Cross-Fitted Weighted AIPW Test of Mean Treatment Heterogeneity}
\usage{
bcf_aipw_score_test(
  X_train, Z_train, y_train, score_basis = NULL,
  nuisance_method = c("bart", "glm", "supplied", "custom"),
  m0_hat = NULL, m1_hat = NULL, e_hat = NULL,
  nuisance_source = "cross_fitted", known_propensity = NULL,
  fold_id = NULL, num_folds = 5L, training_indices = NULL,
  nuisance_fitter = NULL, weighting = c("overlap", "uniform"),
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
\item{score_basis}{Fixed, prespecified finite numeric modifier matrix, defaulting to X_train. Nuisance covariates and score directions may differ.}
\item{nuisance_method}{Armwise BART, armwise linear outcome and logistic propensity GLM, supplied external predictions, or a custom training-fold callback.}
\item{m0_hat,m1_hat,e_hat}{Supplied predictions of E[Y|X,Z=0], E[Y|X,Z=1], and P(Z=1|X). Allowed only with nuisance_method="supplied".}
\item{nuisance_source}{For supplied predictions, a common string "cross_fitted" or "oracle", or a named character vector/list declaring separate sources with keys m0, m1, e. Both outcome sources must be declared; e may be omitted when known_propensity is given. Oracle explicitly asserts the true nuisance function, never a full-sample fitted prediction. Known propensity is recorded separately as source "known".}
\item{known_propensity}{Known design propensity, scalar or length n, strictly inside (0,1). Do not use this argument for an estimated full-sample propensity. Cannot be combined with e_hat.}
\item{fold_id}{One common held-out fold vector, or a named list of vectors for supplied cross-fitted nuisances. Internally fitted nuisances use common folds generated independently of all observations when absent. Supplied oracle/known nuisances do not require folds.}
\item{num_folds}{Number of generated balanced folds, an integer between 2 and n. Training folds must contain both treatment arms; the function never repairs folds using held-out labels. Armwise BART requires at least two training observations in each arm.}
\item{training_indices}{For supplied cross-fitted inputs, either one common list of training row indices in unique(fold_id) order, or a named list of such lists, keyed by nuisance. Every training list must exclude its corresponding held-out fold. Arm-specific training subsets are allowed. Declared provenance is audited; actual external training honesty cannot be inferred from predictions. Internally fitted models generate these records automatically.}
\item{nuisance_fitter}{Custom function with exactly the data interface X_train, Z_train, y_train, X_test. It receives only training-fold labels and held-out covariates. Return numeric held-out vectors named m0_hat, m1_hat, and e_hat; e_hat may be omitted when known_propensity is supplied. Internally recorded provenance conservatively records all rows passed to the callback.}
\item{weighting}{"overlap" uses propensity weights e(1-e); "uniform" uses unit weights. These generally target different finite-basis moments under alternatives.}
\item{variance_correction}{HC0 (default) centres unadjusted score contributions for the empirical sandwich. HC3 divides contributions by one minus leverage of sqrt(omega)[1,H-c] and recentres. This is a sensitivity heuristic, not unrestricted-regression HC3 or a validated finite-sample correction.}
\item{num_null_draws}{At least 99 Gaussian multiplier draws, or zero to omit simulated p-values. Draws keep nuisance predictions fixed.}
\item{seed}{Nonnegative integer within the R integer range. Per-learner R and C++ seeds are deterministic and recorded. The caller's R random-number state is restored on success or failure, including when none existed initially.}
\item{propensity_clip}{Optional explicit bound in (0,.5), with no clipping by default. Raw estimated propensities must be probabilities in [0,1]. Clipping is reported, warns if used, and may introduce bias. True propensities must lie inside the clipping interval for usual rate arguments.}
\item{rank_action}{Error on constant/redundant score columns or explicitly drop them. Rank detection uses centred, scaled columns; retained score columns remain in their original units. Outcome-dependent basis selection is unsupported.}
\item{num_gfr,num_burnin,num_mcmc}{BART sampling budgets. Nonnegative integers, with at least one retained MCMC draw. They are prediction budgets, not posterior-convergence certification.}
\item{bart_params}{General BART parameters. The function controls random_seed and excludes GFR/burn-in draws from prediction means. Probit propensity variance is fixed as required by the probit model.}
\item{mean_forest_params}{BART mean-forest parameters, defaulting to 50 trees. The function sets continuous outcome models for m0/m1 and a probit model for unknown e.}
}
\description{
An experimental finite-basis test against constant conditional average treatment effect, using arm-specific outcome means and a propensity in a weighted augmented inverse-probability pseudo-outcome. Bayesian nuisance fitting is combined with frequentist score inference; this is neither fully Bayesian inference nor a claim of a new general AIPW principle.
}
\details{
Let D=m1-m0+Z(Y-m1)/e-(1-Z)(Y-m0)/(1-e), let omega be e(1-e) or 1, and set t=omega D. The implementation estimates alpha=sum(t)/sum(omega), c=sum(omega H)/sum(omega), and contributions psi=(H-c)(t-alpha omega). Overlap weighting evaluates the algebraically equivalent expression without inverse propensity factors. The HC0 covariance is the cross-product of centred contributions. The quadratic statistic uses the inverse covariance and a chi-square reference with dimension equal to the retained score basis. Singular or severely ill-conditioned covariance produces an error; no ridge-adjusted p-value is returned.

Constant CATE implies zero population moments, but the converse need not hold for a finite basis. AIPW has a double-robust conditional-mean identity when either the propensity or both arm-specific outcome means are correct. That identity alone does not guarantee valid p-values if one learned nuisance is inconsistent. Fixed-dimensional inference requires independent observations, overlap, suitable moments, nonsingular covariance, and nuisance convergence/product-error conditions appropriate to the weighting and training scheme. This implementation establishes no BART nuisance convergence rate and makes no general small-sample calibration claim. Estimated weights can change the target under alternatives.

Per-nuisance provenance allows, for example, oracle m0, cross-fitted m1, and known design e. Every cross-fitted supplied prediction must have its own valid fold and training records, which may be shared or supplied separately. A training subset need not include all non-held-out observations, because arm-specific models use only their treatment arm. Full-sample estimated predictions must never be relabelled as oracle or known. Generated folds are independent of X,Y,Z and never treatment-stratified or repaired after inspecting held-out labels. Supplied fold independence and fixed basis selection remain the caller's responsibilities. Distinct nuisance partitions additionally require a suitable common refinement (for example, a fixed number of observation-independent refinement cells) and nuisance-rate control on those cells; per-vector held-out declarations alone do not prove an independent-fold asymptotic theorem.

BART and GLM fit outcome regressions separately in the two treatment arms within each training fold. Unknown propensity uses all training rows and either a probit BART or logistic GLM. BART propensity predictions average the normal CDF of latent draws, not the CDF of their posterior mean. Custom callbacks see training outcomes/treatments and held-out covariates only. Fixed-score Gaussian multiplier draws do not refit nuisances or capture their finite-sample estimation uncertainty. No clustered-data, outcome-selected-basis, binary-link or distributional-testing extension is claimed.
}
\value{
A bcfaipwscore list with p_value_quad (chi-square), p_value_quad_multiplier and p_value_max (plus-one Gaussian multiplier probabilities), statistic_quad, statistic_max, df, alpha, projection, score_sum, covariance, scores, multiplier_scores, basis, projected_basis, weighted pseudo-outcome t, weights omega, and nuisance_predictions. nuisance_provenance contains a source, learner, fold vector and training-index list for each of m0, m1, e, plus seeds when stochastic learners are fitted internally. training_indices is a named list of those per-nuisance lists; it should not be confused with a single common training-fold list. Nuisance predictions include raw_e_hat as well as the propensity actually used. Diagnostics report overlap, clipping, weight effective sample size, leverage, rank and covariance conditioning. known_propensity describes the raw design input; known_propensity_modified explicitly flags clipping of that input, used_propensity_source becomes "clipped_known", and the propensity provenance records clipped row indices. Simulated p-values are NA when num_null_draws is zero.
}
\references{Chernozhukov et al. (2018), Double/debiased machine learning for treatment and structural parameters, The Econometrics Journal 21:C1--C68, doi:10.1111/ectj.12097.}
\examples{
set.seed(10)
X <- matrix(rnorm(360), 120, 3)
e <- rep(0.5, 120)
z <- rbinom(120, 1, e)
m0 <- sin(X[,1])
m1 <- m0 + 0.5
y <- (1-z)*m0 + z*m1 + rnorm(120)
fit <- bcf_aipw_score_test(X, z, y,
  nuisance_method="supplied", m0_hat=m0, m1_hat=m1,
  nuisance_source=c(m0="oracle", m1="oracle"),
  known_propensity=e, num_null_draws=99)
fit$p_value_quad
}