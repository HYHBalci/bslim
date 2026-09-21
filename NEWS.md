# stochtree 0.1.1.9004

## Breaking Changes

* `bcf_linear_probit(..., general_params = list(orthogonalize = TRUE))`
  now uses fixed, externally supplied propensity scores to code treatment as
  `Z - propensity`. It no longer fits or updates an internal propensity model;
  internal propensity coefficients (`theta`) are no longer sampled.
* This mode requires raw binary treatment values (`0` or `1`) and one finite
  propensity score in `[0, 1]` per observation. `propensity_train` is required,
  as is `propensity_test` when test covariates are supplied. Prediction requires
  a corresponding `propensity` vector.
* Orthogonalisation disables adaptive treatment coding and excludes propensity
  scores from the prognostic forest (`propensity_covariate = "none"`).
  `propensity_seperate` must be `"none"`. The treatment-effect contrast is
  unchanged; Gaussian outcome predictions use `mu + (Z - propensity) * tau`.
* Orthogonalisation requires `standardize_cov = FALSE` and `keep_burnin = FALSE`.
  If covariate scaling is needed, transform training and evaluation covariates
  externally using moments estimated on the training data. Retained GFR draws,
  including GFR-only fits, are supported; retained burn-in coefficient histories
  are not.
* Fixed external-propensity orthogonalisation supports the corrected Gaussian
  hc-hs residual-variance update. The former guard for the internally sampled
  propensity branch has been removed. Non-orthogonalised models retain their
  existing behaviour.

# stochtree 0.1.2

## New Features

* Support for binary outcomes in BART and BCF with a probit link ([#164](https://github.com/StochasticTree/stochtree/pull/164))

## Bug Fixes

* Fixed indexing bug in cleanup of grow-from-root (GFR) samples in BART and BCF models
* Avoid using covariate preprocessor in `computeForestLeafIndices` function when a `ForestSamples` object is provided (rather than a `bartmodel` or `bcfmodel` object)

# stochtree 0.1.1

* Fixed initialization bug in several R package code examples for random effects models

# stochtree 0.1.0

* Initial release on CRAN.
* Support for sampling stochastic tree ensembles using two algorithms: MCMC and Grow-From-Root (GFR)
* High-level model types supported:
    * Supervised learning with constant leaves or user-specified leaf regression models
    * Causal effect estimation with binary or continuous treatments
* Additional high-level modeling features:
    * Forest-based variance function estimation (heteroskedasticity)
    * Additive (univariate or multivariate) group random effects
    * Multi-chain sampling and support for parallelism
    * "Warm-start" initialization of MCMC forest samplers via the Grow-From-Root (GFR) algorithm
    * Automated preprocessing / handling of categorical variables
* Low-level interface:
    * Ability to combine a forest sampler with other (additive) model terms, without using C++
    * Combine and sample an arbitrary number of forests or random effects terms
