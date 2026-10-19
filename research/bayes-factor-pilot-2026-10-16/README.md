# Bayes-factor feasibility pilot, 16 October 2026

For the GitHub checkpoint and commands to continue on another computer, see [LAPTOP.md](LAPTOP.md). The reviewed research package and saved synthetic results are checkpointed in the private `bcf-linear` repository; the original D: paths below record this machine's run locations.

This is an opt-in research sidecar for the local `stochtree` checkout. No package function, export, compiled code, or existing test was changed. Source the new file or run the scripts explicitly. Results live in `D:/heterogeneity-bayes-research-2026-10-16/pilot-results`; the source stays in this checkout.

The pilot asks whether the existing restricted forest fit can be reused for evidence about treatment-effect modifiers. It implements an analytically integrated Gaussian interaction alternative and a finite sparse subset mixture, plus deliberately separate plug-in and score-experiment diagnostics. **It does not implement a Bayes factor for the full horseshoe model in `R/bcf_linear_probit.R` or `R/bcf_linear.R`.**

## Connection to the current model

The draft's restricted working model is `Y = f(X,e) + alpha (Z-e) + error`. In `bcf_linear_probit`, the treatment-effect surface is `tau(X)=alpha+X beta+X_interactions beta_int`: **both** `beta` (linear modifiers) and `beta_int` (modifier-by-modifier products) describe heterogeneity. The `hc-hs` global shrinkage applies to `beta_int`; testing only this block would test nonlinear interaction structure, not a constant CATE. A constant-CATE model must set every nonconstant treatment-effect coefficient to zero. `unlink` gives the interaction block independent local scales. This sidecar retains the prognostic forest and constant treatment effect, and adds `D delta`, where `D = diag(Z-e) H` and the columns of `H` are prespecified effect-modifier functions. It replaces the local/global shrinkage hierarchy with tractable proper Gaussian slabs; this is a controlled starting model, not a claim that the two priors are equivalent. Its quadratic stress basis also differs from the package's usual pairwise-product expansion.

The current `R/bcf_score_test.R` and `R/bcf_orthogonal_score_test.R` are sourced into an isolated environment whose parent is the installed `stochtree` namespace. The installed native library is reused, without rebuilding it. Current-source wrappers and installed binary are separately fingerprinted; this does not establish that the installed binary was compiled from all present-day C++ sources.

The centered treatment `W=Z-e` follows the restricted score model. The full sparse routines have other treatment-coding and propensity-updating options. Rewriting an uncentered model `mu(X)+Z tau(X)` as `m(X)+(Z-e)tau(X)` requires `m=mu+e tau`; retaining independent priors for the newly named nuisance and treatment functions generally changes the joint prior. Equal fitted mean functions therefore do not by themselves imply equal Bayes factors. A full-model implementation must lock down coding and compatible priors explicitly.

For the actual forest experiment:

* `standardize=FALSE`: no outcome-dependent shift or scale is used in the prior.
* Shared constant-effect prior is `alpha ~ N(0,4)`.
* Noise variance is known and fixed at 1; its update and numerical floor are inactive.
* There are 25 prognostic trees, each with fixed Gaussian leaf variance 0.04; leaf-variance updating is disabled. Other tree restrictions are the current wrapper defaults, including branching parameters 0.95 and 2, maximum depth 10, and minimum leaf size 5. Observed-design restrictions are interpreted as a conditional-on-design prior.
* The propensity is supplied in the primary experiments and enters the forest as a predictor.
* Two chains run 5 GFR initialization iterations, 300 burn-in iterations, and 750 retained iterations each. Short pilot chains do not prove adequate posterior exploration.

These choices intentionally differ from the manuscript's standardization, `alpha ~ N(0,100)` on its standardized scale, sampled leaf/noise variances, and longer chains. They remove several barriers to defining a coherent evidence calculation. They do not validate the actual source sampler or establish that the finite native chain has exactly sampled the desired posterior.

## The evidence identity

Write `W = Z-e`, `r_f = Y-f`, and

\[
M_0:\quad r_f = W\alpha+\epsilon,\quad
\alpha\sim N(0,v_\alpha),\quad\epsilon\sim N(0,\sigma^2 I).
\]

Under a nonempty subset `S` of modifier columns,

\[
M_S:\quad r_f=W\alpha+D_S\delta_S+\epsilon,\qquad
\delta_S\sim N(0,s^2 I).
\]

The prior for the forest and the prior for `alpha` are identical under both models; all priors in this reduced experiment are proper. Integrating alpha gives `C0 = sigma² I + v_alpha W W'`. Let `A_S = I+s² D_S' C0^{-1}D_S` and `u_S=D_S'C0^{-1}r_f`. Then

\[
\log R_S(f)=-\tfrac12\log|A_S|
 +\tfrac12s^2u_S^\top A_S^{-1}u_S.
\]

The sidecar uses the Woodbury identity and Cholesky solves; it does not construct or invert an `n × n` covariance matrix. The full-model evidence identity for this *reduced* model is

\[
BF_{S0}=\frac{\int m_S(Y\mid f)\pi(f)df}
 {\int m_0(Y\mid f)\pi(f)df}
 = E_{f\mid Y,M_0}[R_S(f)].
\]

Thus an arithmetic mean of the ratios, implemented using log-sum-exp, estimates the BF. Averaging log ratios gives a different quantity. Integrating `alpha` before forming the ratio improves on conditioning separately on every sampled alpha; the stored forest draws still have the null marginal posterior needed by the identity.

The dense alternative includes every column, `delta ~ N(0,s² I)`. The sparse alternative puts independent inclusion probability `rho=1/q` on each column, **conditional on at least one being included**, and averages `R_S` over all `2^q-1` nonempty subsets with normalized prior weights. The null model mass is external to this alternative. The pilot compares slab standard deviations 0.25, 0.5, and 1 on the original outcome scale after covariate-only basis standardization. It does not choose the best slab after observing outcomes and report that maximum as a BF.

Equal slab SD does not match total prior signal energy. The dense prior has `E||delta||²=q s²`; the sparse prior has `E||delta||²=s² q rho/{1-(1-rho)^q}`. With six modifiers and SD 0.5 these are 1.5 and approximately 0.37588. Thus comparisons change both sparsity structure and total prior scale. The existing dense-SD-0.25 sensitivity (energy 0.375) almost matches sparse-SD-0.5 energy, but a future controlled comparison should match these quantities exactly and still check alternative scales.

Enumeration is affordable for the six/seven modifiers here; it is exponential and not a scalable sparse inference algorithm. The function rejects dimensions above 12. Extending the existing horseshoe model requires integration over proper local/global scales, or another defensible marginal-likelihood/model-switching calculation, with compatible nuisance priors and careful treatment of noise-scale coupling.

## Four quantities that must not be conflated

1. **Oracle conditional Gaussian BF:** the nuisance function is truly known in the simulation. The displayed covariance formula is exact for that conditional Gaussian model. The constant effect is integrated under its proper shared prior. Sampling-frequency statements at a fixed constant effect are not prior-predictive calibrations for every member of a composite null.
2. **Null-posterior integrated forest BF estimate:** the identity above is exact for the stipulated compatible models, while its finite-chain estimate inherits posterior approximation and potentially severe importance-weight problems. This is the actual forest-linked prototype.
3. **Plug-in conditional evidence:** replacing `f` with its posterior mean, or using cross-fitted BART nuisance predictions, changes the calculation. It is not the forest-integrated BF. Cross-fitting does not by itself turn a residual Gaussian likelihood into a genuine marginal likelihood.
4. **Approximate score BF:** for a stipulated standardized score experiment `S ~ N(theta,I_q)`, with `theta ~ N(0,gI_q)`, the normal-mixture ratio is `log BF = -q/2 log(1+g) + g Q/[2(1+g)]`. For the empirical sandwich score it is an asymptotic evidence diagnostic, not a forest BF. Fixed `g=1,10` and `g=n` impose different alternative scales and have different null evidence behavior.

## Simulation and comparisons

Each case has 300 independent observations and six independent standard-normal covariates. The marginal outcome nuisance is `m(X)=sin(X1)+0.5(X2²-1)+0.5 X3 X4`, and

`Y = m(X)+(Z-e(X)) {0.5+h(X)}+N(0,1)`.

The scenarios are constant effect (`h=0`), a sparse linear modifier (`h=.5 X1`), a dense linear modifier (`h=.5 sum(X)/sqrt(6)`), and a quadratic modifier (`h=.6 (X1²-1)/sqrt(2)`). Propensity is either 0.5 or `0.15+0.7 logistic(.6X1-.6X2+.4X3²)`. The primary basis contains the six linear covariates. A prespecified augmented basis adds the quadratic term in the quadratic scenario. This is a basis-sensitivity experiment; it is not outcome-adaptive discovery. Linear score directions have a population blind spot for the symmetric quadratic alternative under randomization. The confounded design changes the weighting and need not have exactly the same orthogonality.

There are 300 replications per oracle cell, 40 per forest cell, and 25 per cross-fit sensitivity cell (confounded null/sparse cases with known versus BART-estimated propensity). Replications share seeds across oracle and forest analyses. The cross-fit API uses three folds, 25 trees, 75 burn-in iterations, and 150 retained draws for nuisance estimation. This part uses its usual nuisance prior calibration and is explicitly only a plug-in/score sensitivity analysis.

The primary existing-model comparator is the manuscript-style pooled restricted-fit residual moment endpoint: center `W(Y-fbar-alphabar W)` and combine it with centered modifiers and its un-ridged empirical covariance. The source wrapper's per-chain internal p-values are also recorded for audit; they are not averaged into a pooled p-value. The restricted endpoint under confounding is an exploratory comparison, not a newly proved validity result. The oracle and cross-fitted comparisons use the current `bcf_orthogonal_score_test` implementation.

The summaries report frequencies of `p<.05`, `BF10>3`, `BF10>10`, `BF01>3`, `BF01>10`, and inconclusive evidence, with 95% Wilson binomial intervals. **These thresholds have different meanings and are not calibrated to equal Type-I error; comparing their raw frequencies is not a fair power ranking.** The simulation is too small for definitive calibration claims, especially in the 25/40-replication cells. No heteroskedastic/noise-misspecification, binary-outcome, unknown-noise, heavy-tail, high-dimensional or real-data validation is claimed.

## Diagnostics and reproducibility

`importance_ess = 1/sum(normalized_ratio_weight²)` is a concentration diagnostic, **not MCMC effective sample size**. The largest weight, top 1% weight share, chain-specific log-BFs, and batch-means relative SE reveal instability. Weight ESS ignores serial dependence and cannot detect a posterior region that every chain missed. Batch estimates have the same limitation. The code flags weight ESS below 100, a weight above 10%, or chain evidence differing by a factor above 3. These are diagnostic heuristics, not a certification threshold. Basic split-Rhat for alpha, average forest prediction, and forest projection onto W is also recorded; it is not the modern rank-normalized Rhat and does not guarantee a stable BF tail estimate.

`runtime` on integrated-BF rows measures the shared computation of both prior choices **and** their plug-in comparisons for that slab; it is repeated across the two prior rows and must not be summed. `fit_seconds` in the fit-diagnostics file measures fitting separately. Oracle BF runtime similarly covers both priors. Timings are descriptive on this machine.

Run from the repository root in PowerShell:

```powershell
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla research/bayes-factor-pilot-2026-10-16/test_sidecar.R
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla research/bayes-factor-pilot-2026-10-16/run_pilot.R 'D:/heterogeneity-bayes-research-2026-10-16/pilot-results' smoke
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla research/bayes-factor-pilot-2026-10-16/run_pilot.R 'D:/heterogeneity-bayes-research-2026-10-16/pilot-results' full
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla research/bayes-factor-pilot-2026-10-16/run_diagnostics.R 'D:/heterogeneity-bayes-research-2026-10-16/pilot-results'
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla research/bayes-factor-pilot-2026-10-16/analyze_results.R 'D:/heterogeneity-bayes-research-2026-10-16/pilot-results'
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla research/bayes-factor-pilot-2026-10-16/publish_pilot_report.R 'D:/heterogeneity-bayes-research-2026-10-16/pilot-results'
& ./research/bayes-factor-pilot-2026-10-16/finalize_manifest.ps1
```

The mathematical checks compare against direct multivariate-normal densities, a one-coordinate formula, all-zero prior/design limits, row/column/sign invariance, explicitly summed subset marginal densities, the zero-g score identity, and independent numerical quadrature of the null-posterior expectation identity in a fully integrated Gaussian nuisance special case. They also reject invalid inputs. The smoke run checks actual native fitting and every output path.

Full output includes raw per-replication results, event summaries with intervals, warnings/errors, fit diagnostics, two-chain likelihood-ratio traces and fitted objects for representative replications, session information, configuration and run logs. Source/native hashes are recorded separately. The R installation emitted four startup locale warnings (`LC_COLLATE`, `LC_CTYPE`, `LC_MONETARY`, `LC_TIME`); the run uses locale `C`. These occurred before the R logging sink and are recorded here explicitly.

The optional long-chain diagnostic selects the dataset with the largest dense-prior chain log-BF disagreement in each null/sparse-by-design cell. It reruns these four datasets with four chains under a fresh call seed, 1,000 burn-in and 3,000 retained draws each. The wrapper advances the call's R/native random streams across chains; it does not assign a separate user seed to every chain. This is a deliberately diagnostic-selected stress check, not an independent set of operating-characteristic replications. Its case selection, chain values, and ratio/parameter traces are saved separately.

## Most useful next steps

Use the actual finite-chain results to decide whether null-posterior averaging is numerically viable near the null; do not assume it remains viable under strong alternatives. A trustworthy continuation should implement bridge sampling or an augmented model indicator when the null posterior and alternative posterior barely overlap, use longer independently initialized chains, audit the native target, integrate unknown noise and shrinkage scales under explicitly fixed proper priors, and compare methods at matched null error rates using much larger simulations. The nonlinear basis check should remain part of that protocol: evidence for zero coefficients in a selected basis is not proof that the CATE is globally constant.
