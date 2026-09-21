# Fixed external-propensity orthogonalisation (0.1.1.9004)

`bcf_linear_probit()` now interprets `general_params$orthogonalize = TRUE`
as fixed external treatment coding. Supply raw binary treatment and an external
propensity per row; the sampler uses `Z - propensity`, disables adaptive coding,
and excludes propensity from the prognostic forest. No propensity model is
fitted internally. Predictions use the same coding, and `tau_hat` is the
treatment-effect contrast, not its product with the residualised treatment.

Use the usual supported Gaussian hc-hs settings (`sample_global_prior = "hc-hs"`,
`unlink = TRUE`, `gibbs = TRUE`, `hn_scale = 1`). This changes neither the
coefficient-update kernel nor the corrected residual-variance formula. The
fixed basis simply changes the design and corresponding residuals. The other
unsupported variance-update variants remain guarded.

`propensity_train` is required. With `X_test`, supply binary `Z_test` and
`propensity_test`; for `predict()` or `predict_linear_bcf()`, supply binary `Z`
and `propensity`. Scores must be finite, in [0, 1], with exactly one per row.
`propensity_seperate` must be `"none"`. Keep `standardize_cov = FALSE` and
`keep_burnin = FALSE`; pre-transform covariates using training-data moments if
needed. Outcome standardization remains supported.

## Verification performed on 21 September 2026

Run against an installed package from the repository root:

```r
# Shell: Rscript tools/check-fixed-external-orthogonalization.R [R_LIBRARY]
```

- 163 new fixed-orthogonalisation assertions and 12 existing prediction
  assertions passed against an isolated installed package (175 total).
- The new mode matched ordinary fitting with manually residualised treatment
  to 1e-12 for coefficient, shrinkage, residual-variance and prognostic draws,
  both with and without GFR. Training and held-out predictions agreed.
- Tests covered invalid inputs, no internal propensity model, forest exclusion,
  caller-input preservation, pure GFR, retained GFR, multiple thinned chains,
  coefficient/forest draw alignment, outcome scaling, endpoint probabilities,
  seed reproducibility, and both prediction entry points.
- Updated R and Rd files parsed, and `git diff --check` passed.

The local isolated installation used the updated R package and the unchanged,
previously compiled repaired native DLL from version 0.1.1.9003. Native source
was not modified. This was not a fresh C++ rebuild or a full R CMD check.
The user's default R library was not replaced. Tiny test fits are not evidence
of convergence, frequentist coverage, or completed scientific reruns.

The broader existing `test-bcf.R` suite had 13 passing assertions and one error:
its warm-start test calls the absent `saveBCFModelToJsonString()` function.
The identical error was reproduced on the unchanged 0.1.1.9003 installation;
it is not introduced by this change and is not included in the scoped checker.

## Migration

This intentionally replaces the old internally sampled propensity branch.
There are no `theta_samples`. Reinstall this version before running studies;
scripts expecting internal propensity draws must be updated. Do not mix old
internal-propensity fits with new fixed-external fits in one analysis or resume
directory. Archived fitted objects are not relabelled as fixed-external models.
