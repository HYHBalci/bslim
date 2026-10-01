# Portable heterogeneity-research entry points

The analysis and paper live in the sibling bcf-linear/heterogeneity-project checkout. This folder provides a dependency checker, opt-in source installer, and a fresh synthetic diagnostic smoke. It changes no sampler.

From the stochtree checkout:

    Rscript research/heterogeneity/bootstrap.R --check

The checker reads Imports and LinkingTo from this source's DESCRIPTION, checks research dependencies, checks required vendored headers, and verifies that the installed package exports bcf_score_diagnostics. It exits 2 if the environment is incomplete. Reports are written to local-output.

To install missing R dependencies and compile this checkout:

    Rscript research/heterogeneity/bootstrap.R --install

The isolated default library is .research-library/R-<version>-<platform>. Set STOCHTREE_R_LIBRARY to choose another library. Set R_REPOSITORY to choose an R package repository; the default is https://cloud.r-project.org. The installer uses versions available there at installation time and records actual versions. It is not a historical dependency lock.

A source build needs a C++17 compiler. On Windows use the Rtools release corresponding to the installed R version; on macOS a configured compiler toolchain; on Linux an R development toolchain. Populate the repository submodules when obtaining a clean checkout, for example using git clone --recurse-submodules. The R build uses Eigen, fmt, and fast_double_parser headers under deps; the checker reports missing headers. A fresh source build was validated on 1 October 2026 in .research-library/fresh-2026-10-01 using R 4.4.3 and the available C++17 toolchain. All 226 targeted diagnostic/AIPW/orthogonal/probit assertions, 29 live horseshoe checks, and 32 extracted-source horseshoe checks passed. Five installed API bodies and formals matched source; the installed-namespace 250-tree example passed. The historical private DLL was not changed. Evidence is in the sibling project update fresh_build_validation.json and fresh_r_packages.csv.

Then run:

    Rscript research/heterogeneity/diagnostics_demo.R

This creates one deterministic synthetic dataset with known propensity 1/2, fits a short two-chain restricted forest, reconstructs its score with the public diagnostic helper, and checks the reconstruction and covariance rank. Its short MCMC budget is a functionality check and is not evidence of convergence, calibration, or power.

For an existing compatible binary plus current R source, explicitly use:

    Rscript research/heterogeneity/diagnostics_demo.R --source-overlay

This mode loads the current sampler and diagnostic R functions into a private environment over the chosen installed namespace. It labels that runtime mode and does not mutate the installed package. The retained September 30 binary was used to validate this smoke, yielding reconstruction error below 1e-16. The default installed-namespace mode correctly refuses an older package without the exported helper.

The historical study used an earlier private Windows stochtree DLL plus frozen R source. A fresh compile, even with version 0.1.1, is a distinct runtime. See the [portability guide](../../../bcf-linear/heterogeneity-project/research-update-2026-10-01/PORTABILITY.md) in a sibling layout for the recorded binary hash, observed dependency versions, and external artifact requirements. The simpler base-R companion and Python CSV verifier there require no stochtree DLL.

Do not commit .research-library, local-output, compiled binaries, or full simulation records. Source and seed identities should be recorded before a new study; these examples do not overwrite or append to old experiments.
