# Continue the Bayes-factor research on another computer

Checkpoint prepared on 19 October 2026. The code repository is now named `HYHBalci/spbart` on GitHub; the former `HYHBalci/stochtree` URL redirects there. The R package remains named `stochtree`.

## Get both repositories

Clone these branches into sibling folders, using the GitHub account with access to the private research repository:

```powershell
git clone --recurse-submodules --branch codex/heterogeneity-evidence-2026-10-01 https://github.com/HYHBalci/spbart.git stochtree
git clone --branch codex/heterogeneity-competitor-repair-2026-10-12 https://github.com/HYHBalci/bcf-linear.git bcf-linear
```

For existing clones, fetch origin and switch to the named branch, then pull with `--ff-only`; preserve any existing local changes before switching. Initialize the code submodules with `git submodule update --init --recursive`.

The reviews, proposal, four PDFs, source records, Grok scouting and saved synthetic pilot results are in:

```text
bcf-linear/heterogeneity-project/bayes-evidence-2026-10-16/
```

Start with that folder's `CONTINUE_ON_LAPTOP.md` and `DISCUSSION_2026-10-19.md`. Its `COPIED_ARTIFACT_SHA256.csv` records the byte-identical copies of the original D: research artifacts. Per-page rendering previews are reproducible and are not part of the checkpoint.

## Check the mathematics

From the `stochtree` checkout, using an Rscript installation on PATH:

```powershell
Rscript --vanilla research/bayes-factor-pilot-2026-10-16/test_sidecar.R
```

These checks use base R and do not require an installed forest package. They passed on the checkpoint machine using R 4.4.3. R emitted the same four startup locale warnings documented in the original run.

## Set up a fresh forest runtime

A source build requires R and a compatible C++17 toolchain (Rtools on Windows), along with the initialized source submodules. The existing bootstrap installs missing dependencies and compiles this checkout into an isolated library:

```powershell
$env:STOCHTREE_R_LIBRARY = Join-Path (Get-Location) '.research-library/laptop'
Rscript research/heterogeneity/bootstrap.R --install
$env:R_LIBS_USER = $env:STOCHTREE_R_LIBRARY
Rscript research/heterogeneity/bootstrap.R --check
Rscript --vanilla research/bayes-factor-pilot-2026-10-16/run_pilot.R research/bayes-factor-pilot-2026-10-16/local-output/laptop-smoke smoke
```

Use a new output directory for new experiments. The pilot takes the output directory as its first argument, so a D: drive is unnecessary. For a full rerun, replace `smoke` with `full` and choose a separate output directory.

The bootstrap uses dependencies available at installation time; it is not a historical version lock. The saved `session_info.txt`, `provenance.json` and source/binary hashes describe the historical experiment. A fresh native build is a new runtime and should be recorded as such. Compiled libraries are not transferred through Git.

## Research state

The pilot implements proper Gaussian slabs and sparse subset averaging, with analytic integration of the common effect. It does not implement a full horseshoe Bayes factor. Forest-integrated evidence was unstable in selected diagnostic reruns; establish an independent evidence estimator before treating the numbers as reliable evidence or launching a large study.

Both `beta` and `beta_int` describe treatment-effect heterogeneity. A global constant-effect test must remove both blocks. A fixed modifier basis yields a global test within that basis. A broader omnibus claim requires a rich treatment-effect alternative and additional theory; flexibility of the prognostic forest alone does not establish it.
