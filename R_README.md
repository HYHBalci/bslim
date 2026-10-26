# B-SLIM R package

BART Sparse Linear Interaction Model (B-SLIM) for causal inference,
based on [StochTree](https://github.com/StochasticTree/stochtree). The original
authors, licenses and copyright notices are retained.

## Getting started

Clone `bslim` with its source dependencies, then install it with R and a C++17
compiler (Rtools on Windows). Install the R package `remotes` first if needed.

```
git clone --recurse-submodules https://github.com/HYHBalci/bslim.git
cd bslim
Rscript -e 'remotes::install_deps(".", dependencies = TRUE, upgrade = "never")'
R CMD INSTALL .
```

The package can be installed alongside `stochtree`. Use explicit namespaces when
working with both packages in one R session:

```
library(bslim)
bslim::bcf_linear_probit
stochtree::bcf
```
