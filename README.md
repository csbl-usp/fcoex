# fcoex

A package for FCBF-based coexpression analysis of single cell data.

The fcoex package implements an easy-to use interface to co-expression analysis based on the FCBF (Fast Correlation-Based Filter) algorithm. The package structure is adapted from the CEMiTool package, relying on visualizations and code designed and written by CEMiTool's authors.

## Installing the package

fcoex is now on Bioconductor (https://bioconductor.org/packages/fcoex/)

To install it, just run:`

```

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel')

BiocManager::install("fcoex")

```

It is also possible to install from source code using devtools

```

install.packages(devtools)
devtools::install_github("csbl-usp/fcoex")

```




