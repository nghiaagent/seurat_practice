# Analysis of the 10X PBMC dataset with Seurat


This repository contains the code for a basic clustering and annotation
analysis of the 10X PBMC scRNA-seq dataset with Seurat.

# Prerequisites:

You need these software to run the code:

- R 4.5.1
- Rtools45

# Quick start

- Clone the repository.
- Run the following code in R:

``` r
renv::restore()
renv::install("presto")
# or renv::install("immunogenomics/presto") if the above fails

# Run project pipeline
source("R/99_pipe.R")
```
