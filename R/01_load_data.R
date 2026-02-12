here::i_am("R/01_load_data.R")

# Import packages
library(Seurat)
library(patchwork)
library(tidyverse)

# This file will download the 10X PBMC dataset, if it's not already downloaded.
# It will then create a Seurat object from the dataset.

# Download the 10X PBMC dataset if not already downloaded
if (
  !file.exists(here::here("input", "pbmc3k_filtered_gene_bc_matrices.tar.gz"))
) {
  download.file(
    url = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
    destfile = here::here("input", "pbmc3k_filtered_gene_bc_matrices.tar.gz")
  )
}

# Extract the dataset
untar(
  here::here("input", "pbmc3k_filtered_gene_bc_matrices.tar.gz"),
  exdir = here::here("input", "pbmc3k_filtered_gene_bc_matrices")
)

# Create a Seurat object from the dataset.
# Include genes detected in at least 3 cells.
# Include cells featuring at least 200 genes.
pbmc <- "input/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19" |>
  here::here() |>
  Read10X(
    data.dir = _,
    gene.column = 2,
    cell.column = 1
  ) |>
  CreateSeuratObject(
    counts = _,
    project = "pbmc3k",
    min.cells = 3,
    min.features = 200
  )

# Print the Seurat object
show(pbmc)

# Show the structure of the Seurat object
str(pbmc)

# Show the fitst 20x20 section of the dataset
message("Printing first 20x20 section of the dataset")
pbmc |>
  GetAssayData(assay = "RNA", layer = "counts") |>
  _[1:20, 1:20] |>
  print()

# Save data
saveRDS(pbmc, file = here::here("output", "pbmc.rds"))
