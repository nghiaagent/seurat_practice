here::i_am("R/02_preprocess.R")

# Import packages
library(Seurat)
library(GGally)
library(ggforce)
library(glmGamPoi)
library(patchwork)
library(tidyverse)

# In this file, we perform Seurat pre-processing
# QC to remove low-quality cells (empty, dying, multiplets) using gene counts
# QC to remove low-quality cell with mitochondrial contamination
# Visualise QC metrics
# Normalise using SCTransform
# Calculate PCA, UMAP, clusters
# Save data

# Load data
pbmc <- readRDS(here::here("output", "pbmc.rds"))

# QC
## Add percentage of mitochondrial genes to the object metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## Select filtering metrics.
## We filter out cells with >2500 or <200 unique genes
## and >5% mitochondrial contamination
### Set threshold
qc_nfeature_low <- 200
qc_nfeature_high <- 2500
qc_percentmt <- 5

### Filter cells
pbmc_filtered <- pbmc |>
  subset(
    subset = nFeature_RNA > qc_nfeature_low &
      nFeature_RNA < qc_nfeature_high &
      percent.mt < qc_percentmt
  )

## Visualise QC and filtermetrics
### Create vector for labelling metrics
plot_qc_labeller <- as_labeller(
  c(
    "nFeature_RNA" = "Number of genes per cell",
    "nCount_RNA" = "Number of counts per cell",
    "percent.mt" = "Percentage of mitochondrial genes"
  )
)

### Create faceted violin plot
plot_qc_violins <- pbmc@meta.data |>
  # Mark cells that are kept
  mutate(
    keep = nFeature_RNA > qc_nfeature_low &
      nFeature_RNA < qc_nfeature_high &
      percent.mt < qc_percentmt
  ) |>
  # Make table longer for visualisation
  pivot_longer(
    cols = c(nCount_RNA, nFeature_RNA, percent.mt),
    names_to = "metric",
    values_to = "value"
  ) |>
  # Maintain column order and label columns
  mutate(
    metric = metric |>
      factor(levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  ) |>
  # Perform plotting
  ggplot(aes(x = metric, y = value)) +
  geom_violin() +
  geom_jitter(aes(colour = keep), alpha = 0.1, size = 0.8) +
  scale_colour_manual(
    values = c(
      "TRUE" = "black",
      "FALSE" = "red"
    )
  ) +
  theme_bw() +
  facet_wrap(~metric, scales = "free", labeller = plot_qc_labeller)

### Create scattermatrix of QC metrics
plot_qc_scattermatrix <- pbmc@meta.data |>
  # Mark cells that are kept
  mutate(
    keep = nFeature_RNA > qc_nfeature_low &
      nFeature_RNA < qc_nfeature_high &
      percent.mt < qc_percentmt
  ) |>
  ggpairs(
    columns = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    lower = list(
      mapping = aes(colour = keep),
      continuous = wrap(
        ggally_autopoint,
        alpha = 0.1,
        size = 0.8
      )
    ),
    labeller = plot_qc_labeller
  ) +
  scale_colour_manual(
    values = c(
      "TRUE" = "black",
      "FALSE" = "red"
    )
  ) +
  theme_bw()

### Combine and save plots
wrap_plots(
  plot_qc_violins,
  ggmatrix_gtable(plot_qc_scattermatrix),
  ncol = 1,
  heights = c(1, 2)
) |>
  ggsave(
    filename = here::here("output", "plot_qc.png"),
    width = 8,
    height = 12
  )

# Normalisation using SCTransform
pbmc_norm <- pbmc_filtered |>
  SCTransform(vars.to.regress = "percent.mt")

# Run PCA and UMAP
pbmc_norm <- pbmc_norm |>
  ## Run PCA, find neightbors using computed PCs
  RunPCA() |>
  FindNeighbors(dims = 1:30) |>
  ## Run UMAP using computed PCs, find clusters
  RunUMAP(dims = 1:30) |>
  FindClusters(resolution = 0.5)

# Save data
saveRDS(pbmc_norm, file = here::here("output", "pbmc_norm.rds"))
