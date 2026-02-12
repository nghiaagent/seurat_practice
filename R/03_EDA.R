here::i_am("R/03_EDA.R")

# Import packages
library(ggforce)
library(Seurat)
library(patchwork)
library(tidyverse)

# Load data
pbmc_norm <- readRDS(here::here("output", "pbmc_norm.rds"))

# In this file, we perform exploratory data analysis
# Visualise PCA and UMAP
# Save data

# Visualise PCA and UMAP
## Visualise PCA
plot_pca <- pbmc_norm |>
  DimPlot(reduction = "pca")

plot_pca_heat <- pbmc_norm |>
  DimHeatmap(dims = 1:10, cells = 500, balanced = TRUE)

## Visualise UMAP
plot_umap <- pbmc_norm |>
  DimPlot(reduction = "umap", label = TRUE)

# Save plots
saveRDS(
  list(
    "PCA" = plot_pca,
    "PCA_heatmaps" = plot_pca_heat,
    "UMAP" = plot_umap
  ),
  file = here::here("output", "plots_eda.rds")
)

# Export images
plot_output <- wrap_plots(plot_pca, plot_umap, ncol = 2) +
  plot_annotation(tag_levels = list(c("PCA", "UMAP"))) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 12, hjust = 0, vjust = 0)
  )

ggsave(
  plot = plot_output,
  filename = here::here("output", "plots_pca_umap.png"),
  width = 10,
  height = 5
)
