here::i_am("R/04_cluster_annotation.R")

# Import packages
library(presto)
library(Seurat)
library(tidyverse)

# Load data
pbmc_norm <- readRDS(here::here("output", "pbmc_norm.rds"))

# In this file, we perform cluster annotation
# Find marker genes for each cluster
# Annotate clusters
# Save data

# Find marker genes for each cluster
pbmc_markers <- pbmc_norm |>
  FindAllMarkers()

pbmc_markers_filtered <- pbmc_markers |>
  group_by(cluster) |>
  filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) |>
  slice_min(n = 20, order_by = p_val_adj)

# Plot heatmap of top 20 markers for each cluster
top_markers <- pbmc_markers |>
  group_by(cluster) |>
  filter(avg_log2FC > 1 & p_val_adj < 0.05) |>
  slice_min(n = 20, order_by = p_val_adj)

plot_heat_top_markers <- DoHeatmap(pbmc_norm, features = top_markers$gene)

# Manually annotate clusters
## Create cluster ID mapping, derived from PanglaoDB search
cluster_ids <- c(
  "Naive CD4",
  "CD14+ Mono",
  "Memory CD4 T",
  "B cells",
  "CD8+ T",
  "FCGR3A+ Mono",
  "NK",
  "DC",
  "Macrophages",
  "Platelet"
) |>
  set_names(levels(pbmc_norm))

## Name clusters in dataset
pbmc_annotated <- pbmc_norm |>
  RenameIdents(cluster_ids)

## Visualise annotated clusters
plot_umap_annotated <- pbmc_annotated |>
  DimPlot(reduction = "umap", label = TRUE) +
  NoLegend()

# Save data
saveRDS(
  pbmc_markers,
  file = here::here("output", "pbmc_markers.rds")
)

saveRDS(
  pbmc_annotated,
  file = here::here("output", "pbmc_annotated.rds")
)

ggsave(
  plot = plot_umap_annotated,
  filename = here::here("output", "plot_umap_annotated.png"),
  width = 10,
  height = 10
)
