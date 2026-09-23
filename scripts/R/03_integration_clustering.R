#!/usr/bin/env Rscript
# =============================================================================
# 03_integration_clustering.R — TF-IDF normalization, LSI, Harmony, clustering
#
# Purpose : Normalize the merged scATAC-seq object, reduce dimensionality with
#           LSI (Latent Semantic Indexing = TF-IDF + SVD), correct batch effects
#           across the 4 samples with Harmony, then cluster and embed with UMAP.
#
# Key decisions:
#   - Component 1 of LSI is ALWAYS excluded — it correlates with sequencing
#     depth (not biology) and contaminates downstream clustering if kept.
#   - Harmony operates on LSI embeddings (not raw counts) to correct sample
#     batch effects without distorting the accessibility signal itself.
#   - SLM clustering algorithm (algorithm=3) is Signac's recommendation for
#     the sparse, near-binary count matrices typical of chromatin accessibility.
#   - clustree sweeps multiple resolutions to help select the stable one.
#
# Input   : atac_merged_consensus.rds
# Output  : atac_integrated_clustered.rds, UMAP plots, metadata CSV
#
# Tutorial: https://ngs101.com/scatac-seq-complete-beginners-guide-part-3
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(harmony)
  library(clustree)
  library(ggplot2)
  library(patchwork)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--merged_rds",    required = TRUE)
parser$add_argument("--lsi_dims",      default = "2:30",
                    help = "LSI dimensions to use, e.g. '2:30' (1 excluded)")
parser$add_argument("--harmony_dims",  default = "2:30")
parser$add_argument("--resolutions",   default = "0.2,0.4,0.6,0.8,1.0",
                    help = "Comma-separated clustering resolutions to sweep")
parser$add_argument("--min_cells",     type = "integer", default = 20)
parser$add_argument("--outdir",        default = ".")
parser$add_argument("--plot_dir",      default = "plots")
args <- parser$parse_args()

set.seed(1234)
dir.create(args$outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

# Parse dimension and resolution arguments
parse_dims <- function(s) {
  if (grepl(":", s)) {
    parts <- as.integer(strsplit(s, ":")[[1]])
    seq(parts[1], parts[2])
  } else {
    as.integer(strsplit(s, ",")[[1]])
  }
}
lsi_dims     <- parse_dims(args$lsi_dims)
harmony_dims <- parse_dims(args$harmony_dims)
resolutions  <- as.numeric(strsplit(args$resolutions, ",")[[1]])

cat("=== Integration & Clustering ===\n")
cat("LSI dims:", paste(range(lsi_dims), collapse = "-"), "\n")
cat("Harmony dims:", paste(range(harmony_dims), collapse = "-"), "\n")
cat("Resolutions:", paste(resolutions, collapse = ", "), "\n")

# ── 1. Load merged object ─────────────────────────────────────────────────────
atac <- readRDS(args$merged_rds)
cat("Loaded:", ncol(atac), "cells,", nrow(atac), "peaks\n")

# ── 2. TF-IDF normalization ───────────────────────────────────────────────────
# Method 1: log-TF × log-IDF — the published standard for scATAC-seq
# FindTopFeatures selects peaks detected in >= min_cells cells to reduce noise
atac <- RunTFIDF(atac, method = 1)
atac <- FindTopFeatures(atac, min.cutoff = args$min_cells)

cat("Variable features:", length(VariableFeatures(atac)), "\n")

# ── 3. Singular value decomposition (LSI) ────────────────────────────────────
# LSI = TF-IDF + SVD — the standard dimensionality reduction for scATAC-seq
atac <- RunSVD(atac, n = 50, reduction.name = "lsi", reduction.key = "LSI_")

# ── 4. Diagnose depth correlation — confirm component 1 should be excluded ────
p_depth_cor <- DepthCor(atac, n = 30)
ggsave(file.path(args$plot_dir, "integration_01_depth_correlation.png"),
       p_depth_cor, width = 8, height = 4)

# The correlation plot should show component 1 with |r| > 0.9 (depth-driven).
# Components 2+ should have low depth correlation and represent biology.
depth_cors <- cor(Embeddings(atac, "lsi"), log1p(atac$nCount_peaks))
cat("Component 1 depth correlation:", round(depth_cors[1, 1], 3),
    "(should be > 0.8 to justify exclusion)\n")

# ── 5. Harmony batch correction ───────────────────────────────────────────────
# Corrects for sample-of-origin differences within LSI space.
# project.dim=FALSE skips re-projecting peak loadings (not needed here).
cat("Running Harmony...\n")
atac <- RunHarmony(
  object       = atac,
  group.by.vars = "sample_id",
  reduction.use = "lsi",
  dims.use      = harmony_dims,         # component 1 excluded here
  reduction.save = "harmony",
  project.dim   = FALSE,
  verbose       = FALSE
)

# ── 6. UMAP embedding ─────────────────────────────────────────────────────────
atac <- RunUMAP(
  object         = atac,
  reduction      = "harmony",
  dims           = seq_along(harmony_dims),
  reduction.name = "umap"
)

# ── 7. Neighbor graph + multi-resolution clustering ──────────────────────────
atac <- FindNeighbors(
  object    = atac,
  reduction = "harmony",
  dims      = seq_along(harmony_dims),
  k.param   = 20
)

# Sweep resolutions with SLM algorithm (algorithm=3, best for sparse data)
for (res in resolutions) {
  atac <- FindClusters(
    object     = atac,
    algorithm  = 3,
    resolution = res,
    verbose    = FALSE
  )
  cat(sprintf("  Resolution %.1f: %d clusters\n", res, length(unique(atac@active.ident))))
}

# ── 8. Clustree — choose the most stable resolution ──────────────────────────
p_clustree <- clustree(atac, prefix = "peaks_snn_res.")
ggsave(file.path(args$plot_dir, "integration_02_clustree.png"),
       p_clustree, width = 10, height = 10)

# Set active identity to the middle resolution (0.6) as a reasonable default;
# inspect the clustree plot and adjust params.cluster_resolution if needed.
default_res_col <- paste0("peaks_snn_res.", resolutions[ceiling(length(resolutions)/2)])
Idents(atac)    <- atac[[default_res_col, drop = TRUE]]
atac$cluster_final <- Idents(atac)

cat("Active resolution column:", default_res_col, "\n")
cat("Clusters:", nlevels(Idents(atac)), "\n")

# ── 9. Diagnostic UMAP plots ─────────────────────────────────────────────────
p_clusters <- DimPlot(atac, group.by = "cluster_final", label = TRUE, repel = TRUE) +
  ggtitle("Clusters") + NoLegend()

p_samples <- DimPlot(atac, group.by = "sample_id") +
  ggtitle("Sample")

p_condition <- DimPlot(atac, group.by = "condition") +
  ggtitle("Condition")

# QC metrics projected onto UMAP (technical artifacts appear as separate "blobs")
p_tss <- FeaturePlot(atac, features = "TSS.enrichment", max.cutoff = "q95") +
  ggtitle("TSS enrichment")
p_frip <- FeaturePlot(atac, features = "pct_reads_in_peaks") +
  ggtitle("FRiP (%)")
p_depth <- FeaturePlot(atac, features = "nCount_peaks", max.cutoff = "q95") +
  ggtitle("Fragment count")

p_overview <- (p_clusters | p_samples | p_condition) / (p_tss | p_frip | p_depth)
ggsave(file.path(args$plot_dir, "integration_03_umap_overview.png"),
       p_overview, width = 18, height = 10)

# Per-sample batch mixing: each cluster should contain cells from all 4 samples
cluster_sample_table <- table(atac$cluster_final, atac$sample_id)
write.csv(as.data.frame.matrix(cluster_sample_table),
          file.path(args$outdir, "cluster_sample_composition.csv"))

# ── 10. Save outputs ──────────────────────────────────────────────────────────
write.csv(atac@meta.data, file.path(args$outdir, "integrated_cell_metadata.csv"))

out_rds <- file.path(args$outdir, "atac_integrated_clustered.rds")
saveRDS(atac, out_rds)
cat("Saved:", out_rds, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(args$outdir, "sessionInfo_integration.txt"))
cat("=== Integration complete ===\n")
