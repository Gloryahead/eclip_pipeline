#!/usr/bin/env Rscript
# =============================================================================
# 04_cell_type_annotation.R — cell type identification for scATAC-seq clusters
#
# Purpose : Assign cell type labels to scATAC-seq clusters using three
#           independent evidence streams and motif enrichment analysis.
#
# Evidence streams (converging agreement = high confidence):
#   1. Gene activity + label transfer from matched scRNA-seq (CCA method)
#   2. SingleR bulk immune reference (Monaco, 29 sorted populations)
#   3. Curated marker gene panel scoring (AddModuleScore)
#   4. JASPAR 2020 motif enrichment on cluster marker peaks
#
# Key technical notes:
#   - Gene activity is computed by counting fragments near gene promoters
#     (2 kb upstream); it's a proxy for expression, not direct measurement.
#   - CCA (canonical correlation analysis) is required for cross-modality
#     label transfer; standard PCA fails when RNA and ATAC occupy different
#     feature spaces.
#   - Cluster-level aggregation is necessary because individual scATAC-seq
#     cells are too sparse for reliable single-cell annotation.
#   - Motif families (e.g., T-box, PAX) are more interpretable than individual
#     factors because TF families share nearly identical binding sequences.
#
# Input   : atac_integrated_clustered.rds, pbmc RNA reference RDS
# Output  : atac_annotated.rds, CSV tables, PNG plots
#
# Tutorial: https://ngs101.com/scatac-seq-complete-beginners-guide-part-4
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(SingleR)
  library(celldex)
  library(motifmatchr)
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--clustered_rds",  required = TRUE)
parser$add_argument("--rna_ref_rds",    required = TRUE,
                    help = "Annotated scRNA-seq Seurat object (pbmc_10k_v3)")
parser$add_argument("--harmony_dims",   default = "2:30")
parser$add_argument("--gene_upstream",  type = "integer", default = 2000)
parser$add_argument("--pred_score_min", type = "double",  default = 0.5)
parser$add_argument("--marker_fdr",     type = "double",  default = 0.01)
parser$add_argument("--genome",         default = "hg38")
parser$add_argument("--outdir",         default = ".")
parser$add_argument("--plot_dir",       default = "plots")
args <- parser$parse_args()

set.seed(1234)
dir.create(args$outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

parse_dims <- function(s) {
  if (grepl(":", s)) { parts <- as.integer(strsplit(s, ":")[[1]]); seq(parts[1], parts[2]) }
  else as.integer(strsplit(s, ",")[[1]])
}
harmony_dims <- parse_dims(args$harmony_dims)

cat("=== Cell Type Annotation ===\n")

# ── 1. Load objects ───────────────────────────────────────────────────────────
atac    <- readRDS(args$clustered_rds)
pbmc_rna <- readRDS(args$rna_ref_rds)

cat("ATAC cells:", ncol(atac), "| Clusters:", nlevels(Idents(atac)), "\n")
cat("RNA reference cells:", ncol(pbmc_rna), "\n")

# ── 2. Gene activity scoring ──────────────────────────────────────────────────
# Counts ATAC fragments within gene bodies + promoter regions as a proxy
# for transcriptional activity. Normalized with median depth scaling.
cat("Computing gene activity matrix...\n")
gene_activities <- GeneActivity(
  object            = atac,
  extend.upstream   = args$gene_upstream,
  extend.downstream = 0
)

atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
atac <- NormalizeData(
  object               = atac,
  assay                = "ACTIVITY",
  normalization.method = "LogNormalize",
  scale.factor         = median(atac$nCount_ACTIVITY)
)
atac <- ScaleData(atac, assay = "ACTIVITY",
                  features = rownames(atac[["ACTIVITY"]]))

# ── 3. Label transfer from scRNA-seq (CCA cross-modality) ─────────────────────
# CCA finds correlated gene-activity/expression axes across modalities.
# Standard PCA cannot do this because the two modalities have different features.
cat("Finding transfer anchors (CCA)...\n")

DefaultAssay(pbmc_rna) <- "RNA"
shared_features <- intersect(
  rownames(atac[["ACTIVITY"]]),
  rownames(pbmc_rna)
)
cat("Shared features for transfer:", length(shared_features), "\n")

transfer_anchors <- FindTransferAnchors(
  reference       = pbmc_rna,
  query           = atac,
  features        = shared_features,
  reference.assay = "RNA",
  query.assay     = "ACTIVITY",
  reduction       = "cca",
  dims            = 1:30
)

predicted_labels <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = pbmc_rna$celltype,
  weight.reduction = atac[["harmony"]],
  dims             = seq_along(harmony_dims)
)

atac <- AddMetaData(atac, metadata = predicted_labels)
write.csv(predicted_labels,
          file.path(args$outdir, "label_transfer_scores.csv"))

# Cluster-level majority vote for predicted labels
cluster_labels <- atac@meta.data %>%
  group_by(cluster_final) %>%
  summarise(
    predicted_celltype = names(sort(table(predicted.id), decreasing = TRUE))[1],
    median_score       = median(predicted.score, na.rm = TRUE),
    high_conf_fraction = mean(predicted.score > args$pred_score_min, na.rm = TRUE),
    n_cells            = n()
  )
write.csv(cluster_labels, file.path(args$outdir, "cluster_label_transfer.csv"), row.names = FALSE)
cat("Label transfer complete. High-confidence fraction per cluster:\n")
print(cluster_labels[, c("cluster_final", "predicted_celltype", "high_conf_fraction")])

# ── 4. SingleR bulk immune reference annotation ───────────────────────────────
# Independent cross-check using Monaco bulk RNA-seq reference.
# Aggregated at the cluster level because single cells are too sparse.
cat("Running SingleR (Monaco immune reference)...\n")
monaco_ref <- celldex::fetchReference("monaco_immune", "2024-02-26")

activity_mat <- GetAssayData(atac, assay = "ACTIVITY", layer = "data")
singler_res <- SingleR(
  test     = activity_mat,
  ref      = monaco_ref,
  labels   = monaco_ref$label.main,
  clusters = atac$cluster_final
)

singler_df <- data.frame(
  cluster_final  = rownames(singler_res),
  singler_label  = singler_res$labels,
  singler_score  = apply(singler_res$scores, 1, max)
)
write.csv(singler_df, file.path(args$outdir, "singler_cluster_labels.csv"), row.names = FALSE)

# ── 5. Marker gene panel scoring ──────────────────────────────────────────────
marker_panels <- list(
  T_cell   = c("CD3D", "CD3E", "CD3G", "LCK", "IL7R", "THEMIS"),
  CD4_T    = c("CD4", "CD40LG", "MAL", "TRAT1"),
  CD8_T    = c("CD8A", "CD8B", "GZMK", "LINC02446"),
  NK_cell  = c("NCAM1", "KLRD1", "NKG7", "GNLY", "KLRF1", "PRF1"),
  B_cell   = c("MS4A1", "CD79A", "CD79B", "BANK1", "PAX5", "EBF1"),
  Monocyte = c("LYZ", "CD14", "VCAN", "S100A8", "CSF1R", "FCN1"),
  DC       = c("CLEC9A", "CLEC4C", "IL3RA", "FLT3", "LILRA4"),
  Platelet = c("PPBP", "PF4", "GP9", "ITGA2B")
)

# Filter to genes actually present in the ACTIVITY assay
marker_panels <- lapply(marker_panels, function(genes) {
  intersect(genes, rownames(atac[["ACTIVITY"]]))
})
marker_panels <- Filter(function(x) length(x) >= 2, marker_panels)

atac <- AddModuleScore(
  object   = atac,
  features = marker_panels,
  assay    = "ACTIVITY",
  name     = "panel",
  ctrl     = 50,
  seed     = 1234
)

# ── 6. Plot annotation overview ───────────────────────────────────────────────
p_labels <- DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
  ggtitle("Predicted cell types (label transfer)") + NoLegend()
ggsave(file.path(args$plot_dir, "annotation_01_label_transfer_umap.png"),
       p_labels, width = 10, height = 8)

p_score <- FeaturePlot(atac, features = "predicted.score") +
  ggtitle("Label transfer confidence score")
ggsave(file.path(args$plot_dir, "annotation_02_prediction_scores.png"),
       p_score, width = 7, height = 6)

# ── 7. Cluster marker peak discovery ─────────────────────────────────────────
# LR test with depth covariate eliminates sequencing-depth confounding.
# 500-cell cap per cluster speeds up computation without losing power.
cat("Finding cluster marker peaks (logistic regression)...\n")
DefaultAssay(atac) <- "peaks"

all_markers <- FindAllMarkers(
  object              = atac,
  test.use            = "LR",
  latent.vars         = "nCount_peaks",
  min.pct             = 0.1,
  only.pos            = TRUE,
  max.cells.per.ident = 500,
  verbose             = FALSE
)

# Annotate each marker peak with the closest gene
all_markers_sig <- all_markers[all_markers$p_val_adj < args$marker_fdr, ]
if (nrow(all_markers_sig) > 0) {
  closest_genes <- ClosestFeature(atac, regions = all_markers_sig$gene)
  all_markers_sig <- cbind(all_markers_sig, closest_genes)
}
write.csv(all_markers_sig,
          file.path(args$outdir, "marker_peaks_annotated.csv"), row.names = FALSE)
cat("Significant marker peaks:", nrow(all_markers_sig), "\n")

# ── 8. Motif enrichment analysis ─────────────────────────────────────────────
cat("Adding JASPAR 2020 motifs...\n")
pfm <- getMatrixSet(
  x    = JASPAR2020,
  opts = list(species = 9606, collection = "CORE", all_versions = FALSE)
)

atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm    = pfm
)

# Find enriched motifs in each cluster's marker peaks
motif_list <- list()
for (clust in unique(all_markers_sig$cluster)) {
  clust_peaks <- all_markers_sig$gene[all_markers_sig$cluster == clust]
  if (length(clust_peaks) >= 10) {   # need enough peaks for reliable enrichment
    enr <- tryCatch(
      FindMotifs(object = atac, features = clust_peaks),
      error = function(e) { cat("Motif enrichment failed for cluster", clust, ":", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(enr)) {
      enr$cluster <- clust
      motif_list[[as.character(clust)]] <- enr
    }
  }
}
if (length(motif_list) > 0) {
  motif_all <- do.call(rbind, motif_list)
  write.csv(motif_all, file.path(args$outdir, "motif_enrichment.csv"), row.names = FALSE)
  cat("Motif enrichment saved for", length(motif_list), "clusters\n")
}

# ── 9. Combine annotations and finalize cell type labels ─────────────────────
# Merge transfer and SingleR labels; use transfer as primary, SingleR as backup
cluster_annotation <- merge(cluster_labels, singler_df, by = "cluster_final", all = TRUE)
write.csv(cluster_annotation,
          file.path(args$outdir, "cluster_annotation_combined.csv"), row.names = FALSE)

# Add final label to cell metadata (edit this table after reviewing the plots)
label_map <- setNames(cluster_annotation$predicted_celltype,
                      cluster_annotation$cluster_final)
atac$cell_type <- label_map[as.character(atac$cluster_final)]

p_final <- DimPlot(atac, group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Final cell type annotation") + NoLegend()
ggsave(file.path(args$plot_dir, "annotation_03_final_cell_types.png"),
       p_final, width = 10, height = 8)

# ── 10. Save annotated object ─────────────────────────────────────────────────
out_rds <- file.path(args$outdir, "atac_annotated.rds")
saveRDS(atac, out_rds)
cat("Saved annotated object:", out_rds, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(args$outdir, "sessionInfo_annotation.txt"))
cat("=== Cell type annotation complete ===\n")
