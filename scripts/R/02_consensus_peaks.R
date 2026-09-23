#!/usr/bin/env Rscript
# =============================================================================
# 02_consensus_peaks.R — build consensus peak set and requantify all samples
#
# Purpose : Merge per-sample MACS3 peaks into one union peak set, filter it,
#           then recount each sample's fragments against those identical regions
#           with Signac's FeatureMatrix(). Merge all samples into one object.
#
# Why this step is necessary:
#   Each sample's MACS3 call produces slightly different peak coordinates.
#   Naively merging count matrices from different feature spaces produces
#   incomparable columns. By quantifying all samples against the SAME peaks,
#   every zero is a true "zero detected" rather than "not measured here."
#
# Input   : per-sample *_qc.rds, *_peaks.narrowPeak, fragments.tsv.gz, .tbi
# Output  : atac_merged_consensus.rds — merged Seurat/Signac object
#
# Tutorial: https://ngs101.com/scatac-seq-complete-beginners-guide-part-3
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(ggplot2)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--qc_rds_list",    required = TRUE, help = "Comma-separated list of *_qc.rds")
parser$add_argument("--peak_bed_list",  required = TRUE, help = "Comma-separated list of narrowPeak files")
parser$add_argument("--frag_list",      required = TRUE, help = "Comma-separated fragments.tsv.gz")
parser$add_argument("--frag_idx_list",  required = TRUE, help = "Comma-separated fragments.tsv.gz.tbi")
parser$add_argument("--blacklist_bed",  required = TRUE)
parser$add_argument("--peak_min_width", type = "integer", default = 20)
parser$add_argument("--peak_max_width", type = "integer", default = 10000)
parser$add_argument("--min_cells",      type = "integer", default = 20)
parser$add_argument("--outdir",         default = ".")
parser$add_argument("--plot_dir",       default = "plots")
args <- parser$parse_args()

set.seed(1234)
dir.create(args$outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

# Parse comma-separated lists
qc_rds_files   <- trimws(strsplit(args$qc_rds_list,   ",")[[1]])
peak_bed_files  <- trimws(strsplit(args$peak_bed_list, ",")[[1]])
frag_files      <- trimws(strsplit(args$frag_list,     ",")[[1]])
frag_idx_files  <- trimws(strsplit(args$frag_idx_list, ",")[[1]])

cat("Samples found:", length(qc_rds_files), "\n")

# ── 1. Load blacklist ─────────────────────────────────────────────────────────
blacklist <- rtracklayer::import(args$blacklist_bed, format = "bed")
blacklist <- keepStandardChromosomes(blacklist, pruning.mode = "coarse")
genome(blacklist) <- "hg38"

# ── 2. Build consensus peak set ───────────────────────────────────────────────
cat("Reading per-sample peak files...\n")
peak_list <- lapply(peak_bed_files, function(f) {
  gr <- rtracklayer::import(f, format = "narrowPeak")
  keepStandardChromosomes(gr, pruning.mode = "coarse")
})

# Merge all peaks into a union set with reduce() (merges overlapping intervals)
all_peaks <- Reduce(c, peak_list)
consensus  <- reduce(all_peaks)

cat("Total peaks before filters:", length(consensus), "\n")

# Filter 1: Remove peaks outside 20–10,000 bp (fusion artifacts and tiny noise)
w <- width(consensus)
consensus <- consensus[w >= args$peak_min_width & w <= args$peak_max_width]
cat("After width filter:", length(consensus), "\n")

# Filter 2: Standard chromosomes only
consensus <- keepStandardChromosomes(consensus, pruning.mode = "coarse")

# Filter 3: Remove blacklist-overlapping peaks
consensus <- subsetByOverlaps(consensus, blacklist, invert = TRUE)
cat("After blacklist filter:", length(consensus), "\n")

saveRDS(consensus, file.path(args$outdir, "consensus_peaks.rds"))

# ── 3. Load per-sample QC objects and requantify against consensus peaks ───────
cat("Requantifying samples against consensus peaks...\n")
sample_objects <- list()

for (i in seq_along(qc_rds_files)) {
  atac_i <- readRDS(qc_rds_files[i])
  sample_id <- unique(atac_i$sample_id)
  cat("  Processing:", sample_id, "\n")

  # Rebuild fragment object pointing to this sample's fragments file
  frag_obj <- CreateFragmentObject(
    path  = frag_files[i],
    cells = colnames(atac_i)
  )

  # Count fragments in consensus peaks for this sample's QC-passed barcodes
  new_counts <- FeatureMatrix(
    fragments = frag_obj,
    features  = consensus,
    cells     = colnames(atac_i)
  )

  # Rebuild ChromatinAssay with new counts, carrying over annotations
  chrom_assay <- CreateChromatinAssay(
    counts    = new_counts,
    sep       = c(":", "-"),
    fragments = frag_files[i]
  )
  Annotation(chrom_assay) <- Annotation(atac_i)

  # Build new Seurat object, preserving all QC metadata columns
  atac_new <- CreateSeuratObject(
    counts    = chrom_assay,
    assay     = "peaks",
    meta.data = atac_i@meta.data
  )

  sample_objects[[sample_id]] <- atac_new
  cat("    Cells:", ncol(atac_new), "| Peaks:", nrow(atac_new), "\n")
}

# ── 4. Merge all samples into one object ─────────────────────────────────────
# add.cell.ids prefixes barcodes with sample names to prevent collisions
cat("Merging samples...\n")
atac_merged <- Reduce(function(a, b) merge(a, y = b), sample_objects)
cat("Merged object — Cells:", ncol(atac_merged), "| Peaks:", nrow(atac_merged), "\n")

# ── 5. Remove rare peaks after merging (fewer cells → higher threshold) ───────
peak_counts    <- LayerData(atac_merged, assay = "peaks", layer = "counts")
cells_per_peak <- Matrix::rowSums(peak_counts > 0)
peaks_keep     <- cells_per_peak >= args$min_cells
atac_merged    <- atac_merged[as.vector(peaks_keep), ]
cat("Peaks after rare-peak filter:", nrow(atac_merged), "\n")

# ── 6. Cell composition plot ─────────────────────────────────────────────────
cell_counts <- table(atac_merged$sample_id)
cat("Cells per sample:\n"); print(cell_counts)

df_cells <- as.data.frame(cell_counts)
colnames(df_cells) <- c("sample_id", "n_cells")
p_cells <- ggplot(df_cells, aes(x = sample_id, y = n_cells, fill = sample_id)) +
  geom_col() +
  labs(title = "Cells per sample after QC", x = NULL, y = "Cells") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(args$plot_dir, "consensus_01_cells_per_sample.png"),
       p_cells, width = 6, height = 4)

# ── 7. Save merged object ─────────────────────────────────────────────────────
out_rds <- file.path(args$outdir, "atac_merged_consensus.rds")
saveRDS(atac_merged, out_rds)
cat("Saved:", out_rds, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(args$outdir, "sessionInfo_consensus.txt"))
cat("=== Consensus peaks complete ===\n")
