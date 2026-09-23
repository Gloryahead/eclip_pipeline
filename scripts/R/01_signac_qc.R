#!/usr/bin/env Rscript
# =============================================================================
# 01_signac_qc.R — per-sample scATAC-seq quality control with Signac
#
# Purpose : Load Cell Ranger ATAC outputs, compute QC metrics, filter low-quality
#           barcodes, remove doublets with AMULET, filter blacklist peaks.
#
# Input   : Cell Ranger ATAC outs/ directory (one sample)
# Output  : <sample_id>_qc.rds — filtered Seurat/Signac object
#
# Reference genome : GRCh38 / hg38 (Ensembl 110 annotations via AnnotationHub)
# Tutorial source  : https://ngs101.com/scatac-seq-complete-beginners-guide-part-2
#
# Packages (versions pinned in containers/r_signac.def):
#   Signac 1.17.1, Seurat 5.5.1, Bioconductor 3.22, scDblFinder (AMULET)
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(AnnotationHub)
  library(EnsDb.Hsapiens.v86)   # fallback; primary annotation loaded from AnnotationHub
  library(rtracklayer)
  library(Rsamtools)
  library(scDblFinder)
  library(ggplot2)
  library(patchwork)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--sample_id",      required = TRUE)
parser$add_argument("--cr_outs",        required = TRUE, help = "Path to Cell Ranger outs/")
parser$add_argument("--fragments",      required = TRUE, help = "fragments.tsv.gz")
parser$add_argument("--frag_index",     required = TRUE, help = "fragments.tsv.gz.tbi")
parser$add_argument("--blacklist_bed",  required = TRUE)
parser$add_argument("--ensdb_cache",    required = TRUE, help = "AnnotationHub cache directory")
parser$add_argument("--min_fragments",  type = "integer", default = 3000)
parser$add_argument("--max_fragments",  type = "integer", default = 100000)
parser$add_argument("--min_tss",        type = "double",  default = 2.0)
parser$add_argument("--min_frip",       type = "double",  default = 15.0)
parser$add_argument("--max_nucleosome", type = "double",  default = 4.0)
parser$add_argument("--max_bl_ratio",   type = "double",  default = 0.05)
parser$add_argument("--amulet_q",       type = "double",  default = 0.05)
parser$add_argument("--outdir",         default = ".")
parser$add_argument("--plot_dir",       default = "plots")
args <- parser$parse_args()

set.seed(1234)
dir.create(args$outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(args$plot_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== QC:", args$sample_id, "===\n")

# ── 1. Load Cell Ranger outputs ───────────────────────────────────────────────
counts_h5 <- file.path(args$cr_outs, "filtered_peak_bc_matrix.h5")
metadata   <- read.csv(file.path(args$cr_outs, "singlecell.csv"),
                       header = TRUE, row.names = 1)

counts <- Read10X_h5(counts_h5)

# Remove non-standard chromosomes (scaffolds, patches) before creating assay
peak_ranges     <- GRanges(rownames(counts))
peaks_standard  <- seqnames(peak_ranges) %in% standardChromosomes(peak_ranges)
counts          <- counts[as.vector(peaks_standard), ]
cat("Peaks after standard-chromosome filter:", nrow(counts), "\n")

# ── 2. Create ChromatinAssay and Seurat object ────────────────────────────────
chrom_assay <- CreateChromatinAssay(
  counts      = counts,
  sep         = c(":", "-"),
  fragments   = args$fragments,
  min.cells   = 10,
  min.features = 200
)

atac <- CreateSeuratObject(
  counts    = chrom_assay,
  assay     = "peaks",
  meta.data = metadata
)

cat("Cells after initial load:", ncol(atac), "\n")

# ── 3. Genome annotation (Ensembl 110, chromosome names in UCSC format) ───────
setAnnotationHubOption("CACHE", args$ensdb_cache)
ah        <- AnnotationHub()
ens_query <- query(ah, c("EnsDb", "Homo sapiens", "110"))
ensdb     <- ah[[names(ens_query)[1]]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))  # UCSC format
genome(annotations)    <- "hg38"
Annotation(atac)       <- annotations

# ── 4. QC metrics ─────────────────────────────────────────────────────────────
# TSS enrichment — slow step; fast=FALSE gives the diagnostic coverage plot
cat("Computing TSS enrichment...\n")
atac <- TSSEnrichment(object = atac, fast = FALSE)

cat("Computing nucleosome signal...\n")
atac <- NucleosomeSignal(object = atac)

# FRiP: fraction of reads in peaks
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100

# Blacklist ratio
blacklist_hg38 <- rtracklayer::import(args$blacklist_bed, format = "bed")
blacklist_hg38 <- keepStandardChromosomes(blacklist_hg38, pruning.mode = "coarse")
genome(blacklist_hg38) <- "hg38"

atac$blacklist_ratio <- FractionCountsInRegion(
  object  = atac,
  assay   = "peaks",
  regions = blacklist_hg38
)

cat(sprintf("Metric summary for %s:\n", args$sample_id))
cat(sprintf("  passed_filters : %.0f – %.0f (median %.0f)\n",
    min(atac$passed_filters, na.rm=TRUE),
    max(atac$passed_filters, na.rm=TRUE),
    median(atac$passed_filters, na.rm=TRUE)))
cat(sprintf("  TSS enrichment : %.2f – %.2f (median %.2f)\n",
    min(atac$TSS.enrichment, na.rm=TRUE),
    max(atac$TSS.enrichment, na.rm=TRUE),
    median(atac$TSS.enrichment, na.rm=TRUE)))

# ── 5. Diagnostic plots (pre-filter) ─────────────────────────────────────────
# Nucleosome group for fragment histogram
atac$nucleosome_group <- ifelse(is.na(atac$nucleosome_signal), "Undefined",
  ifelse(atac$nucleosome_signal > args$max_nucleosome, "NS > 4", "NS < 4"))

p_frag_hist <- FragmentHistogram(object = atac, group.by = "nucleosome_group")
ggsave(file.path(args$plot_dir, paste0(args$sample_id, "_01_fragment_histogram.png")),
       p_frag_hist, width = 8, height = 4)

atac$tss_group <- ifelse(atac$TSS.enrichment > args$min_tss, "High TSS", "Low TSS")
p_tss <- TSSPlot(atac, group.by = "tss_group") + NoLegend()
ggsave(file.path(args$plot_dir, paste0(args$sample_id, "_02_tss_enrichment.png")),
       p_tss, width = 8, height = 4)

p_violin_pre <- VlnPlot(
  object   = atac,
  features = c("passed_filters", "TSS.enrichment", "pct_reads_in_peaks",
               "nucleosome_signal", "blacklist_ratio"),
  pt.size  = 0,
  ncol     = 5
)
ggsave(file.path(args$plot_dir, paste0(args$sample_id, "_03_qc_violins_prefilter.png")),
       p_violin_pre, width = 20, height = 5)

p_density <- DensityScatter(atac, x = "nCount_peaks", y = "TSS.enrichment",
                             log_x = TRUE, quantiles = TRUE)
ggsave(file.path(args$plot_dir, paste0(args$sample_id, "_04_density_scatter.png")),
       p_density, width = 6, height = 5)

# ── 6. Cell-level filtering ───────────────────────────────────────────────────
n_before <- ncol(atac)

# Report how many cells fail each criterion before filtering
filter_criteria <- list(
  low_depth      = atac$passed_filters <= args$min_fragments,
  high_depth     = atac$passed_filters >= args$max_fragments,
  low_tss        = atac$TSS.enrichment <= args$min_tss,
  low_frip       = atac$pct_reads_in_peaks <= args$min_frip,
  high_nucleosome = atac$nucleosome_signal >= args$max_nucleosome,
  high_blacklist = atac$blacklist_ratio >= args$max_bl_ratio
)
filter_summary <- data.frame(
  criterion   = names(filter_criteria),
  cells_failed = sapply(filter_criteria, sum, na.rm = TRUE)
)
cat("Cells failing each criterion:\n")
print(filter_summary)

atac_filtered <- subset(
  x      = atac,
  subset = passed_filters    >  args$min_fragments &
           passed_filters    <  args$max_fragments &
           TSS.enrichment    >  args$min_tss       &
           pct_reads_in_peaks > args$min_frip      &
           nucleosome_signal  < args$max_nucleosome &
           blacklist_ratio    < args$max_bl_ratio
)
n_after_qc <- ncol(atac_filtered)
cat(sprintf("Cells: %d → %d after QC filters (%.1f%% retained)\n",
    n_before, n_after_qc, 100 * n_after_qc / n_before))

# ── 7. AMULET doublet detection ───────────────────────────────────────────────
cat("Running AMULET doublet detection...\n")

frag_contigs <- seqnamesTabix(TabixFile(args$fragments))
scaffolds    <- setdiff(frag_contigs, paste0("chr", c(1:22, "X", "Y")))

exclude_regions <- suppressWarnings(c(
  GRanges(c("chrX", "chrY"),   IRanges(1L, width = 10^9)),
  GRanges(scaffolds,            IRanges(1L, width = 10^9)),
  blacklist_hg38
))

amulet_res <- amulet(
  args$fragments,
  regionsToExclude = exclude_regions,
  barcodes         = colnames(atac_filtered)
)

atac_filtered$amulet_q <- amulet_res[colnames(atac_filtered), "q.value"]
atac_filtered$amulet_q[is.na(atac_filtered$amulet_q)] <- 1   # NA → keep

atac_filtered <- subset(atac_filtered, subset = amulet_q >= args$amulet_q)
n_after_amulet <- ncol(atac_filtered)
cat(sprintf("Cells after AMULET: %d (removed %d doublets)\n",
    n_after_amulet, n_after_qc - n_after_amulet))

# ── 8. Feature-level filtering ────────────────────────────────────────────────
# Remove peaks overlapping blacklist regions
peaks_not_bl <- countOverlaps(granges(atac_filtered), blacklist_hg38) == 0
atac_filtered <- atac_filtered[as.vector(peaks_not_bl), ]

# Remove peaks detected in too few cells (rare peaks add noise, not signal)
peak_counts    <- LayerData(atac_filtered, assay = "peaks", layer = "counts")
cells_per_peak <- Matrix::rowSums(peak_counts > 0)
min_cells      <- max(10, ceiling(0.005 * ncol(atac_filtered)))
atac_filtered  <- atac_filtered[cells_per_peak >= min_cells, ]

cat(sprintf("Peaks remaining: %d\n", nrow(atac_filtered)))

# ── 9. Post-filter plots ──────────────────────────────────────────────────────
p_violin_post <- VlnPlot(
  object   = atac_filtered,
  features = c("passed_filters", "TSS.enrichment", "pct_reads_in_peaks",
               "nucleosome_signal", "blacklist_ratio"),
  pt.size  = 0,
  ncol     = 5
)
ggsave(file.path(args$plot_dir, paste0(args$sample_id, "_05_qc_violins_postfilter.png")),
       p_violin_post, width = 20, height = 5)

# ── 10. Add sample metadata and save ─────────────────────────────────────────
atac_filtered$sample_id  <- args$sample_id
atac_filtered$condition  <- ifelse(grepl("^severe", args$sample_id), "Severe", "Healthy")

out_rds <- file.path(args$outdir, paste0(args$sample_id, "_qc.rds"))
saveRDS(atac_filtered, out_rds)
cat("Saved:", out_rds, "\n")

# Save session info for reproducibility audit
writeLines(capture.output(sessionInfo()),
           file.path(args$outdir, paste0(args$sample_id, "_sessionInfo.txt")))

cat("=== QC complete:", args$sample_id, "===\n")
