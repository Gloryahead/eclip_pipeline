/*
 * CELL_TYPE_ANNOTATION — gene activity scoring, label transfer, motif enrichment (Part 4)
 *
 * Three independent evidence streams per cluster:
 *   1. Label transfer from matched scRNA-seq (CCA cross-modality)
 *   2. SingleR bulk immune reference (Monaco, 29 sorted populations)
 *   3. Motif enrichment on cluster marker peaks (JASPAR 2020 core)
 *
 * Using multiple methods guards against reference bias. Agreement between
 * streams gives high confidence; disagreement flags uncertain clusters.
 *
 * Input : clustered RDS, matched RNA reference RDS (pbmc_10k_v3 from Satija lab)
 * Output: annotated RDS + CSV tables + PNG figures
 */
process CELL_TYPE_ANNOTATION {
    tag "all samples"

    container params.r_signac_sif

    publishDir "${params.outdir}/06_annotation", mode: 'copy'

    input:
    path clustered_rds
    path rna_ref_rds

    output:
    path "atac_annotated.rds",              emit: annotated_rds
    path "label_transfer_scores.csv"
    path "marker_peaks_annotated.csv"
    path "motif_enrichment.csv"
    path "plots/annotation_*.png"

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/scripts/R/04_cell_type_annotation.R \\
        --clustered_rds      ${clustered_rds} \\
        --rna_ref_rds        ${rna_ref_rds} \\
        --harmony_dims       "${params.harmony_dims}" \\
        --gene_upstream      ${params.gene_activity_upstream} \\
        --pred_score_min     ${params.prediction_score_min} \\
        --marker_fdr         ${params.marker_peak_fdr} \\
        --genome             ${params.genome} \\
        --outdir             . \\
        --plot_dir           plots
    """
}
