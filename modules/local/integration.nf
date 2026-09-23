/*
 * INTEGRATION — TF-IDF normalization, LSI, Harmony batch correction, clustering (Part 3)
 *
 * Key decisions made here:
 *   - TF-IDF method 1 (log-transformed) — published standard for scATAC-seq
 *   - SVD component 1 is EXCLUDED from all downstream steps because it correlates
 *     almost perfectly with sequencing depth (not biology)
 *   - Harmony corrects for sample-of-origin batch effects on LSI embeddings
 *   - SLM clustering algorithm (algorithm=3) — recommended for sparse chromatin graphs
 *   - clustree sweeps multiple resolutions so you can pick the stable one
 *
 * Input : merged consensus-peak RDS (from CONSENSUS_PEAKS)
 * Output: clustered RDS with UMAP embeddings and cluster assignments
 */
process INTEGRATION {
    tag "all samples"

    container params.r_signac_sif

    publishDir "${params.outdir}/05_integration", mode: 'copy'

    input:
    path merged_rds

    output:
    path "atac_integrated_clustered.rds",  emit: clustered_rds
    path "integrated_cell_metadata.csv"
    path "plots/integration_*.png"

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/scripts/R/03_integration_clustering.R \\
        --merged_rds        ${merged_rds} \\
        --lsi_dims          "${params.lsi_dims}" \\
        --harmony_dims      "${params.harmony_dims}" \\
        --resolutions       "${params.cluster_resolution}" \\
        --min_cells         ${params.min_cells_per_peak} \\
        --outdir            . \\
        --plot_dir          plots
    """
}
