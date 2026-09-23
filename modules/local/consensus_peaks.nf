/*
 * CONSENSUS_PEAKS — build a common peak set and requantify all samples against it (Part 3 phase 1)
 *
 * Why? Each sample's MACS3 peaks are slightly different. If you merge samples
 * before building a shared feature space, the count matrices are incomparable
 * (a "zero" in one sample might mean the peak was never measured, not closed).
 * This step creates one union peak set, then re-counts each sample's fragments
 * against those exact same regions using Signac's FeatureMatrix().
 *
 * Input : lists of [qc_rds, peak_narrowPeak, fragments.tsv.gz, fragments.tbi], blacklist
 * Output: merged Seurat/Signac RDS with consensus peak counts for all samples
 */
process CONSENSUS_PEAKS {
    tag "all samples"

    container params.r_signac_sif

    publishDir "${params.outdir}/05_integration", mode: 'copy'

    input:
    path qc_rds_list     // list of per-sample _qc.rds files
    path peak_bed_list   // list of per-sample MACS3 narrowPeak files
    path frag_list       // list of fragments.tsv.gz files
    path frag_idx_list   // list of fragments.tsv.gz.tbi files
    path blacklist_bed

    output:
    path "atac_merged_consensus.rds", emit: merged_rds
    path "consensus_peaks.rds"
    path "plots/consensus_*.png"

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/scripts/R/02_consensus_peaks.R \\
        --qc_rds_list   "${qc_rds_list.join(',')}" \\
        --peak_bed_list "${peak_bed_list.join(',')}" \\
        --frag_list     "${frag_list.join(',')}" \\
        --frag_idx_list "${frag_idx_list.join(',')}" \\
        --blacklist_bed ${blacklist_bed} \\
        --peak_min_width ${params.peak_min_width} \\
        --peak_max_width ${params.peak_max_width} \\
        --min_cells      ${params.min_cells_per_peak} \\
        --outdir         . \\
        --plot_dir       plots
    """
}
