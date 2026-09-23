/*
 * SIGNAC_QC — per-sample quality control and doublet removal (Part 2)
 *
 * Runs scripts/R/01_signac_qc.R inside the R/Signac container.
 * Outputs a filtered Seurat/Signac RDS object ready for integration.
 *
 * QC metrics computed:
 *   - TSS enrichment (signal-to-noise; threshold > 2)
 *   - Nucleosome signal (transposition efficiency; threshold < 4)
 *   - FRiP — fraction of reads in peaks (threshold > 15%)
 *   - Blacklist ratio (artifact peaks; threshold < 5%)
 *   - AMULET doublet detection (q-value threshold 0.05)
 *
 * Input : [sample_id, outs/, fragments.tsv.gz.tbi], blacklist_bed, ensdb_cache/
 * Output: [sample_id, <sample_id>_qc.rds]
 */
process SIGNAC_QC {
    tag "${sample_id}"

    container params.r_signac_sif

    publishDir "${params.outdir}/04_qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(cr_outs), path(frag_index)
    path blacklist_bed
    path ensdb_cache

    output:
    tuple val(sample_id), path("${sample_id}_qc.rds"), emit: qc_rds
    path "plots/${sample_id}_*.png"

    script:
    """
    mkdir -p plots

    Rscript ${projectDir}/scripts/R/01_signac_qc.R \\
        --sample_id        ${sample_id} \\
        --cr_outs          ${cr_outs} \\
        --fragments        ${cr_outs}/fragments.tsv.gz \\
        --frag_index       ${frag_index} \\
        --blacklist_bed    ${blacklist_bed} \\
        --ensdb_cache      ${ensdb_cache} \\
        --min_fragments    ${params.min_fragments} \\
        --max_fragments    ${params.max_fragments} \\
        --min_tss          ${params.min_tss} \\
        --min_frip         ${params.min_frip} \\
        --max_nucleosome   ${params.max_nucleosome} \\
        --max_bl_ratio     ${params.max_blacklist_ratio} \\
        --amulet_q         ${params.amulet_qvalue} \\
        --outdir           . \\
        --plot_dir         plots
    """
}
