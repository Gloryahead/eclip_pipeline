/*
 * CELLRANGER_ATAC — align reads, correct barcodes, call peaks per sample
 *
 * Cell Ranger ATAC is proprietary (10x Genomics license required).
 * It cannot be distributed as a public container image.
 *
 * Setup: download the tarball from https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads
 * then set params.cellranger_atac_path in nextflow.config to the full path
 * of the cellranger-atac binary.
 *
 * Input : [sample_id, condition, fastq_dir/], reference_path
 * Output:
 *   cr_outs          — [sample_id, outs/]  (full Cell Ranger output directory)
 *   fragments        — [sample_id, fragments.tsv.gz]
 */
process CELLRANGER_ATAC {
    tag "${sample_id}"

    // No container — Cell Ranger ATAC is run from the binary on the filesystem.
    // Disable the container directive so Nextflow uses the system PATH.
    container null

    publishDir "${params.outdir}/02_cellranger/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(fastq_dir)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}/outs/"),              emit: cr_outs
    tuple val(sample_id), path("${sample_id}/outs/fragments.tsv.gz"), emit: fragments

    script:
    """
    # Run Cell Ranger ATAC alignment + peak calling
    ${params.cellranger_atac_path} count \\
        --id=${sample_id} \\
        --reference=${reference} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_id} \\
        --localcores=${task.cpus} \\
        --localmem=${(task.memory.toGiga() * 0.9).toInteger()}

    echo "Cell Ranger finished for ${sample_id}. Web summary:"
    ls -lh ${sample_id}/outs/web_summary.html
    echo "Key outputs:"
    ls -lh ${sample_id}/outs/fragments.tsv.gz
    ls -lh ${sample_id}/outs/filtered_peak_bc_matrix.h5
    ls -lh ${sample_id}/outs/singlecell.csv
    """
}
