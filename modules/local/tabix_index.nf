/*
 * TABIX_INDEX — index the fragments.tsv.gz file for random-access queries
 *
 * Signac requires the .tbi index to sit alongside the fragments file.
 * The fragments file must be block-compressed (bgzip) — Cell Ranger produces
 * it this way already, so we only need to index it.
 *
 * Input : [sample_id, fragments.tsv.gz]
 * Output: [sample_id, fragments.tsv.gz.tbi]  (index only — fragments already published)
 */
process TABIX_INDEX {
    tag "${sample_id}"

    container 'quay.io/biocontainers/htslib:1.21--h5efdd21_0'

    publishDir "${params.outdir}/02_cellranger/${sample_id}/outs", mode: 'copy'

    input:
    tuple val(sample_id), path(fragments)

    output:
    tuple val(sample_id), path("${fragments}.tbi"), emit: indexed_frags

    script:
    """
    tabix -p bed ${fragments}
    echo "Indexed: ${fragments}.tbi"
    """
}
