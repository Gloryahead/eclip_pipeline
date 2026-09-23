/*
 * MACS3_CALLPEAK — call open-chromatin peaks from the Cell Ranger fragments file
 *
 * Why re-call peaks with MACS3?
 *   Cell Ranger's wavelet peak caller tends to produce wider peaks; MACS3 gives
 *   narrower, more precise peaks useful for building a tight consensus peak set
 *   across samples in Part 3.
 *
 * -f FRAG reads the 10x fragment format directly (no BAM needed).
 * --keep-dup all is correct for scATAC-seq because UMIs are not used; each
 * fragment line is already a deduplicated molecule from Cell Ranger.
 *
 * Input : [sample_id, fragments.tsv.gz]
 * Output: [sample_id, <sample_id>_peaks.narrowPeak]
 */
process MACS3_CALLPEAK {
    tag "${sample_id}"

    container 'quay.io/biocontainers/macs3:3.0.2--py312h66b675f_0'

    publishDir "${params.outdir}/03_macs3/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fragments)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), emit: peak_bed

    script:
    """
    macs3 callpeak \\
        -f FRAG \\
        -t ${fragments} \\
        -g hs \\
        -n ${sample_id} \\
        --outdir . \\
        --keep-dup all \\
        -q 0.01

    echo "MACS3 peaks called for ${sample_id}:"
    wc -l ${sample_id}_peaks.narrowPeak
    """
}
