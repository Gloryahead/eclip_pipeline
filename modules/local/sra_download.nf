/*
 * SRA_DOWNLOAD — download FASTQ files from NCBI SRA
 *
 * Input : [sample_id, srr_accession, condition]
 * Output: [sample_id, srr_accession, condition, raw_fastq_dir/]
 *
 * --include-technical is REQUIRED for 10x scATAC-seq: it preserves the 16 bp
 * barcode read that SRA would otherwise silently discard as "technical."
 */
process SRA_DOWNLOAD {
    tag "${sample_id} (${srr})"

    container 'quay.io/biocontainers/sra-tools:3.1.1--h9f5acd7_0'

    publishDir "${params.outdir}/01_raw_fastqs/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), val(srr), val(condition)

    output:
    tuple val(sample_id), val(srr), val(condition), path("${sample_id}_raw/"), emit: raw_fastqs

    script:
    """
    mkdir -p ${sample_id}_raw

    # Step 1: Download the SRA archive (up to 100 GB per file)
    prefetch --max-size 100G --output-directory . ${srr}

    # Step 2: Convert to FASTQ — keep ALL reads including the barcode (technical) read
    fasterq-dump \\
        --split-files \\
        --include-technical \\
        --threads ${task.cpus} \\
        --outdir ${sample_id}_raw \\
        ${srr}/${srr}.sra

    # Step 3: Compress with pigz (parallel gzip, much faster than gzip)
    pigz -p ${task.cpus} ${sample_id}_raw/*.fastq

    echo "Downloaded and compressed FASTQs for ${sample_id}:"
    ls -lh ${sample_id}_raw/
    """
}
