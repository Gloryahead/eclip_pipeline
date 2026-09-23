/*
 * RENAME_FASTQS — detect read lengths and rename SRA files to Cell Ranger format
 *
 * 10x scATAC-seq read structure from SRA:
 *   8 bp  → I1 (sample index)
 *   16 bp → R2 (cell barcode — SRA stores I2/barcode as R2)
 *   50 bp → R1 (genomic read 1)
 *   50 bp → R3 (genomic read 2)
 *
 * Cell Ranger ATAC requires: <sample>_S1_L001_{I1,R1,R2,R3}_001.fastq.gz
 *
 * Input : [sample_id, srr, condition, raw_fastq_dir/]
 * Output: [sample_id, condition, renamed_fastq_dir/]
 */
process RENAME_FASTQS {
    tag "${sample_id}"

    container 'quay.io/biocontainers/pigz:2.8'

    publishDir "${params.outdir}/01_raw_fastqs/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(srr), val(condition), path(raw_dir)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_fastqs/"), emit: named_fastqs

    script:
    """
    #!/usr/bin/env python3
    import os, subprocess, shutil

    raw_dir   = "${raw_dir}"
    sample_id = "${sample_id}"
    srr       = "${srr}"
    out_dir   = f"{sample_id}_fastqs"
    os.makedirs(out_dir, exist_ok=True)

    # Measure the read length of each FASTQ by looking at the second line of the file
    def read_length(fq_gz):
        result = subprocess.run(
            f"zcat {fq_gz} | head -2 | tail -1",
            shell=True, capture_output=True, text=True
        )
        return len(result.stdout.strip())

    # Gather all FASTQ files for this SRR
    files = sorted(f for f in os.listdir(raw_dir) if f.startswith(srr) and f.endswith(".fastq.gz"))
    print(f"Found {len(files)} FASTQ files: {files}")

    lengths = {}
    for f in files:
        path = os.path.join(raw_dir, f)
        ln = read_length(path)
        lengths[f] = ln
        print(f"  {f}: {ln} bp")

    # Sort files into their roles by length
    # 8 bp = Illumina sample index → I1
    # 16 bp = 10x cell barcode     → R2 (Cell Ranger ATAC calls it R2)
    # 50 bp = genomic reads         → R1 and R3 (sorted by filename)
    i1_src = r1_src = r2_src = r3_src = None
    genomic_50 = []

    for f, ln in lengths.items():
        src = os.path.join(raw_dir, f)
        if ln == 8:
            i1_src = src
        elif ln == 16:
            r2_src = src
        else:
            genomic_50.append(src)

    genomic_50.sort()   # alphabetical order gives _1 before _2 which is R1 before R3
    if len(genomic_50) >= 2:
        r1_src, r3_src = genomic_50[0], genomic_50[1]

    # Copy with Cell Ranger naming
    targets = {
        "I1": i1_src,
        "R1": r1_src,
        "R2": r2_src,
        "R3": r3_src,
    }
    for read, src in targets.items():
        if src:
            dst = os.path.join(out_dir, f"{sample_id}_S1_L001_{read}_001.fastq.gz")
            shutil.copy2(src, dst)
            print(f"  {os.path.basename(src)} → {os.path.basename(dst)}")
        else:
            print(f"  WARNING: no file identified for {read}")

    print("\\nFinal FASTQs for Cell Ranger:")
    for f in sorted(os.listdir(out_dir)):
        print(f"  {f}")
    """
}
