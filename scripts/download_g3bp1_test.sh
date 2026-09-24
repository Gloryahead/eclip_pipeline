#!/usr/bin/env bash
# =============================================================================
# download_g3bp1_test.sh
# Downloads 4 FASTQ files for the G3BP1 pipeline test run.
#
# Dataset : He et al. 2021, Nucleic Acids Research (GSE168943)
# Samples : G3BP1 IP + SMInput, Ctrl/DMSO condition, 2 replicates each
# Total   : ~6.1 GB compressed
#
# Usage:
#   bash scripts/download_g3bp1_test.sh
#
# Prerequisites:
#   eclip conda environment activated (has sra-tools + pigz)
# =============================================================================
set -euo pipefail

FASTQ_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq/g3bp1_test"
THREADS="${THREADS:-8}"
TMP_PREFETCH="/tmp/sra_prefetch_g3bp1_$$"

mkdir -p "${FASTQ_DIR}" "${TMP_PREFETCH}"
trap 'rm -rf "${TMP_PREFETCH}"' EXIT

download_one() {
    local label="$1" srr="$2"
    local outfile="${FASTQ_DIR}/${label}.fastq.gz"

    if [[ -f "${outfile}" ]]; then
        echo "[SKIP] ${label}: already downloaded."
        return 0
    fi

    echo "[DOWN] ${label}  (${srr})"
    prefetch --max-size 50G -p "${srr}" -O "${TMP_PREFETCH}/"

    fasterq-dump "${TMP_PREFETCH}/${srr}/${srr}.sra" \
        --outdir "${FASTQ_DIR}" \
        --temp  "${TMP_PREFETCH}" \
        --threads "${THREADS}" \
        --split-3

    if [[ -f "${FASTQ_DIR}/${srr}.fastq" ]]; then
        mv "${FASTQ_DIR}/${srr}.fastq" "${FASTQ_DIR}/${label}.fastq"
    elif [[ -f "${FASTQ_DIR}/${srr}_1.fastq" ]]; then
        mv "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${label}.fastq"
        rm -f "${FASTQ_DIR}/${srr}_2.fastq"
    fi

    pigz -p "${THREADS}" "${FASTQ_DIR}/${label}.fastq"
    rm -rf "${TMP_PREFETCH}/${srr}"
    echo "[DONE] ${label}"
}

echo "=== Downloading G3BP1 test data (GSE168943) ==="
echo "Output: ${FASTQ_DIR}"
echo ""

# Ctrl condition — SMInput
download_one ctrl_input_rep1 SRR13966667  # SRX10344592  ~1.2 GB
download_one ctrl_input_rep2 SRR13966669  # SRX10344594  ~2.5 GB

# Ctrl condition — IP
download_one ctrl_ip_rep1    SRR13966668  # SRX10344593  ~1.1 GB
download_one ctrl_ip_rep2    SRR13966670  # SRX10344595  ~1.3 GB

echo ""
echo "=== Download complete ==="
echo "Files in: ${FASTQ_DIR}"
echo ""
echo "NEXT: sbatch g3bp1_test.slurm"
