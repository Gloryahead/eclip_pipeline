#!/usr/bin/env bash
# =============================================================================
# download_geo_fastqs.sh
# Downloads all 31 eCLIP IP FASTQ files for Rhine et al. 2025 from NCBI SRA.
# GEO accession: GSE276985 | BioProject: PRJNA1159891
#
# Usage:
#   bash scripts/download_geo_fastqs.sh [human|mouse|all]   (default: all)
#
# ⚠  SMInput (size-matched input) samples are NOT deposited in GSE276985.
#    Only the 31 IP samples are publicly available. Before you can run Skipper
#    you need paired input FASTQ files. Options:
#      a) Email the corresponding author to request input FASTQ files.
#      b) Search the ENCODE portal (encodeproject.org) for TDP-43/G3BP1/Caprin1
#         eCLIP in NGN2 or K562 cells — Yeo Lab ENCODE deposits may include
#         size-matched inputs from the same protocol.
#      c) Generate your own size-matched input library using ENCODE4 eCLIP SOP.
#    Fill in Input_fastq paths in config/manifest_human.csv and
#    config/manifest_mouse.csv before running the pipeline.
#
# SRX → SRR mapping verified 2026-09-15 against NCBI SRA runinfo API.
# =============================================================================
set -euo pipefail

MODE="${1:-all}"

FASTQ_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq"
THREADS="${THREADS:-8}"
TMP_PREFETCH="/tmp/sra_prefetch_$$"

# Activate the eclip environment (has sra-tools + pigz)
# env already activated by eclip_pipeline.slurm; prefetch/fasterq-dump/pigz on PATH

mkdir -p "${FASTQ_DIR}" "${TMP_PREFETCH}"
trap 'rm -rf "${TMP_PREFETCH}"' EXIT

# ── Helper ────────────────────────────────────────────────────────────────────
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
        --threads "${THREADS}" \
        --split-3

    # fasterq-dump names single-end reads <SRR>.fastq; rename to label
    if [[ -f "${FASTQ_DIR}/${srr}.fastq" ]]; then
        mv "${FASTQ_DIR}/${srr}.fastq" "${FASTQ_DIR}/${label}.fastq"
    elif [[ -f "${FASTQ_DIR}/${srr}_1.fastq" ]]; then
        # Unexpected paired-end: use read 1 (informative read for ENCODE4 eCLIP)
        mv "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${label}.fastq"
        rm -f "${FASTQ_DIR}/${srr}_2.fastq"
    fi

    pigz -p "${THREADS}" "${FASTQ_DIR}/${label}.fastq"
    rm -rf "${TMP_PREFETCH}/${srr}"
    echo "[DONE] ${label}"
}

# ── Human eCLIP IP samples (GEO: GSE276985) ──────────────────────────────────
download_human() {
    echo "=== Downloading human eCLIP IP samples ==="

    # TDP-43 | NGN2 (iPSC-derived neurons)
    download_one TDP43_NGN2_unstress_ip_rep1 SRR30639407  # SRX26061190
    download_one TDP43_NGN2_unstress_ip_rep2 SRR30639406  # SRX26061191
    download_one TDP43_NGN2_stress_ip_rep1   SRR30639405  # SRX26061192
    download_one TDP43_NGN2_stress_ip_rep2   SRR30639404  # SRX26061193

    # TDP-43 | NK (transdifferentiated neurons)
    download_one TDP43_NK_unstress_ip_rep1   SRR30639403  # SRX26061194
    download_one TDP43_NK_unstress_ip_rep2   SRR30639402  # SRX26061195

    # Caprin1 | NGN2
    download_one Caprin1_NGN2_unstress_ip_rep1 SRR30639393  # SRX26061204
    download_one Caprin1_NGN2_unstress_ip_rep2 SRR30639392  # SRX26061205
    download_one Caprin1_NGN2_stress_ip_rep1   SRR30639391  # SRX26061206
    download_one Caprin1_NGN2_stress_ip_rep2   SRR30639390  # SRX26061207

    # Caprin1 | NK
    download_one Caprin1_NK_unstress_ip_rep1   SRR30639389  # SRX26061208
    download_one Caprin1_NK_unstress_ip_rep2   SRR30639388  # SRX26061209

    # G3BP1 | NGN2
    download_one G3BP1_NGN2_unstress_ip_rep1   SRR30639387  # SRX26061210
    download_one G3BP1_NGN2_unstress_ip_rep2   SRR30639386  # SRX26061211
    download_one G3BP1_NGN2_stress_ip_rep1     SRR30639385  # SRX26061212
    download_one G3BP1_NGN2_stress_ip_rep2     SRR30639384  # SRX26061213

    # G3BP1 | NK (1 replicate only)
    download_one G3BP1_NK_unstress_ip_rep1     SRR30639383  # SRX26061214

    # G3BP1 | Human frontal cortex — mid-age donors (30–50 yr)
    download_one G3BP1_brain_midage_ip_rep1    SRR30639382  # SRX26061215
    download_one G3BP1_brain_midage_ip_rep2    SRR30639381  # SRX26061216
    download_one G3BP1_brain_midage_ip_rep3    SRR30639380  # SRX26061217

    # G3BP1 | Human frontal cortex — old-age donors (80–90 yr)
    download_one G3BP1_brain_oldage_ip_rep1    SRR30639379  # SRX26061218
    download_one G3BP1_brain_oldage_ip_rep2    SRR30639378  # SRX26061219
    download_one G3BP1_brain_oldage_ip_rep3    SRR30639377  # SRX26061220
}

# ── Mouse eCLIP IP samples (GEO: GSE276985) ───────────────────────────────────
download_mouse() {
    echo "=== Downloading mouse eCLIP IP samples ==="

    # TDP-43 | Mouse cerebellum — young-age (1.5 months), 4 replicates
    download_one TDP43_mouse_youngAge_ip_rep1 SRR30639401  # SRX26061196
    download_one TDP43_mouse_youngAge_ip_rep2 SRR30639400  # SRX26061197
    download_one TDP43_mouse_youngAge_ip_rep3 SRR30639399  # SRX26061198
    download_one TDP43_mouse_youngAge_ip_rep4 SRR30639398  # SRX26061199

    # TDP-43 | Mouse cerebellum — mid-age (6 months), 2 replicates
    download_one TDP43_mouse_midAge_ip_rep1   SRR30639397  # SRX26061200
    download_one TDP43_mouse_midAge_ip_rep2   SRR30639396  # SRX26061201

    # TDP-43 | Mouse cerebellum — old-age (22–25 months), 2 replicates
    download_one TDP43_mouse_oldAge_ip_rep1   SRR30639395  # SRX26061202
    download_one TDP43_mouse_oldAge_ip_rep2   SRR30639394  # SRX26061203
}

# ── Run ───────────────────────────────────────────────────────────────────────
case "${MODE}" in
    human) download_human ;;
    mouse) download_mouse ;;
    all)   download_human; download_mouse ;;
    *) echo "Usage: $0 [human|mouse|all]"; exit 1 ;;
esac

echo ""
echo "=== Download complete ==="
echo "FASTQ files are in: ${FASTQ_DIR}"
echo ""
echo "⚠  NEXT STEP — obtain SMInput (size-matched input) FASTQ files."
echo "   They are NOT available from GSE276985. Contact the authors or"
echo "   check the ENCODE portal for Yeo Lab eCLIP input samples."
echo "   Once you have input FASTQs, fill in Input_fastq paths in:"
echo "     config/manifest_human.csv"
echo "     config/manifest_mouse.csv"
echo "   Then run:  bash run_analysis.sh human"
