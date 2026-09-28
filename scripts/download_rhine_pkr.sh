#!/usr/bin/env bash
# =============================================================================
# download_rhine_pkr.sh
# Downloads eCLIP FASTQs for Rhine et al. 2025 (GSE301408).
#
# Paper: "Mitochondrial double-stranded RNA triggers chronic stress in aged neurons"
# GEO:   GSE301408  BioProject: PRJNA1241817
#
# Targets: PKR, PNPT1, SUV3 (2 IP reps + 2 SMInput reps per condition)
# Cell types:
#   tdiff          = transdifferentiated neurons (aged model)
#   ipsc_untreated = iPSC-derived neurons, no treatment  (PAIRED, 51 bp R1)
#   ipsc_10uM      = iPSC-derived neurons + 10 µM AntimycinA  (SINGLE, 101 bp)
#   ipsc_5nM       = iPSC-derived neurons + 5 nM AntimycinA   (SINGLE, 101 bp)
#   brain_mid      = frontal cortex, mid-age donors  (PAIRED)
#   brain_old      = frontal cortex, old-age donors  (PAIRED)
#
# Layout notes:
#   PAIRED samples → we keep R1 only.  In the ENCODE4 / seCLIP protocol the
#   10 nt UMI is embedded at the 5′ end of R1 (the cDNA read), which is what
#   Skipper's se_preprocess.smk expects.  R2 is discarded.
#   SINGLE samples → downloaded directly, no splitting needed.
#
# Usage:
#   bash scripts/download_rhine_pkr.sh
#
# Prerequisites:
#   eclip conda environment activated (sra-tools, pigz)
# =============================================================================
set -euo pipefail

FASTQ_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq/rhine_pkr"
THREADS="${THREADS:-8}"
TMP_PREFETCH="/tmp/sra_prefetch_rhine_pkr_$$"

mkdir -p "${FASTQ_DIR}" "${TMP_PREFETCH}"
trap 'rm -rf "${TMP_PREFETCH}"' EXIT

# ---------------------------------------------------------------------------
# download_se  label srr
#   For SINGLE-end SRR: prefetch → fasterq-dump → compress → rename.
# ---------------------------------------------------------------------------
download_se() {
    local label="$1" srr="$2"
    local outfile="${FASTQ_DIR}/${label}.fastq.gz"

    if [[ -f "${outfile}" ]]; then
        echo "[SKIP] ${label}: already present."
        return 0
    fi

    echo "[DOWN-SE] ${label}  (${srr})"
    prefetch --max-size 50G -p "${srr}" -O "${TMP_PREFETCH}/"

    fasterq-dump "${TMP_PREFETCH}/${srr}/${srr}.sra" \
        --outdir "${FASTQ_DIR}" \
        --temp   "${TMP_PREFETCH}" \
        --threads "${THREADS}"

    mv "${FASTQ_DIR}/${srr}.fastq" "${FASTQ_DIR}/${label}.fastq"
    pigz -p "${THREADS}" "${FASTQ_DIR}/${label}.fastq"
    rm -rf "${TMP_PREFETCH}/${srr}"
    echo "[DONE] ${label}"
}

# ---------------------------------------------------------------------------
# download_pe_r1  label srr
#   For PAIRED-end SRR: prefetch → fasterq-dump --split-files → keep R1 only.
#   R2 is discarded — not needed for the single-end Skipper pipeline.
# ---------------------------------------------------------------------------
download_pe_r1() {
    local label="$1" srr="$2"
    local outfile="${FASTQ_DIR}/${label}.fastq.gz"

    if [[ -f "${outfile}" ]]; then
        echo "[SKIP] ${label}: already present."
        return 0
    fi

    echo "[DOWN-PE→R1] ${label}  (${srr})"
    prefetch --max-size 50G -p "${srr}" -O "${TMP_PREFETCH}/"

    fasterq-dump "${TMP_PREFETCH}/${srr}/${srr}.sra" \
        --outdir  "${FASTQ_DIR}" \
        --temp    "${TMP_PREFETCH}" \
        --threads "${THREADS}" \
        --split-files

    # Keep R1 (cDNA + 5′ UMI); discard R2
    mv "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${label}.fastq"
    rm -f "${FASTQ_DIR}/${srr}_2.fastq"

    pigz -p "${THREADS}" "${FASTQ_DIR}/${label}.fastq"
    rm -rf "${TMP_PREFETCH}/${srr}"
    echo "[DONE] ${label}"
}

echo "=== Downloading Rhine et al. 2025 eCLIP data (GSE301408) ==="
echo "Output: ${FASTQ_DIR}"
echo ""

# ── iPSC-derived neurons + PNPT1 (PAIRED, 51 bp R1) ─────────────────────────
echo "--- PNPT1 iPSC ---"
download_pe_r1 PNPT1_ipsc_ip_rep1    SRR34324053
download_pe_r1 PNPT1_ipsc_ip_rep2    SRR34324054
download_pe_r1 PNPT1_ipsc_input_rep1 SRR34324055
download_pe_r1 PNPT1_ipsc_input_rep2 SRR34324056

# ── iPSC-derived neurons + SUV3 (PAIRED, 51 bp R1) ───────────────────────────
echo "--- SUV3 iPSC ---"
download_pe_r1 SUV3_ipsc_ip_rep1    SRR34324057
download_pe_r1 SUV3_ipsc_ip_rep2    SRR34324058
download_pe_r1 SUV3_ipsc_input_rep1 SRR34324059
download_pe_r1 SUV3_ipsc_input_rep2 SRR34324060

# ── iPSC-derived neurons + PKR + 10 µM AntimycinA (SINGLE, 101 bp) ───────────
echo "--- PKR iPSC 10uM AntimycinA ---"
download_se PKR_ipsc_10uM_ip_rep1    SRR34324061
download_se PKR_ipsc_10uM_ip_rep2    SRR34324062
download_se PKR_ipsc_10uM_input_rep1 SRR34324063
download_se PKR_ipsc_10uM_input_rep2 SRR34324064

# ── iPSC-derived neurons + PKR + 5 nM AntimycinA (SINGLE, 101 bp) ────────────
echo "--- PKR iPSC 5nM AntimycinA ---"
download_se PKR_ipsc_5nM_ip_rep1    SRR34324065
download_se PKR_ipsc_5nM_ip_rep2    SRR34324066
download_se PKR_ipsc_5nM_input_rep1 SRR34324067
download_se PKR_ipsc_5nM_input_rep2 SRR34324068

# ── iPSC-derived neurons + PKR untreated (PAIRED, 51 bp R1) ─────────────────
echo "--- PKR iPSC untreated ---"
download_pe_r1 PKR_ipsc_untr_ip_rep1    SRR34324069
download_pe_r1 PKR_ipsc_untr_ip_rep2    SRR34324070
download_pe_r1 PKR_ipsc_untr_input_rep1 SRR34324071
download_pe_r1 PKR_ipsc_untr_input_rep2 SRR34324072

# ── Transdifferentiated neurons (Tdiff) + PNPT1 (PAIRED, 51 bp R1) ───────────
echo "--- PNPT1 Tdiff ---"
download_pe_r1 PNPT1_tdiff_ip_rep1    SRR34324073
download_pe_r1 PNPT1_tdiff_ip_rep2    SRR34324074
download_pe_r1 PNPT1_tdiff_input_rep1 SRR34324075
download_pe_r1 PNPT1_tdiff_input_rep2 SRR34324076

# ── Transdifferentiated neurons (Tdiff) + SUV3 (PAIRED, 51 bp R1) ────────────
echo "--- SUV3 Tdiff ---"
download_pe_r1 SUV3_tdiff_ip_rep1    SRR34324077
download_pe_r1 SUV3_tdiff_ip_rep2    SRR34324078
download_pe_r1 SUV3_tdiff_input_rep1 SRR34324079
download_pe_r1 SUV3_tdiff_input_rep2 SRR34324080

# ── Transdifferentiated neurons (Tdiff) + PKR (PAIRED, 101 bp R1) ────────────
# These are the central samples — PKR activation in aged neurons.
echo "--- PKR Tdiff (primary condition) ---"
download_pe_r1 PKR_tdiff_ip_rep1    SRR34324081
download_pe_r1 PKR_tdiff_ip_rep2    SRR34324082
download_pe_r1 PKR_tdiff_input_rep1 SRR34324083
download_pe_r1 PKR_tdiff_input_rep2 SRR34324084

# ── Frontal cortex — mid-age brain donors (PAIRED) ───────────────────────────
# TODO: verify which SRRs are mid-age vs old-age from individual GSM records.
# Best current guess from GEO ordering: SRR34324123-126 = mid, 127-130 = old.
echo "--- PKR brain mid-age ---"
download_pe_r1 PKR_brain_mid_ip_rep1    SRR34324123
download_pe_r1 PKR_brain_mid_input_rep1 SRR34324124
download_pe_r1 PKR_brain_mid_ip_rep2    SRR34324125
download_pe_r1 PKR_brain_mid_input_rep2 SRR34324126

# ── Frontal cortex — old-age brain donors (PAIRED) ───────────────────────────
echo "--- PKR brain old-age ---"
download_pe_r1 PKR_brain_old_ip_rep1    SRR34324127
download_pe_r1 PKR_brain_old_input_rep1 SRR34324128
download_pe_r1 PKR_brain_old_ip_rep2    SRR34324129
download_pe_r1 PKR_brain_old_input_rep2 SRR34324130

echo ""
echo "=== Download complete ==="
echo "Files in: ${FASTQ_DIR}"
echo ""
echo "NEXT: sbatch rhine_pkr.slurm"
