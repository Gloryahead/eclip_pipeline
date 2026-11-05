#!/usr/bin/env bash
# =============================================================================
# download_gse230717.sh
# Downloads eCLIP FASTQs for Gosztyla et al. 2023 (GSE230717).
#
# Paper: "Large-scale map of RNA binding protein interactomes across the
#         mRNA life-cycle" (Mol Cell, 2023)
# GEO:   GSE230717
#
# Targets: G3BP1, CAPRIN1, SNRNP200
# Cell type: human epithelial cells
# Conditions (SNRNP200 only):
#   wc  = whole cell, unstressed
#   cy  = cytoplasmic, unstressed
#   wcs = whole cell + stress (sodium arsenite)
#   cys = cytoplasmic + stress (sodium arsenite)
#
# Protocol: ENCODE4 / seCLIP — 10 nt UMI at 5′ of R1 (informative read).
#   Samples are paired-end; R1 only is kept for the Skipper single-end pipeline.
#
# SRR accessions verified from NCBI SRA (SRX20116040–SRX20116063):
#   G3BP1 IP:          SRR24321513 (rep1), SRR24321512 (rep2)
#   G3BP1 SMInput:     SRR24321525 (rep1), SRR24321524 (rep2)
#   CAPRIN1 IP:        SRR24321511 (rep1), SRR24321510 (rep2)
#   CAPRIN1 SMInput:   SRR24321523 (rep1), SRR24321522 (rep2)
#   SNRNP200 WC IP:    SRR24321517 (rep1), SRR24321516 (rep2)
#   SNRNP200 WC SMIn:  SRR24321529 (rep1), SRR24321528 (rep2)
#   SNRNP200 CY IP:    SRR24321515 (rep1), SRR24321514 (rep2)
#   SNRNP200 CY SMIn:  SRR24321527 (rep1), SRR24321526 (rep2)
#   SNRNP200 WCS IP:   SRR24321509 (rep1), SRR24321508 (rep2)
#   SNRNP200 WCS SMIn: SRR24321521 (rep1), SRR24321520 (rep2)
#   SNRNP200 CYS IP:   SRR24321507 (rep1), SRR24321506 (rep2)
#   SNRNP200 CYS SMIn: SRR24321519 (rep1), SRR24321518 (rep2)
#
# Usage:
#   bash scripts/download_gse230717.sh
#
# Prerequisites:
#   eclip conda environment activated (sra-tools, pigz)
# =============================================================================
set -euo pipefail

FASTQ_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq/gse230717"
THREADS="${THREADS:-8}"
TMP_PREFETCH="/tmp/sra_prefetch_gse230717_$$"

mkdir -p "${FASTQ_DIR}" "${TMP_PREFETCH}"
trap 'rm -rf "${TMP_PREFETCH}"' EXIT

# ---------------------------------------------------------------------------
# download_pe_r1  label srr
#   PAIRED-end: prefetch → fasterq-dump --split-files → keep R1 only.
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

echo "=== Downloading Gosztyla et al. 2023 eCLIP data (GSE230717) ==="
echo "Output: ${FASTQ_DIR}"
echo ""

# ── G3BP1 ────────────────────────────────────────────────────────────────────
echo "--- G3BP1 ---"
download_pe_r1 G3BP1_ip_rep1    SRR24321513
download_pe_r1 G3BP1_ip_rep2    SRR24321512
download_pe_r1 G3BP1_input_rep1 SRR24321525
download_pe_r1 G3BP1_input_rep2 SRR24321524

# ── CAPRIN1 ───────────────────────────────────────────────────────────────────
echo "--- CAPRIN1 ---"
download_pe_r1 CAPRIN1_ip_rep1    SRR24321511
download_pe_r1 CAPRIN1_ip_rep2    SRR24321510
download_pe_r1 CAPRIN1_input_rep1 SRR24321523
download_pe_r1 CAPRIN1_input_rep2 SRR24321522

# ── SNRNP200 whole cell, unstressed ──────────────────────────────────────────
echo "--- SNRNP200 whole cell (unstressed) ---"
download_pe_r1 SNRNP200_wc_ip_rep1    SRR24321517
download_pe_r1 SNRNP200_wc_ip_rep2    SRR24321516
download_pe_r1 SNRNP200_wc_input_rep1 SRR24321529
download_pe_r1 SNRNP200_wc_input_rep2 SRR24321528

# ── SNRNP200 cytoplasmic, unstressed ─────────────────────────────────────────
echo "--- SNRNP200 cytoplasmic (unstressed) ---"
download_pe_r1 SNRNP200_cy_ip_rep1    SRR24321515
download_pe_r1 SNRNP200_cy_ip_rep2    SRR24321514
download_pe_r1 SNRNP200_cy_input_rep1 SRR24321527
download_pe_r1 SNRNP200_cy_input_rep2 SRR24321526

# ── SNRNP200 whole cell, stress ───────────────────────────────────────────────
echo "--- SNRNP200 whole cell (stress) ---"
download_pe_r1 SNRNP200_wcs_ip_rep1    SRR24321509
download_pe_r1 SNRNP200_wcs_ip_rep2    SRR24321508
download_pe_r1 SNRNP200_wcs_input_rep1 SRR24321521
download_pe_r1 SNRNP200_wcs_input_rep2 SRR24321520

# ── SNRNP200 cytoplasmic, stress ──────────────────────────────────────────────
echo "--- SNRNP200 cytoplasmic (stress) ---"
download_pe_r1 SNRNP200_cys_ip_rep1    SRR24321507
download_pe_r1 SNRNP200_cys_ip_rep2    SRR24321506
download_pe_r1 SNRNP200_cys_input_rep1 SRR24321519
download_pe_r1 SNRNP200_cys_input_rep2 SRR24321518

echo ""
echo "=== Download complete ==="
echo "Files in: ${FASTQ_DIR}"
echo ""
echo "NEXT: sbatch gse230717.slurm"
