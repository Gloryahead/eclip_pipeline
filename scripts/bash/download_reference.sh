#!/usr/bin/env bash
# =============================================================================
# download_reference.sh — download all reference files needed by the pipeline
#
# Run this ONCE before submitting the pipeline. It fetches:
#   1. Cell Ranger ARC GRCh38-2024-A reference (used by Cell Ranger ATAC)
#   2. ENCODE hg38 blacklist (ENCFF356LFX)
#   3. Cell Ranger ATAC binary (you will be prompted to accept the 10x license)
#
# Usage:
#   bash /home/u11/maarowosegbe/scATAC-seq_pipeline/scripts/bash/download_reference.sh
# =============================================================================

set -euo pipefail

REF_DIR="/xdisk/haining/maarowosegbe/scATAC-seq/reference"
TOOLS_DIR="/xdisk/haining/maarowosegbe/scATAC-seq/tools"

mkdir -p "${REF_DIR}" "${TOOLS_DIR}"

# ── 1. Cell Ranger ARC GRCh38-2024-A reference ───────────────────────────────
# This is the ATAC+RNA joint reference but Cell Ranger ATAC can use it too.
# File size: ~15 GB. Takes 20-40 min to download.
CELLRANGER_REF_URL="https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2024-A.tar.gz"
CELLRANGER_REF_TARBALL="${REF_DIR}/refdata-cellranger-arc-GRCh38-2024-A.tar.gz"

if [[ -d "${REF_DIR}/refdata-cellranger-arc-GRCh38-2024-A" ]]; then
    echo "Cell Ranger reference already exists — skipping download."
else
    echo "Downloading Cell Ranger ARC GRCh38-2024-A reference (~15 GB)..."
    curl -L -o "${CELLRANGER_REF_TARBALL}" "${CELLRANGER_REF_URL}"
    echo "Extracting..."
    tar -xzf "${CELLRANGER_REF_TARBALL}" -C "${REF_DIR}"
    rm "${CELLRANGER_REF_TARBALL}"
    echo "Reference ready: ${REF_DIR}/refdata-cellranger-arc-GRCh38-2024-A"
fi

# ── 2. ENCODE hg38 blacklist (ENCFF356LFX) ────────────────────────────────────
# 910 regions of known artifact signal in the hg38 genome.
BLACKLIST_URL="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
BLACKLIST_BED="${REF_DIR}/ENCFF356LFX.bed"

if [[ -f "${BLACKLIST_BED}" ]]; then
    echo "Blacklist already exists — skipping."
else
    echo "Downloading ENCODE hg38 blacklist..."
    curl -L -o "${REF_DIR}/ENCFF356LFX.bed.gz" "${BLACKLIST_URL}"
    gunzip -k "${REF_DIR}/ENCFF356LFX.bed.gz"
    echo "Blacklist ready: ${BLACKLIST_BED}"
fi

# ── 3. Cell Ranger ATAC binary ────────────────────────────────────────────────
# Cell Ranger ATAC requires a free license from 10x Genomics.
# This section checks if you already have it installed.
CELLRANGER_ATAC_DIR="${TOOLS_DIR}/cellranger-atac-2.1.0"

if [[ -f "${CELLRANGER_ATAC_DIR}/cellranger-atac" ]]; then
    echo "Cell Ranger ATAC already installed at ${CELLRANGER_ATAC_DIR}"
else
    echo ""
    echo "=========================================================="
    echo "  Cell Ranger ATAC is NOT installed."
    echo ""
    echo "  To install:"
    echo "  1. Go to: https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads"
    echo "  2. Download cellranger-atac-2.1.0.tar.gz (requires free 10x account)"
    echo "  3. Upload it to the HPC: scp cellranger-atac-2.1.0.tar.gz netid@hpc.arizona.edu:${TOOLS_DIR}/"
    echo "  4. Extract it:"
    echo "       tar -xzf ${TOOLS_DIR}/cellranger-atac-2.1.0.tar.gz -C ${TOOLS_DIR}/"
    echo "  5. Verify: ${CELLRANGER_ATAC_DIR}/cellranger-atac --version"
    echo ""
    echo "  The path in nextflow.config (params.cellranger_atac_path) already"
    echo "  points to ${CELLRANGER_ATAC_DIR}/cellranger-atac"
    echo "=========================================================="
fi

# ── 4. AnnotationHub cache directory ─────────────────────────────────────────
# Pre-create this so all samples share the same Ensembl annotation cache,
# avoiding repeated downloads inside each SIGNAC_QC process.
CACHE_DIR="${REF_DIR}/annotationhub_cache"
mkdir -p "${CACHE_DIR}"
echo "AnnotationHub cache directory ready: ${CACHE_DIR}"

echo ""
echo "=== Reference setup complete ==="
echo "Paths:"
echo "  Cell Ranger reference : ${REF_DIR}/refdata-cellranger-arc-GRCh38-2024-A"
echo "  Blacklist BED         : ${BLACKLIST_BED}"
echo "  Cell Ranger ATAC bin  : ${CELLRANGER_ATAC_DIR}/cellranger-atac"
echo "  AnnotationHub cache   : ${CACHE_DIR}"
