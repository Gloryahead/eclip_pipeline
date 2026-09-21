#!/usr/bin/env bash
# =============================================================================
# 01_download_references.sh
# Downloads and prepares all reference files needed by the Skipper eCLIP
# pipeline for both human (hg38/GRCh38) and mouse (mm10/GRCm38) genomes.
#
# What it does:
#   1. Downloads GENCODE genome FASTA and GFF3 annotation
#   2. Filters the GFF3 with Skipper's transcript quality filter
#   3. Builds STAR genome indices for each genome
#   4. Downloads ENCODE eCLIP blacklist regions (optional but recommended)
#
# Reference versions used (match paper's analysis timeframe, June 2025 pub):
#   Human : GENCODE v47 (GRCh38 / hg38)
#   Mouse : GENCODE M35 (GRCm38 / mm10)
#
# CPU/memory: STAR indexing requires ~32 GB RAM and 8 CPUs. Submit as a
# SLURM job on the HPC, not as an interactive session.
#
# Usage:
#   sbatch setup/01_download_references.sh
#   OR on an interactive high-memory node:
#   bash setup/01_download_references.sh
# =============================================================================
#SBATCH --job-name=eclip_build_refs
#SBATCH --partition=standard
#SBATCH --account=haining
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=06:00:00
#SBATCH --output=logs/build_refs_%j.out
#SBATCH --error=logs/build_refs_%j.err

set -euo pipefail

# ── User-configurable ──────────────────────────────────────────────────────
BASE_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline"
REF_DIR="${BASE_DIR}/references"
SKIPPER_DIR="/home/u11/maarowosegbe/eclip_pipeline/skipper/gff_utils"
THREADS=16

# env already activated by eclip_pipeline.slurm; tools are on PATH

mkdir -p "${REF_DIR}"/{hg38,mm10} logs

# ── Human (hg38 / GRCh38, GENCODE v47) ────────────────────────────────────
HG38_DIR="${REF_DIR}/hg38"
GENCODE_HUMAN_VER="47"
GENCODE_HUMAN_FASTA="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_HUMAN_VER}/GRCh38.primary_assembly.genome.fa.gz"
GENCODE_HUMAN_GFF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_HUMAN_VER}/gencode.v${GENCODE_HUMAN_VER}.annotation.gff3.gz"

echo "=== [hg38 1/4] Downloading human genome FASTA (GENCODE v${GENCODE_HUMAN_VER}) ==="
if [[ ! -f "${HG38_DIR}/GRCh38.primary_assembly.genome.fa.gz" ]]; then
    wget -c -P "${HG38_DIR}" "${GENCODE_HUMAN_FASTA}"
else
    echo "Human FASTA already present, skipping."
fi

echo "=== [hg38 2/4] Downloading human annotation GFF3 ==="
if [[ ! -f "${HG38_DIR}/gencode.v${GENCODE_HUMAN_VER}.annotation.gff3.gz" ]]; then
    wget -c -P "${HG38_DIR}" "${GENCODE_HUMAN_GFF}"
else
    echo "Human GFF3 already present, skipping."
fi

echo "=== [hg38 3/4] Filtering GFF3 for transcript quality (Skipper requirement) ==="
FILTERED_GFF="${HG38_DIR}/gencode.v${GENCODE_HUMAN_VER}.filtered.gff3.gz"
if [[ ! -f "${FILTERED_GFF}" ]]; then
    Rscript \
        "${SKIPPER_DIR}/gff_transcript_quality_filter.R" \
        gencode \
        "${HG38_DIR}/gencode.v${GENCODE_HUMAN_VER}.annotation.gff3.gz" \
        "${FILTERED_GFF}"
else
    echo "Filtered GFF3 already present, skipping."
fi

echo "=== [hg38 4/4] Building STAR genome index ==="
STAR_INDEX_HG38="${HG38_DIR}/star_index"
if [[ ! -d "${STAR_INDEX_HG38}/Genome" ]]; then
    mkdir -p "${STAR_INDEX_HG38}"
    HG38_FASTA_TMP="${HG38_DIR}/GRCh38.primary_assembly.genome.fa"
    echo "Decompressing hg38 FASTA for STAR (STAR requires uncompressed input)..."
    zcat "${HG38_DIR}/GRCh38.primary_assembly.genome.fa.gz" > "${HG38_FASTA_TMP}"
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX_HG38}" \
        --outFileNamePrefix "${STAR_INDEX_HG38}/" \
        --genomeFastaFiles "${HG38_FASTA_TMP}" \
        --sjdbGTFfile "${FILTERED_GFF}" \
        --sjdbGTFfeatureExon CDS \
        --runThreadN "${THREADS}" \
        --genomeSAindexNbases 14 \
        --limitGenomeGenerateRAM 60000000000
    rm -f "${HG38_FASTA_TMP}"
else
    echo "STAR index for hg38 already exists, skipping."
fi

# ── Mouse (mm10 / GRCm38, GENCODE M35) ────────────────────────────────────
MM10_DIR="${REF_DIR}/mm10"
GENCODE_MOUSE_VER="M35"
GENCODE_MOUSE_FASTA="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_MOUSE_VER}/GRCm38.primary_assembly.genome.fa.gz"
GENCODE_MOUSE_GFF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_MOUSE_VER}/gencode.${GENCODE_MOUSE_VER}.annotation.gff3.gz"

echo "=== [mm10 1/4] Downloading mouse genome FASTA (GENCODE ${GENCODE_MOUSE_VER}) ==="
if [[ ! -f "${MM10_DIR}/GRCm38.primary_assembly.genome.fa.gz" ]]; then
    wget -c -P "${MM10_DIR}" "${GENCODE_MOUSE_FASTA}"
else
    echo "Mouse FASTA already present, skipping."
fi

echo "=== [mm10 2/4] Downloading mouse annotation GFF3 ==="
if [[ ! -f "${MM10_DIR}/gencode.${GENCODE_MOUSE_VER}.annotation.gff3.gz" ]]; then
    wget -c -P "${MM10_DIR}" "${GENCODE_MOUSE_GFF}"
else
    echo "Mouse GFF3 already present, skipping."
fi

echo "=== [mm10 3/4] Filtering mouse GFF3 ==="
FILTERED_GFF_MM10="${MM10_DIR}/gencode.${GENCODE_MOUSE_VER}.filtered.gff3.gz"
if [[ ! -f "${FILTERED_GFF_MM10}" ]]; then
    Rscript \
        "${SKIPPER_DIR}/gff_transcript_quality_filter.R" \
        gencode \
        "${MM10_DIR}/gencode.${GENCODE_MOUSE_VER}.annotation.gff3.gz" \
        "${FILTERED_GFF_MM10}"
else
    echo "Filtered mouse GFF3 already present, skipping."
fi

echo "=== [mm10 4/4] Building STAR genome index ==="
STAR_INDEX_MM10="${MM10_DIR}/star_index"
if [[ ! -d "${STAR_INDEX_MM10}/Genome" ]]; then
    mkdir -p "${STAR_INDEX_MM10}"
    MM10_FASTA_TMP="${MM10_DIR}/GRCm38.primary_assembly.genome.fa"
    echo "Decompressing mm10 FASTA for STAR (STAR requires uncompressed input)..."
    zcat "${MM10_DIR}/GRCm38.primary_assembly.genome.fa.gz" > "${MM10_FASTA_TMP}"
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX_MM10}" \
        --outFileNamePrefix "${STAR_INDEX_MM10}/" \
        --genomeFastaFiles "${MM10_FASTA_TMP}" \
        --sjdbGTFfile "${FILTERED_GFF_MM10}" \
        --sjdbGTFfeatureExon CDS \
        --runThreadN "${THREADS}" \
        --genomeSAindexNbases 14 \
        --limitGenomeGenerateRAM 60000000000
    rm -f "${MM10_FASTA_TMP}"
else
    echo "STAR index for mm10 already exists, skipping."
fi

# ── ENCODE eCLIP blacklist (optional) ─────────────────────────────────────
echo "=== [opt] Downloading ENCODE eCLIP blacklist regions ==="
BLACKLIST="${REF_DIR}/encode_eclip_blacklist.bed"
if [[ ! -f "${BLACKLIST}" ]]; then
    wget -c -O "${BLACKLIST}.gz" \
        "https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz" \
        && gunzip "${BLACKLIST}.gz" || true
else
    echo "Blacklist already present, skipping."
fi

echo ""
echo "=== Reference preparation complete ==="
echo "Human refs : ${HG38_DIR}"
echo "Mouse refs : ${MM10_DIR}"
echo "Next step  : fill in config/manifest_human.csv and config/manifest_mouse.csv"
echo "Then run   : bash run_analysis.sh human  OR  bash run_analysis.sh mouse"
