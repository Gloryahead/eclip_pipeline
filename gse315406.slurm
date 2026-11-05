#!/bin/bash
# =============================================================================
# gse315406.slurm
# ENCORE (GSE315406) eCLIP-seq pipeline — 57 RBPs in K562 and HepG2 cells.
#
# Paper: "Comprehensive RNA-binding protein analyses using enhanced CLIP (ENCORE)"
#
# What it does:
#   0. Builds config/manifest_gse315406.csv from NCBI metadata (Step 0)
#   1. Generates download_gse315406.sh from the manifest SRR accessions
#   2. Downloads all FASTQs from NCBI SRA (skips existing)
#   3. Runs Skipper eCLIP pipeline for all 57 experiments
#
# Usage:
#   cd ~/eclip_pipeline && git pull
#   sbatch gse315406.slurm
#
# Results will appear in:
#   /xdisk/haining/maarowosegbe/eclip_pipeline/results/gse315406/
#
# NOTE: This job downloads ~300+ GB of FASTQ data. The 120-hour walltime is set
# accordingly. Monitor disk usage on xdisk before submitting.
# =============================================================================

#SBATCH --job-name=gse315406
#SBATCH --output=/xdisk/haining/maarowosegbe/eclip_pipeline/logs/%x-%j.out
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=120:00:00

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu
source ~/.bashrc
mamba-haining
micromamba activate eclip
set -eu
export PATH="${CONDA_PREFIX}/bin:/groups/haining/maarowosegbe/micromamba/envs/snakemake_env/bin:${PATH}"

SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"
FASTQ_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq/gse315406"
THREADS=8

# Patch Skipper submodule — idempotent; safe to re-run.
python3 "${SCRIPT_BASE}/scripts/patch_se_preprocess_rule.py" \
  "${SCRIPT_BASE}/skipper/rules/se_preprocess.smk"
python3 "${SCRIPT_BASE}/scripts/patch_annotate_finemap_rule.py" \
  "${SCRIPT_BASE}/skipper/rules/finemap.smk" \
  "${SCRIPT_BASE}/scripts/annotate_finemapped_regions.py"
python3 "${SCRIPT_BASE}/scripts/patch_homer_rule.py" \
  "${SCRIPT_BASE}/skipper/rules/analysis.smk"

# ── Step 0: Build manifest from NCBI metadata ───────────────────────────────
echo "--- Step 0: Building manifest for GSE315406 ---"
python3 "${SCRIPT_BASE}/scripts/build_manifest_gse315406.py" \
    > /tmp/gse315406_srr_list.txt 2>&1
cat /tmp/gse315406_srr_list.txt

# Extract download commands and generate the download script
grep "^download_se " /tmp/gse315406_srr_list.txt \
    > /tmp/gse315406_downloads.txt || true

echo "    Manifest written to ${SCRIPT_BASE}/config/manifest_gse315406.csv"

# ── Step 1: Download FASTQs ──────────────────────────────────────────────────
echo "--- Step 1: Downloading ENCORE FASTQs (GSE315406) ---"

TMP_PREFETCH="/tmp/sra_prefetch_gse315406_$$"
mkdir -p "${FASTQ_DIR}" "${TMP_PREFETCH}"

download_se() {
    local label="$1" srr="$2"
    local outfile="${FASTQ_DIR}/${label}.fastq.gz"
    if [[ -f "${outfile}" ]]; then
        echo "[SKIP] ${label}"
        return 0
    fi
    echo "[DOWN] ${label}  (${srr})"
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

# Run each download_se line extracted from the manifest builder output
while IFS=' ' read -r cmd label srr; do
    [[ "${cmd}" == "download_se" ]] && download_se "${label}" "${srr}"
done < /tmp/gse315406_downloads.txt

rm -rf "${TMP_PREFETCH}"

# ── Step 2: Run Skipper pipeline ─────────────────────────────────────────────
echo "--- Step 2: Running Skipper pipeline ---"
bash "${SCRIPT_BASE}/run_analysis.sh" gse315406

# ── Step 3: Negate minus-strand bigwigs ──────────────────────────────────────
echo "--- Step 3: Negating minus-strand bigwigs ---"
bash "${SCRIPT_BASE}/scripts/negate_minus_bigwigs.sh" \
    /xdisk/haining/maarowosegbe/eclip_pipeline/results/gse315406 \
    /xdisk/haining/maarowosegbe/eclip_pipeline/references/hg38/star_index/chrNameLength.txt \
    8

echo "========== Job finished at: $(date) =========="
