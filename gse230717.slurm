#!/bin/bash
# =============================================================================
# gse230717.slurm
# Gosztyla et al. 2023 (GSE230717) eCLIP-seq pipeline — G3BP1, CAPRIN1, SNRNP200.
#
# Paper: "Large-scale map of RNA binding protein interactomes across the mRNA life-cycle"
#
# What it does:
#   1. Downloads all 24 FASTQ files from NCBI SRA (skips existing)
#   2. Runs Skipper eCLIP pipeline for 6 experiments:
#        G3BP1, CAPRIN1,
#        SNRNP200_wc, SNRNP200_cy,
#        SNRNP200_wcs, SNRNP200_cys
#
# Usage:
#   cd ~/eclip_pipeline && git pull
#   sbatch gse230717.slurm
#
# Results will appear in:
#   /xdisk/haining/maarowosegbe/eclip_pipeline/results/gse230717/
# =============================================================================

#SBATCH --job-name=gse230717
#SBATCH --output=/xdisk/haining/maarowosegbe/eclip_pipeline/logs/%x-%j.out
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=96:00:00

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu
source ~/.bashrc
mamba-haining
micromamba activate eclip
set -eu
export PATH="${CONDA_PREFIX}/bin:/groups/haining/maarowosegbe/micromamba/envs/snakemake_env/bin:${PATH}"

SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"

# Patch Skipper submodule — idempotent; safe to re-run.
#
# 1. -Xss512m: fixes StackOverflow in recursive Directional DFS at high-depth positions.
# 2. -Xms4000M: decouples initial heap from max heap to avoid SLURM cgroup OOM kill.
python3 "${SCRIPT_BASE}/scripts/patch_se_preprocess_rule.py" \
  "${SCRIPT_BASE}/skipper/rules/se_preprocess.smk"

# 3. Stream annotation file + direct bedtools call to fix annotate_finemap OOM.
python3 "${SCRIPT_BASE}/scripts/patch_annotate_finemap_rule.py" \
  "${SCRIPT_BASE}/skipper/rules/finemap.smk" \
  "${SCRIPT_BASE}/scripts/annotate_finemapped_regions.py"

# 4. HOMER process-substitution fix: replace <(...) with mktemp files.
python3 "${SCRIPT_BASE}/scripts/patch_homer_rule.py" \
  "${SCRIPT_BASE}/skipper/rules/analysis.smk"

# Step 1: Download FASTQ files (skips any already present)
echo "--- Step 1: Downloading Gosztyla et al. 2023 FASTQs (GSE230717) ---"
bash "${SCRIPT_BASE}/scripts/download_gse230717.sh"

# Step 2: Run Skipper pipeline
echo "--- Step 2: Running Skipper pipeline ---"
bash "${SCRIPT_BASE}/run_analysis.sh" gse230717

# Step 3: Negate minus-strand bigwigs for IGV display (peaks face downward)
echo "--- Step 3: Negating minus-strand bigwigs ---"
bash "${SCRIPT_BASE}/scripts/negate_minus_bigwigs.sh" \
    /xdisk/haining/maarowosegbe/eclip_pipeline/results/gse230717 \
    /xdisk/haining/maarowosegbe/eclip_pipeline/references/hg38/star_index/chrNameLength.txt \
    8

echo "========== Job finished at: $(date) =========="
