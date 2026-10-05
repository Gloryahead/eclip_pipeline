#!/bin/bash
# =============================================================================
# rhine_pkr.slurm
# Rhine et al. 2025 (GSE301408) eCLIP-seq pipeline — PKR, PNPT1, SUV3.
#
# Paper: "Mitochondrial double-stranded RNA triggers chronic stress in aged neurons"
#
# What it does:
#   1. Downloads all FASTQ files (~300+ GB total) from NCBI SRA (skips existing)
#   2. Runs Skipper eCLIP pipeline for 10 conditions:
#        PKR_tdiff, PKR_ipsc_untr, PKR_ipsc_10uM, PKR_ipsc_5nM,
#        PNPT1_tdiff, PNPT1_ipsc, SUV3_tdiff, SUV3_ipsc,
#        PKR_brain_mid, PKR_brain_old
#
# Usage:
#   cd ~/eclip_pipeline && git pull
#   sbatch rhine_pkr.slurm
#
# Results will appear in:
#   /xdisk/haining/maarowosegbe/eclip_pipeline/results/rhine_pkr/
# =============================================================================

#SBATCH --job-name=rhine_pkr
#SBATCH --output=/xdisk/haining/maarowosegbe/eclip_pipeline/logs/%x-%j.out
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=96:00:00

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu; source ~/.bashrc; set -eu
mamba-haining
micromamba activate eclip
export PATH="${CONDA_PREFIX}/bin:/groups/haining/maarowosegbe/micromamba/envs/snakemake_env/bin:${PATH}"

SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"

# Patch Skipper submodule — two fixes, both idempotent:
#
# 1. -Xss512m: fixes StackOverflow in recursive Directional DFS at high-depth
#    positions. Regex matches any existing -Xss value.
#
# 2. -Xms4000M: fixes OOM kill (exit 137). Skipper sets Xms=Xmx which pre-commits
#    the full heap at startup, breaching the SLURM cgroup limit. Starting at 4 GB
#    and growing to Xmx on demand keeps RSS well within limits.
sed -i \
  's/-Xss[0-9]*[mMgGkK]/-Xss512m/g;
   s/-Xms{resources\.java_mem_mb}M -Xmx/-Xms4000M -Xmx/g' \
  "${SCRIPT_BASE}/skipper/rules/se_preprocess.smk"

# Step 1: Download FASTQ files (skips any already present)
echo "--- Step 1: Downloading Rhine et al. 2025 FASTQs (GSE301408) ---"
bash "${SCRIPT_BASE}/scripts/download_rhine_pkr.sh"

# Step 2: Run Skipper pipeline
echo "--- Step 2: Running Skipper pipeline ---"
bash "${SCRIPT_BASE}/run_analysis.sh" rhine_pkr

echo "========== Job finished at: $(date) =========="
