#!/bin/bash
# =============================================================================
# g3bp1_test.slurm
# Full pipeline test run using G3BP1 seCLIP data from He et al. 2021 (GSE168943).
# This is a standalone job — independent of the Rhine et al. main pipeline.
#
# What it does:
#   1. Downloads 4 FASTQ files (~6 GB) from NCBI SRA into a separate subdirectory
#   2. Runs Skipper eCLIP pipeline on the G3BP1 Ctrl condition (2 IP + 2 SMInput)
#
# Usage:
#   cd ~/eclip_pipeline && git pull
#   sbatch g3bp1_test.slurm
#
# Results will appear in:
#   /xdisk/haining/maarowosegbe/eclip_pipeline/results/g3bp1_test/
# =============================================================================

#SBATCH --job-name=g3bp1_test
#SBATCH --output=/xdisk/haining/maarowosegbe/eclip_pipeline/logs/%x-%j.out
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=48:00:00

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu; source ~/.bashrc; set -eu
mamba-haining
micromamba activate eclip
export PATH="${CONDA_PREFIX}/bin:${PATH}"

SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"

# Step 1: Download FASTQ files (skips any already present)
echo "--- Step 1: Downloading G3BP1 test FASTQs ---"
bash "${SCRIPT_BASE}/scripts/download_g3bp1_test.sh"

# Step 2: Run Skipper pipeline
echo "--- Step 2: Running Skipper pipeline ---"
bash "${SCRIPT_BASE}/run_analysis.sh" g3bp1_test

echo "========== Job finished at: $(date) =========="
