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
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu
source ~/.bashrc
mamba-haining
micromamba activate eclip
set -eu
export PATH="${CONDA_PREFIX}/bin:/groups/haining/maarowosegbe/micromamba/envs/snakemake_env/bin:${PATH}"

SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"

# Patch Skipper submodule — two fixes, both idempotent:
#
# 1. -Xss: increase Java stack from 40M to 512m.
#    Directional dedup uses recursive DFS over UMI edit-distance graph; the
#    default 40M overflows at high-depth rRNA positions in G3BP1 data.
#    2G was tried but exceeds the JVM signed-32-bit stack parser limit.
#    Regex matches any existing -Xss value (40M, 2G, etc.).
#
# 2. -Xms: decouple initial heap from max heap.
#    Skipper sets Xms=Xmx=32/48/72 GB. That pre-commits the full heap at
#    startup; add thread stacks and the cgroup limit is breached → OOM kill
#    (exit 137) even though actual heap usage is much lower. Fix: start at
#    4 GB and grow on demand up to the retry-scaled Xmx.
sed -i \
  's/-Xss[0-9]*[mMgGkK]/-Xss512m/g;
   s/-Xms{resources\.java_mem_mb}M -Xmx/-Xms4000M -Xmx/g' \
  "${SCRIPT_BASE}/skipper/rules/se_preprocess.smk"

sed -i \
  "s/finemapped_sites_df = finemapped_sites\.to_dataframe()/finemapped_sites_df = finemapped_sites.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'input', 'clip', 'cdf'])/" \
  "${SCRIPT_BASE}/skipper/tools/annotate_finemapped_regions.py"

# Step 1: Download FASTQ files (skips any already present)
echo "--- Step 1: Downloading G3BP1 test FASTQs ---"
bash "${SCRIPT_BASE}/scripts/download_g3bp1_test.sh"

# Step 2: Run Skipper pipeline
echo "--- Step 2: Running Skipper pipeline ---"
bash "${SCRIPT_BASE}/run_analysis.sh" g3bp1_test

# Step 3: Negate minus-strand bigwigs for IGV display (peaks face downward)
echo "--- Step 3: Negating minus-strand bigwigs ---"
bash "${SCRIPT_BASE}/scripts/negate_minus_bigwigs.sh" \
    /xdisk/haining/maarowosegbe/eclip_pipeline/results/g3bp1_test \
    /xdisk/haining/maarowosegbe/eclip_pipeline/references/hg38/star_index/chrNameLength.txt \
    8

echo "========== Job finished at: $(date) =========="
