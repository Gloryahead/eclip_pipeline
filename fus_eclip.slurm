#!/bin/bash
# =============================================================================
# fus_eclip.slurm
# Skipper eCLIP pipeline — FUS pull-down (user's own samples, 3 conditions)
#
# Conditions  : FUS_healthy (010), FUS_disease (004), FUS_remedy (002)
# Protocol    : seCLIP — ENCODE4, 10 nt UMI at 5' of R1
# Input FASTQs: already on Puma at
#               /xdisk/haining/maarowosegbe/Lei_Wang_eCLIP/9040-0003-10012026_153425/
# Results will appear in:
#   /xdisk/haining/maarowosegbe/eclip_pipeline/results/fus_eclip/
#
# Usage:
#   cd ~/eclip_pipeline && git pull
#   sbatch fus_eclip.slurm
# =============================================================================

#SBATCH --job-name=fus_eclip
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
#    default 40M overflows at high-depth rRNA positions.
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

echo "--- Step 1: Running FUS eCLIP Skipper pipeline (3 conditions) ---"
bash "${SCRIPT_BASE}/run_analysis.sh" fus_eclip

# Step 2: Negate minus-strand bigwigs for IGV display (peaks face downward)
echo "--- Step 2: Negating minus-strand bigwigs ---"
bash "${SCRIPT_BASE}/scripts/negate_minus_bigwigs.sh" \
    /xdisk/haining/maarowosegbe/eclip_pipeline/results/fus_eclip \
    /xdisk/haining/maarowosegbe/eclip_pipeline/references/hg38/star_index/chrNameLength.txt \
    8

echo "========== Job finished at: $(date) =========="
