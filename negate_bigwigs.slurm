#!/bin/bash
# =============================================================================
# negate_bigwigs.slurm
# Post-processing: negate all minus-strand bigwigs so peaks face downward
# in IGV, matching the display convention used in eCLIP publications.
#
# Runs scripts/negate_minus_bigwigs.sh for one or all experiments.
# Skips any *_neg.bw files that already exist, so re-runs are safe.
#
# Usage — all experiments (auto-discovers results with minus bigwigs):
#   sbatch negate_bigwigs.slurm
#
# Usage — one specific experiment:
#   sbatch --export=ALL,EXPERIMENT=rhine_pkr negate_bigwigs.slurm
#   sbatch --export=ALL,EXPERIMENT=fus_eclip  negate_bigwigs.slurm
#   sbatch --export=ALL,EXPERIMENT=g3bp1_test negate_bigwigs.slurm
#
# Output files land alongside the originals, e.g.:
#   results/rhine_pkr/bigwig/PKR_tdiff_IP_1.unscaled.minus_neg.bw
# =============================================================================

#SBATCH --job-name=negate_bigwigs
#SBATCH --output=/xdisk/haining/maarowosegbe/eclip_pipeline/logs/%x-%j.out
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu
source ~/.bashrc
mamba-haining
micromamba activate eclip
set -eu
export PATH="${CONDA_PREFIX}/bin:/groups/haining/maarowosegbe/micromamba/envs/snakemake_env/bin:${PATH}"

SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"
RESULTS_BASE="/xdisk/haining/maarowosegbe/eclip_pipeline/results"
CHROM_SIZES="/xdisk/haining/maarowosegbe/eclip_pipeline/references/hg38/star_index/chrNameLength.txt"
THREADS=8

# ── Resolve which experiments to process ──────────────────────────────────
# If EXPERIMENT is set (via --export=ALL,EXPERIMENT=rhine_pkr), process only
# that one. Otherwise auto-discover every results subdirectory that already
# contains minus-strand bigwig files.
if [[ -n "${EXPERIMENT:-}" ]]; then
    EXPERIMENTS=("${EXPERIMENT}")
    echo "Target experiment: ${EXPERIMENT}"
else
    echo "No EXPERIMENT specified — discovering all experiments with minus bigwigs..."
    mapfile -t EXPERIMENTS < <(
        find "${RESULTS_BASE}" -maxdepth 3 -name "*minus*.bw" ! -name "*_neg.bw" \
            -printf '%h\n' \
        | xargs -I{} bash -c 'p="{}"; while [[ "${p}" != "${RESULTS_BASE}" && "${p}" != "/" ]]; do parent=$(dirname "${p}"); if [[ "${parent}" == "${RESULTS_BASE}" ]]; then echo "${p##*/}"; break; fi; p="${parent}"; done' \
        | sort -u
    )
    if [[ ${#EXPERIMENTS[@]} -eq 0 ]]; then
        echo "No experiments with minus bigwigs found under ${RESULTS_BASE}"
        exit 0
    fi
    echo "Found experiments: ${EXPERIMENTS[*]}"
fi

# ── Process each experiment ────────────────────────────────────────────────
for EXP in "${EXPERIMENTS[@]}"; do
    RESULTS_DIR="${RESULTS_BASE}/${EXP}"
    if [[ ! -d "${RESULTS_DIR}" ]]; then
        echo "WARNING: results directory not found for '${EXP}': ${RESULTS_DIR} — skipping"
        continue
    fi
    echo ""
    echo "--- Negating minus bigwigs for: ${EXP} ---"
    bash "${SCRIPT_BASE}/scripts/negate_minus_bigwigs.sh" \
        "${RESULTS_DIR}" \
        "${CHROM_SIZES}" \
        "${THREADS}"
done

echo ""
echo "========== Job finished at: $(date) =========="
echo "Load *minus_neg.bw files in IGV and set data range to e.g. [-250, 0]."
