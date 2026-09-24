#!/usr/bin/env bash
# =============================================================================
# run_analysis.sh
# Main entry point for the Rhine et al. 2025 eCLIP-seq pipeline.
# Runs the Skipper pipeline via Snakemake on the UArizona HPC (SLURM).
#
# Usage:
#   bash run_analysis.sh human    # Run all human eCLIP experiments (hg38)
#   bash run_analysis.sh mouse    # Run mouse age-series eCLIP (mm10)
#   bash run_analysis.sh human --dry-run   # Preview without submitting jobs
#   bash run_analysis.sh human --unlock    # Unlock a stale working directory
#
# Prerequisites (run once before first use):
#   1. bash setup/00_install_micromamba.sh
#   2. sbatch setup/01_download_references.sh
#   3. Fill in FASTQ paths in config/manifest_human.csv or config/manifest_mouse.csv
# =============================================================================
set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR="/xdisk/haining/maarowosegbe/eclip_pipeline"
SCRIPT_BASE="/home/u11/maarowosegbe/eclip_pipeline"
SKIPPER_DIR="${SCRIPT_BASE}/skipper"
LOG_DIR="${BASE_DIR}/logs"

# env already activated by eclip_pipeline.slurm; snakemake and all tools on PATH

GENOME="${1:-}"
shift || true
EXTRA_ARGS=("$@")

if [[ -z "${GENOME}" ]]; then
    echo "Usage: bash run_analysis.sh [human|mouse] [--dry-run|--unlock|...]"
    exit 1
fi

case "${GENOME}" in
    human)
        CONFIG="${SCRIPT_BASE}/config/skipper_config_human.yaml"
        ;;
    mouse)
        CONFIG="${SCRIPT_BASE}/config/skipper_config_mouse.yaml"
        ;;
    g3bp1_test)
        CONFIG="${SCRIPT_BASE}/config/skipper_config_g3bp1_test.yaml"
        ;;
    *)
        echo "Unknown genome '${GENOME}'. Use 'human', 'mouse', or 'g3bp1_test'."
        exit 1
        ;;
esac

mkdir -p "${LOG_DIR}"

echo "=== Running Skipper eCLIP pipeline ==="
echo "Genome  : ${GENOME}"
echo "Config  : ${CONFIG}"
echo "Skipper : ${SKIPPER_DIR}"
echo "Env     : ${CONDA_DEFAULT_ENV:-not set — check mamba-haining activate above}"
echo "Commit  : $(git -C "${SKIPPER_DIR}" rev-parse HEAD 2>/dev/null || echo unknown)"
echo ""

# Unset SLURM_JOB_ID — if this script is called from inside a SLURM job
# (e.g. via eclip_pipeline.slurm), Snakemake sees it and thinks it is
# already inside a cluster job, which prevents it from submitting children.
unset SLURM_JOB_ID

snakemake \
    --snakefile "${SKIPPER_DIR}/Skipper.py" \
    --configfile "${CONFIG}" \
    --profile "${SCRIPT_BASE}/config/profiles/slurm" \
    --rerun-incomplete \
    --keep-going \
    "${EXTRA_ARGS[@]}"
