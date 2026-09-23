#!/bin/bash
# =============================================================================
# run_pipeline.sh — SLURM job to launch the Nextflow scATAC-seq pipeline
#
# This is the one script you submit to SLURM. It starts Nextflow, which then
# submits every pipeline step as its own separate SLURM job automatically.
# You do NOT need to write individual sbatch scripts for each step.
#
# Submit with:
#   sbatch /home/u11/maarowosegbe/scATAC-seq_pipeline/slurm/run_pipeline.sh
#
# To resume after a failure:
#   sbatch /home/u11/maarowosegbe/scATAC-seq_pipeline/slurm/run_pipeline.sh --resume
# =============================================================================

# ── SLURM resource request ────────────────────────────────────────────────────
#SBATCH --job-name=scatac_pipeline
#SBATCH --output=/xdisk/haining/maarowosegbe/scATAC-seq/logs/scatac_pipeline-%j.out
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=240:00:00
# This job just runs Nextflow (lightweight). Nextflow submits the heavy
# compute steps (Cell Ranger, Signac, etc.) as separate SLURM jobs itself.

echo "========== Pipeline started: $(date) =========="

# ── Directory layout ──────────────────────────────────────────────────────────
PIPELINE_DIR=/home/u11/maarowosegbe/scATAC-seq_pipeline
DATA_DIR=/xdisk/haining/maarowosegbe/scATAC-seq

mkdir -p "${DATA_DIR}/logs"
mkdir -p "${DATA_DIR}/results"
mkdir -p "${DATA_DIR}/nextflow_work"

# ── Activate micromamba (matches your ~/.bashrc setup on Puma) ────────────────
set +eu
source ~/.bashrc
mamba-haining                          # sets MAMBA_ROOT_PREFIX to haining group
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow           # the env you created: micromamba create -n nextflow
set -eu

# ── Load Apptainer (available as a module on UA Puma) ─────────────────────────
module load apptainer

echo "Nextflow version: $(nextflow -version 2>&1 | head -1)"
echo "Apptainer version: $(apptainer --version)"

# ── Run Nextflow ──────────────────────────────────────────────────────────────
# -profile slurm  : use SLURM executor + Apptainer (defined in nextflow.config)
# -resume         : re-use completed steps if the pipeline was interrupted
# Pass --resume as a command-line argument to this sbatch script to enable it

RESUME_FLAG=""
if [[ "${1:-}" == "--resume" ]]; then
    RESUME_FLAG="-resume"
    echo "Resume mode enabled — completed steps will be cached."
fi

nextflow run "${PIPELINE_DIR}/main.nf" \
    -c "${PIPELINE_DIR}/nextflow.config" \
    -profile slurm \
    -work-dir "${DATA_DIR}/nextflow_work" \
    ${RESUME_FLAG} \
    2>&1 | tee "${DATA_DIR}/logs/nextflow_$(date +%Y%m%d_%H%M%S).log"

EXIT_CODE=${PIPESTATUS[0]}

echo "========== Pipeline finished: $(date) =========="
echo "Exit code: ${EXIT_CODE}"

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "SUCCESS — results are in: ${DATA_DIR}/results/"
    echo "Pipeline report: ${DATA_DIR}/results/pipeline_report.html"
else
    echo "FAILED — check log above and re-submit with: sbatch run_pipeline.sh --resume"
fi

exit ${EXIT_CODE}
