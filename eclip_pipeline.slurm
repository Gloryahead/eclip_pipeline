#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources to run your job
# --------------------------------------------------------------

### Optional. Set the job name
#SBATCH --job-name=eclip_pipeline
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=/xdisk/haining/maarowosegbe/eclip_pipeline/logs/%x-%j.out
### REQUIRED. Specify the PI group for this job.
#SBATCH --account=haining
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the number of CPUs that will be used for this job.
#SBATCH --ntasks=32
### REQUIRED. Set the memory required for this job.
#SBATCH --mem-per-cpu=8gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=240:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

# Scripts and pipeline code live in home (quota-safe, version-controlled)
SCRIPT_BASE=/home/u11/maarowosegbe/eclip_pipeline
# Data, logs, and results go to xdisk (large storage)
DATA_BASE=/xdisk/haining/maarowosegbe/eclip_pipeline

mkdir -p "${DATA_BASE}/logs"

set +eu; source ~/.bashrc; set -eu
mamba-haining                          # sets MAMBA_ROOT_PREFIX to haining group prefix
micromamba activate eclip              # actually activates the eclip environment
export PATH="${CONDA_PREFIX}/bin:${PATH}"
echo "CONDA_PREFIX: ${CONDA_PREFIX}"
echo "Rscript: $(which Rscript 2>/dev/null || echo NOT FOUND)"

# Function to run script if it exists
run_if_exists() {
    SCRIPT_PATH=$1
    shift
    if [[ -f "$SCRIPT_PATH" ]]; then
        echo "Running: $SCRIPT_PATH $*"
        bash "$SCRIPT_PATH" "$@"
    else
        echo "WARNING: Script not found: $SCRIPT_PATH"
    fi
}

# Step 1 — Download reference genomes and build STAR indices
run_if_exists "${SCRIPT_BASE}/setup/01_download_references.sh"

# Step 3 — Download all 31 eCLIP IP FASTQ files from SRA (GSE276985)
run_if_exists "${SCRIPT_BASE}/scripts/download_geo_fastqs.sh" all

# Step 4 — Run Skipper eCLIP pipeline: human experiments (hg38)
run_if_exists "${SCRIPT_BASE}/run_analysis.sh" human

# Step 5 — Run Skipper eCLIP pipeline: mouse experiments (mm10)
run_if_exists "${SCRIPT_BASE}/run_analysis.sh" mouse

echo "========== Job finished at: $(date) =========="
