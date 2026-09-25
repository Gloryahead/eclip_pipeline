#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources to run your job
# --------------------------------------------------------------

### Optional. Set the job name
#SBATCH --job-name=eclip_build_refs
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
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=06:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

SCRIPT_BASE=/home/u11/maarowosegbe/eclip_pipeline
DATA_BASE=/xdisk/haining/maarowosegbe/eclip_pipeline

mkdir -p "${DATA_BASE}/logs"

set +eu; source ~/.bashrc; set -eu
mamba-haining
micromamba activate eclip
export PATH="${CONDA_PREFIX}/bin:${PATH}"
echo "CONDA_PREFIX: ${CONDA_PREFIX}"
echo "STAR version: $(STAR --version 2>&1 | head -1)"

bash "${SCRIPT_BASE}/setup/01_download_references.sh"

echo "========== Job finished at: $(date) =========="
