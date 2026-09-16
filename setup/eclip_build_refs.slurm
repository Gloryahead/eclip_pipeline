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
### 16 CPUs needed for STAR genome indexing.
#SBATCH --ntasks=16
### REQUIRED. Set the memory required for this job.
### STAR genome indexing requires ~64 GB total; 4 gb x 16 CPUs = 64 GB.
#SBATCH --mem-per-cpu=4gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=08:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu; source ~/.bashrc; set -eu
mamba-haining
micromamba activate eclip
export PATH="${CONDA_PREFIX}/bin:${PATH}"

bash /home/u11/maarowosegbe/eclip_pipeline/setup/01_download_references.sh

echo "========== Job finished at: $(date) =========="
