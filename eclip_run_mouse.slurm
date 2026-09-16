#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources to run your job
# --------------------------------------------------------------

### Optional. Set the job name
#SBATCH --job-name=eclip_run_mouse
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
### This job only runs Snakemake; actual compute is in child SLURM jobs.
#SBATCH --ntasks=2
### REQUIRED. Set the memory required for this job.
#SBATCH --mem-per-cpu=8gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
### Keep alive long enough for all child jobs to finish (~24 hrs typical).
#SBATCH --time=48:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/eclip_pipeline/logs

set +eu; source ~/.bashrc; set -eu
mamba-haining
micromamba activate eclip
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# Snakemake will submit child SLURM jobs automatically.
bash /home/u11/maarowosegbe/eclip_pipeline/run_analysis.sh mouse

echo "========== Job finished at: $(date) =========="
