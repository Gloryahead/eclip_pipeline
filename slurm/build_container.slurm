#!/bin/bash
#SBATCH --job-name=build_sif
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/xdisk/haining/maarowosegbe/scATAC-seq/logs/build_sif-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/scATAC-seq/logs/build_sif-%j.err

echo "========== Container build started: $(date) =========="

module --ignore_cache load apptainer 2>/dev/null || true

mkdir -p /xdisk/haining/maarowosegbe/scATAC-seq/containers

apptainer build \
    /xdisk/haining/maarowosegbe/scATAC-seq/containers/r_signac.sif \
    /home/u11/maarowosegbe/scATAC-seq_pipeline/containers/r_signac.def

echo "========== Container build finished: $(date) =========="

# Quick test
echo "Testing container..."
apptainer exec /xdisk/haining/maarowosegbe/scATAC-seq/containers/r_signac.sif \
    Rscript -e 'library(Signac); cat("Signac version:", as.character(packageVersion("Signac")), "\n")'
