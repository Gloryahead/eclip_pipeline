#!/bin/bash
#SBATCH --job-name=dl_refs
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --time=02:00:00
#SBATCH --output=/xdisk/haining/maarowosegbe/scATAC-seq/logs/dl_refs-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/scATAC-seq/logs/dl_refs-%j.err

echo "========== Reference download started: $(date) =========="

mkdir -p /xdisk/haining/maarowosegbe/scATAC-seq/logs

bash /home/u11/maarowosegbe/scATAC-seq_pipeline/scripts/bash/download_reference.sh

echo "========== Reference download finished: $(date) =========="
