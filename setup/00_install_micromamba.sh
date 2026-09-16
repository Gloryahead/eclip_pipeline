#!/usr/bin/env bash
# =============================================================================
# 00_install_micromamba.sh
# One-time setup: creates the 'eclip' micromamba environment containing
# every tool needed for the Skipper eCLIP pipeline.
#
# Prerequisites:
#   - mamba-haining is available in your shell (UArizona HPC: already set up)
#   - Skipper has been cloned:
#       git clone https://github.com/YeoLab/skipper.git \
#           /home/u11/maarowosegbe/eclip_pipeline/skipper
#
# Usage (run from the login node or an interactive session):
#   bash setup/00_install_micromamba.sh
#
# This script is safe to re-run; it skips creation if the env already exists.
# =============================================================================
set -euo pipefail

set +eu; source ~/.bashrc; set -eu
mamba-haining    # sets MAMBA_ROOT_PREFIX to /groups/haining/maarowosegbe/micromamba

ENV_YAML="$(dirname "$0")/../envs/eclip.yaml"

echo "=== Creating 'eclip' micromamba environment ==="
if micromamba env list | grep -q "^eclip "; then
    echo "'eclip' environment already exists. To rebuild, remove it first:"
    echo "  micromamba env remove -n eclip"
else
    micromamba env create -y -f "${ENV_YAML}"
    echo "'eclip' environment created successfully."
fi

echo ""
echo "=== Verifying key tools ==="
micromamba run -n eclip snakemake --version
micromamba run -n eclip STAR --version 2>&1 | head -1
micromamba run -n eclip umi_tools --version
echo ""
echo "=== Setup complete ==="
echo "Next: sbatch setup/01_download_references.sh"
