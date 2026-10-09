#!/usr/bin/env bash
# =============================================================================
# negate_minus_bigwigs.sh
#
# Post-processing: negate all minus-strand bigwigs in a Skipper results
# directory so that minus-strand peaks face downward in IGV.
#
# For every *minus*.bw (excluding already-negated *_neg.bw), converts to
# bedGraph, multiplies signal by -1, sorts, and writes *minus_neg.bw
# alongside the original.  Both scaled and unscaled bigwigs are processed.
#
# Usage:
#   bash scripts/negate_minus_bigwigs.sh <results_dir> <chrom_sizes> [threads]
#
# Arguments:
#   results_dir  — Skipper output directory (e.g. results/rhine_pkr)
#   chrom_sizes  — two-column chr<TAB>length file
#                  (e.g. references/hg38/star_index/chrNameLength.txt)
#   threads      — parallel jobs (default: 4)
#
# Requires: bigWigToBedGraph, bedGraphToBigWig (UCSC tools; in eclip conda env)
#
# Example:
#   bash scripts/negate_minus_bigwigs.sh \
#     /xdisk/haining/maarowosegbe/eclip_pipeline/results/rhine_pkr \
#     /xdisk/haining/maarowosegbe/eclip_pipeline/references/hg38/star_index/chrNameLength.txt \
#     8
# =============================================================================
set -euo pipefail

RESULTS_DIR="${1:?ERROR: results_dir required. Usage: $0 <results_dir> <chrom_sizes> [threads]}"
CHROM_SIZES="${2:?ERROR: chrom_sizes required. Usage: $0 <results_dir> <chrom_sizes> [threads]}"
THREADS="${3:-4}"

# ── Validate inputs ────────────────────────────────────────────────────────
[[ -d "${RESULTS_DIR}" ]] || { echo "ERROR: not a directory: ${RESULTS_DIR}"; exit 1; }
[[ -f "${CHROM_SIZES}" ]] || { echo "ERROR: chrom sizes file not found: ${CHROM_SIZES}"; exit 1; }

for tool in bigWigToBedGraph bedGraphToBigWig; do
    command -v "${tool}" &>/dev/null \
        || { echo "ERROR: ${tool} not found — activate the eclip conda env first"; exit 1; }
done

# ── Find minus bigwigs (skip already-negated files) ────────────────────────
mapfile -t MINUS_BWS < <(
    find "${RESULTS_DIR}" -name "*minus*.bw" ! -name "*_neg.bw" | sort
)

if [[ ${#MINUS_BWS[@]} -eq 0 ]]; then
    echo "No minus-strand bigwig files found under: ${RESULTS_DIR}"
    exit 0
fi

echo "Found ${#MINUS_BWS[@]} minus-strand bigwig(s) — negating with ${THREADS} parallel job(s)"
echo "Chrom sizes: ${CHROM_SIZES}"
echo

# ── Per-file negation function ─────────────────────────────────────────────
negate_one() {
    local bw="$1"
    local out="${bw%.bw}_neg.bw"

    if [[ -f "${out}" ]]; then
        printf 'SKIP  (exists) %s\n' "$(basename "${out}")"
        return 0
    fi

    printf 'START %s\n' "$(basename "${bw}")"
    bigWigToBedGraph "${bw}" stdout \
        | awk 'BEGIN{OFS="\t"} {$4 = -$4; print}' \
        | LC_COLLATE=C sort -k1,1 -k2,2n \
        | bedGraphToBigWig stdin "${CHROM_SIZES}" "${out}"
    printf 'DONE  %s\n' "$(basename "${out}")"
}

export -f negate_one
export CHROM_SIZES   # inherited by each bash subprocess spawned by xargs

# ── Run in parallel using null-delimited xargs ─────────────────────────────
printf '%s\0' "${MINUS_BWS[@]}" \
    | xargs -0 -P "${THREADS}" -I {} bash -c 'negate_one "$@"' _ {}

echo
echo "Done. Load the *minus_neg.bw files in IGV and set data range to e.g. [-250, 0]."
