# Redirect stderr to stdout FIRST — before any imports — so import errors
# and all tracebacks appear in the Snakemake stdout log file rather than
# the separate .err file that is never shown in the job log.
import sys
sys.stderr = sys.stdout

import gzip
import pandas as pd
from pybedtools import BedTool
from pathlib import Path
from pandas.errors import EmptyDataError

# Column names for the finemapped sites file (Skipper finemap output).
# Used as a label prefix; the actual column count is read from the file
# dynamically so this list never needs to exactly match the file format.
FINEMAPPED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'input', 'clip', 'cdf']


def count_bed_columns(filepath):
    """Return the number of tab-separated columns in the first data line."""
    opener = gzip.open if str(filepath).endswith('.gz') else open
    with opener(str(filepath), 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            return len(line.strip().split('\t'))
    return 0


if __name__ == '__main__':
    finemapped_path, ranking_path, window_path, outf = sys.argv[1:5]

    print(f"[annotate_finemap] args: {sys.argv[1:]}", flush=True)

    n_finemap_cols = count_bed_columns(finemapped_path)
    print(f"[annotate_finemap] finemapped file has {n_finemap_cols} columns", flush=True)

    finemapped_sites = BedTool(finemapped_path)
    ranking = pd.read_csv(ranking_path, sep='\t')
    window = pd.read_csv(window_path, sep='\t')

    print(f"[annotate_finemap] window columns ({len(window.columns)}): {window.columns.tolist()}", flush=True)
    print(f"[annotate_finemap] ranking columns: {ranking.columns.tolist()}", flush=True)

    try:
        finemapped_sites.to_dataframe()
        print("[annotate_finemap] finemapped_sites non-empty, proceeding", flush=True)
    except EmptyDataError:
        print("No sites to parse. Creating an empty file and exiting cleanly.")
        Path(outf).touch()
        sys.exit(0)

    # Build finemap column names to match the actual file, falling back to
    # generic names if the file has more columns than FINEMAPPED_COLS covers.
    finemap_col_names = (
        FINEMAPPED_COLS + [f'finemap_col{i}' for i in range(len(FINEMAPPED_COLS), n_finemap_cols)]
    )[:n_finemap_cols]

    window_bed = BedTool.from_dataframe(window)

    # find windows overlapping finemapped regions
    intersection = finemapped_sites.intersect(window_bed, s=True, wb=True)
    try:
        site_annotation_finemap = intersection.to_dataframe()
        print(f"[annotate_finemap] intersection shape: {site_annotation_finemap.shape}", flush=True)
    except EmptyDataError:
        print("No annotated overlaps. Creating an empty file and exiting cleanly.")
        Path(outf).touch()
        sys.exit(0)

    window_cols = window.columns.tolist()
    expected_cols = finemap_col_names + [c + '.1' for c in window_cols[:6]] + window_cols[6:]

    actual_n = len(site_annotation_finemap.columns)
    expected_n = len(expected_cols)
    print(f"[annotate_finemap] expected {expected_n} cols, got {actual_n}", flush=True)
    print(f"[annotate_finemap] expected_cols: {expected_cols}", flush=True)

    if actual_n != expected_n:
        print(
            f"[annotate_finemap] FATAL: column count mismatch "
            f"(expected {expected_n}, intersection has {actual_n}). "
            f"finemapped has {n_finemap_cols} cols, window has {len(window_cols)} cols.",
            flush=True,
        )
        sys.exit(1)

    site_annotation_finemap.columns = expected_cols

    site_annotation_finemap['overlap_size'] = (
        site_annotation_finemap['end'] - site_annotation_finemap['start']
    )
    site_annotation_finemap['gene_types_rank'] = site_annotation_finemap['gene_types'].map(
        ranking.set_index('accession')['rank']
    )

    # get unique 1-to-1 mapping, prioritize by ranking and overlap size
    site_annotation_uniq = site_annotation_finemap.sort_values(
        by='overlap_size', ascending=False
    ).sort_values(
        by='gene_types_rank', ascending=True
    ).drop_duplicates(subset='name', keep='first')

    site_annotation_uniq.to_csv(outf, sep='\t', index=False)
    print("[annotate_finemap] done.", flush=True)
