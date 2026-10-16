# Redirect stderr to stdout FIRST — before any imports — so import errors
# and all tracebacks appear in the Snakemake stdout log file rather than
# the separate .err file that is never shown in the job log.
import sys
sys.stderr = sys.stdout

import gzip
import os
import shutil
import subprocess
import tempfile
import pandas as pd
from pathlib import Path

# Column names for the finemapped sites file (Skipper finemap output).
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


def is_bed_empty(filepath):
    """Return True if the BED/BED.gz file has no data lines."""
    opener = gzip.open if str(filepath).endswith('.gz') else open
    with opener(str(filepath), 'rt') as fh:
        for line in fh:
            if line.strip() and not line.startswith('#'):
                return False
    return True


if __name__ == '__main__':
    finemapped_path, ranking_path, window_path, outf = sys.argv[1:5]

    print(f"[annotate_finemap] args: {sys.argv[1:]}", flush=True)

    n_finemap_cols = count_bed_columns(finemapped_path)
    print(f"[annotate_finemap] finemapped file has {n_finemap_cols} columns", flush=True)

    if is_bed_empty(finemapped_path):
        print("No sites to parse. Creating an empty file and exiting cleanly.")
        Path(outf).parent.mkdir(parents=True, exist_ok=True)
        Path(outf).touch()
        sys.exit(0)
    print("[annotate_finemap] finemapped_sites non-empty, proceeding", flush=True)

    ranking = pd.read_csv(ranking_path, sep='\t')
    print(f"[annotate_finemap] ranking columns: {ranking.columns.tolist()}", flush=True)

    # Read ONLY the header from the large annotation file.
    # Loading the full 3.45 GB TSV via pd.read_csv() uses 15-20 GB RAM
    # (pandas object dtype for string columns) — enough to trigger OOM kill.
    with open(window_path) as f:
        header_line = f.readline()
    window_cols = header_line.strip().split('\t')
    print(f"[annotate_finemap] window columns ({len(window_cols)}): {window_cols}", flush=True)

    # Stream the annotation file (skip header row) to a headerless BED in TMPDIR.
    # shutil.copyfileobj reads in 16 KB chunks — O(1) memory regardless of file size.
    tmpdir = os.environ.get('TMPDIR', tempfile.gettempdir())
    os.makedirs(tmpdir, exist_ok=True)
    window_bed_path = os.path.join(tmpdir, 'skipper_annotation_window.bed')

    with open(window_path) as src, open(window_bed_path, 'w') as dst:
        src.readline()  # discard header
        shutil.copyfileobj(src, dst)
    window_bed_size = os.path.getsize(window_bed_path)
    print(f"[annotate_finemap] streamed window BED to {window_bed_path} ({window_bed_size:,} bytes)", flush=True)

    # Call bedtools intersect directly via subprocess.
    # pybedtools' .intersect() silently returned (0, 0): the bedtools subprocess
    # was OOM-killed (peak RSS ~21 GB for Python + bedtools combined), leaving an
    # empty temp file that pybedtools returned as a (0, 0) DataFrame with no error.
    # Direct bedtools call confirmed working on Puma (2026-10-16 manual test).
    bedtools_bin = os.path.join(os.environ.get('CONDA_PREFIX', ''), 'bin', 'bedtools')
    if not os.path.exists(bedtools_bin):
        bedtools_bin = 'bedtools'

    intersection_path = os.path.join(tmpdir, 'skipper_intersection.bed')
    print(f"[annotate_finemap] running bedtools intersect (direct subprocess)", flush=True)

    with open(intersection_path, 'w') as intersection_out:
        bt_result = subprocess.run(
            [bedtools_bin, 'intersect',
             '-a', str(finemapped_path),
             '-b', window_bed_path,
             '-s', '-wb'],
            stdout=intersection_out,
            stderr=subprocess.PIPE,
            text=True,
        )

    print(f"[annotate_finemap] bedtools returncode: {bt_result.returncode}", flush=True)
    if bt_result.stderr:
        print(f"[annotate_finemap] bedtools stderr: {bt_result.stderr[:2000]}", flush=True)
    if bt_result.returncode != 0:
        print(f"[annotate_finemap] FATAL: bedtools intersect failed", flush=True)
        sys.exit(1)

    intersection_size = os.path.getsize(intersection_path)
    print(f"[annotate_finemap] intersection output: {intersection_size:,} bytes", flush=True)

    if intersection_size == 0:
        print("No annotated overlaps. Creating an empty file and exiting cleanly.")
        Path(outf).parent.mkdir(parents=True, exist_ok=True)
        Path(outf).touch()
        sys.exit(0)

    site_annotation_finemap = pd.read_csv(intersection_path, sep='\t', header=None)
    print(f"[annotate_finemap] intersection shape: {site_annotation_finemap.shape}", flush=True)

    finemap_col_names = (
        FINEMAPPED_COLS + [f'finemap_col{i}' for i in range(len(FINEMAPPED_COLS), n_finemap_cols)]
    )[:n_finemap_cols]

    expected_cols = finemap_col_names + [c + '.1' for c in window_cols[:6]] + window_cols[6:]
    actual_n = len(site_annotation_finemap.columns)
    expected_n = len(expected_cols)
    print(f"[annotate_finemap] expected {expected_n} cols, got {actual_n}", flush=True)

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

    site_annotation_uniq = site_annotation_finemap.sort_values(
        by='overlap_size', ascending=False
    ).sort_values(
        by='gene_types_rank', ascending=True
    ).drop_duplicates(subset='name', keep='first')

    Path(outf).parent.mkdir(parents=True, exist_ok=True)
    site_annotation_uniq.to_csv(outf, sep='\t', index=False)
    print("[annotate_finemap] done.", flush=True)
