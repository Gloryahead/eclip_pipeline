# to annotate finemapped regions back to window's annotations
import pandas as pd
from pybedtools import BedTool
import sys
from pathlib import Path
from pandas.errors import EmptyDataError

# pybedtools >=0.10.0 regression: to_dataframe(names=...) silently ignores the
# names parameter and still returns integer RangeIndex. Fix: call to_dataframe()
# without names, then rename columns directly on the resulting DataFrame.
FINEMAPPED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'input', 'clip', 'cdf']

if __name__ == '__main__':
    finemapped_sites = BedTool(sys.argv[1])

    ranking = pd.read_csv(sys.argv[2], sep='\t')
    window = pd.read_csv(sys.argv[3], sep='\t')

    outf = sys.argv[4]

    try:
        finemapped_sites.to_dataframe()
    except EmptyDataError:
        print("No sites to parse. Creating an empty file and exiting cleanly with this warning.")
        Path(outf).touch()
        sys.exit(0)

    window_bed = BedTool.from_dataframe(window)

    # find windows overlapping finemapped regions
    intersection = finemapped_sites.intersect(window_bed, s=True, wb=True)
    try:
        site_annotation_finemap = intersection.to_dataframe()
    except EmptyDataError:
        print("No annotated overlaps. Creating an empty file and exiting cleanly with this warning.")
        Path(outf).touch()
        sys.exit(0)

    expected_cols = (FINEMAPPED_COLS +
                     [i + '.1' for i in window.columns.tolist()[:6]] +
                     window.columns.tolist()[6:])
    site_annotation_finemap.columns = expected_cols

    site_annotation_finemap['overlap_size'] = site_annotation_finemap['end'] - site_annotation_finemap['start']
    site_annotation_finemap['gene_types_rank'] = site_annotation_finemap['gene_types'].map(ranking.set_index('accession')['rank'])

    # get unique 1-to-1 mapping, prioritize by ranking and overlap size
    site_annotation_uniq = site_annotation_finemap.sort_values(
        by='overlap_size', ascending=False
    ).sort_values(
        by='gene_types_rank', ascending=True
    ).drop_duplicates(subset='name', keep='first')

    site_annotation_uniq.to_csv(outf, sep='\t', index=False)
