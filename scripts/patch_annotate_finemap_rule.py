#!/usr/bin/env python3
"""
Patch skipper/rules/finemap.smk to use a hardcoded path for
annotate_finemapped_regions.py instead of {TOOL_DIR}.

Why: {TOOL_DIR}/annotate_finemapped_regions.py points into the skipper
submodule on NFS. Compute nodes cache that inode and serve stale content
when the parent job patches the file but the child job starts before the
NFS attribute cache (typically 60 s) expires. The scripts/ version has a
fresh inode that no node has previously cached.

The patch rewrites the shell command so Snakemake compiles the child job's
shell script with the scripts/ path. Since Snakemake compiles shell templates
on the parent node at DAG-build time (after this patch runs on the same
node), close-to-open NFS consistency guarantees the parent node sees the
patched finemap.smk. Child SLURM jobs receive the pre-compiled shell script
with the hardcoded path and open the fresh inode.
"""
import sys

smk_path = sys.argv[1]
script_path = sys.argv[2]  # absolute path to our annotate_finemapped_regions.py

with open(smk_path) as fh:
    src = fh.read()

OLD = '{TOOL_DIR}/annotate_finemapped_regions.py'
NEW = script_path

if OLD not in src:
    if script_path in src:
        print("already patched, skipping.")
        sys.exit(0)
    print("WARNING — pattern not found, file unchanged.", file=sys.stderr)
    sys.exit(0)

patched = src.replace(OLD, NEW)
with open(smk_path, 'w') as fh:
    fh.write(patched)
print(f"patched {smk_path}.")
