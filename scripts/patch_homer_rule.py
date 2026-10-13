#!/usr/bin/env python3
"""
Patch Skipper's run_homer rule in rules/analysis.smk to replace process
substitutions with temp files.

HOMER internally calls mergePeaks as a child subprocess. Child processes
cannot read from /dev/fd/N file descriptors created by bash process
substitution in the parent shell — only the parent shell inherits them.
Result: 0 background peaks, invalid motif normalization, no output files.

Fix: write foreground and background peaks to mktemp files before calling
findMotifsGenome.pl so all subprocesses can open them as regular files.

Usage:
    python3 scripts/patch_homer_rule.py skipper/rules/analysis.smk
"""

import sys

smk_path = sys.argv[1]

with open(smk_path) as fh:
    src = fh.read()

# Already patched check — idempotent
if 'fg_tmp=$(mktemp' in src:
    print("patch_homer_rule.py: already patched, skipping.")
    sys.exit(0)

OLD_PS = (
    '        fg_input=<(zcat {input.finemapped_windows} \\\n'
    '            | awk -v OFS="\\t" \'{{print $4 ":" $9, $1, $2+1, $3, $6}}\')\n'
    '\n'
    '        bg_input=<(zcat {input.background} \\\n'
    '            | awk -v OFS="\\t" \'{{print $4, $1, $2+1, $3, $6}}\')'
)

NEW_TMP = (
    '        fg_tmp=$(mktemp --suffix=.homer_fg.bed)\n'
    '        bg_tmp=$(mktemp --suffix=.homer_bg.bed)\n'
    '        zcat {input.finemapped_windows} \\\n'
    '            | awk -v OFS="\\t" \'{{print $4 ":" $9, $1, $2+1, $3, $6}}\' > "$fg_tmp"\n'
    '        zcat {input.background} \\\n'
    '            | awk -v OFS="\\t" \'{{print $4, $1, $2+1, $3, $6}}\' > "$bg_tmp"'
)

if OLD_PS not in src:
    print("patch_homer_rule.py: WARNING — pattern not found, file unchanged.", file=sys.stderr)
    sys.exit(0)

patched = src.replace(OLD_PS, NEW_TMP, 1)

# Replace variable references in findMotifsGenome.pl call
patched = patched.replace('findMotifsGenome.pl "$fg_input"', 'findMotifsGenome.pl "$fg_tmp"', 1)
patched = patched.replace(
    '-bg "$bg_input" \\\n        >> {log.stdout} 2> {log.stderr}',
    '-bg "$bg_tmp" \\\n        >> {log.stdout} 2> {log.stderr}\n        rm -f "$fg_tmp" "$bg_tmp"',
    1
)

with open(smk_path, 'w') as fh:
    fh.write(patched)

print(f"patch_homer_rule.py: patched {smk_path}.")
