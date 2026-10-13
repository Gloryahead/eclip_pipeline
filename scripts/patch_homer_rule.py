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
import re

smk_path = sys.argv[1]

with open(smk_path) as fh:
    src = fh.read()

OLD = (
    r"fg_input=<\(zcat \{input\.finemapped_windows\}"
    r".*?"
    r'\| awk -v OFS="\\t" \'{{print \$4, \$1, \$2\+1, \$3, \$6}}\'\)'
    r"\s*\nbg_input=<\(zcat \{input\.background\}"
    r".*?"
    r'\| awk -v OFS="\\t" \'{{print \$4, \$1, \$2\+1, \$3, \$6}}\'\)'
)

# Build replacement: temp files for both fg and bg
NEW = (
    'fg_tmp=$(mktemp --suffix=.homer_fg.bed)\n'
    '    bg_tmp=$(mktemp --suffix=.homer_bg.bed)\n'
    '    zcat {input.finemapped_windows} \\\n'
    '        | awk -v OFS="\\t" \'{{print $4 ":" $9, $1, $2+1, $3, $6}}\' > "$fg_tmp"\n'
    '    zcat {input.background} \\\n'
    '        | awk -v OFS="\\t" \'{{print $4, $1, $2+1, $3, $6}}\' > "$bg_tmp"'
)

patched, count = re.subn(OLD, NEW, src, flags=re.DOTALL)

if count == 0:
    # Already patched or pattern changed — check for fg_tmp to avoid double-patch
    if 'fg_tmp=$(mktemp' in src:
        print("patch_homer_rule.py: already patched, skipping.")
        sys.exit(0)
    print("patch_homer_rule.py: WARNING — pattern not found, file unchanged.", file=sys.stderr)
    sys.exit(0)

# Also replace the findMotifsGenome.pl call to use temp files and clean up
patched = patched.replace(
    'findMotifsGenome.pl "$fg_input"',
    'findMotifsGenome.pl "$fg_tmp"'
)
patched = patched.replace(
    '-bg "$bg_input"',
    '-bg "$bg_tmp"\n    rm -f "$fg_tmp" "$bg_tmp"'
)

with open(smk_path, 'w') as fh:
    fh.write(patched)

print(f"patch_homer_rule.py: patched {smk_path} ({count} substitution(s)).")
