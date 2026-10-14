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

if 'fg_tmp=$(mktemp' in src:
    print("patch_homer_rule.py: already patched, skipping.")
    sys.exit(0)

# Use a flexible regex to match the process-substitution block regardless of
# exact indentation or trailing whitespace — the exact string match proved
# fragile across environments.
m = re.search(
    r'( +)fg_input=<\(zcat \{input\.finemapped_windows\} \\\n'
    r'( +)\| awk[^\n]+\)\n'
    r'\n'
    r' +bg_input=<\(zcat \{input\.background\} \\\n'
    r' +\| awk[^\n]+\)',
    src,
)

if m is None:
    print(
        "patch_homer_rule.py: WARNING — process-substitution block not found, file unchanged.",
        file=sys.stderr,
    )
    sys.exit(1)

indent = m.group(1)   # leading spaces on fg_input line (e.g. 8 spaces)
indent2 = m.group(2)  # leading spaces on | awk line (e.g. 12 spaces)

new_block = (
    indent + 'fg_tmp=$(mktemp --suffix=.homer_fg.bed)\n' +
    indent + 'bg_tmp=$(mktemp --suffix=.homer_bg.bed)\n' +
    indent + 'zcat {input.finemapped_windows} \\\n' +
    indent2 + '| awk -v OFS="\\t" \'{{print $4 ":" $9, $1, $2+1, $3, $6}}\' > "$fg_tmp"\n' +
    indent + 'zcat {input.background} \\\n' +
    indent2 + '| awk -v OFS="\\t" \'{{print $4, $1, $2+1, $3, $6}}\' > "$bg_tmp"'
)

patched = src[:m.start()] + new_block + src[m.end():]

# Replace variable references in the findMotifsGenome.pl call
patched = patched.replace('"$fg_input"', '"$fg_tmp"', 1)
patched = patched.replace('"$bg_input"', '"$bg_tmp"', 1)

# Add rm cleanup after the log redirect that ends the HOMER call
log_redirect = '>> {log.stdout} 2> {log.stderr}'
idx_bg = patched.find('-bg "$bg_tmp"')
if idx_bg >= 0:
    idx_log = patched.find(log_redirect, idx_bg)
    if idx_log >= 0:
        insert_at = idx_log + len(log_redirect)
        patched = (
            patched[:insert_at]
            + '\n' + indent + 'rm -f "$fg_tmp" "$bg_tmp"'
            + patched[insert_at:]
        )

with open(smk_path, 'w') as fh:
    fh.write(patched)

print(f"patch_homer_rule.py: patched {smk_path}.")
