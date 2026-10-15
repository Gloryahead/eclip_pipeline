#!/usr/bin/env python3
"""
Patch Skipper's se_preprocess.smk — two JVM flag fixes:

1. -Xss: increase Java stack from 40M to 512m.
   Directional dedup uses recursive DFS over UMI edit-distance graph; the
   default 40M overflows at high-depth rRNA positions.

2. -Xms: decouple initial heap from max heap.
   Skipper sets Xms=Xmx, pre-committing the full heap at startup.
   Starting at 4 GB instead avoids cgroup OOM kills on startup.

Idempotent: skips if both replacements are already in place.

Usage:
    python3 scripts/patch_se_preprocess_rule.py skipper/rules/se_preprocess.smk
"""

import sys
import re

smk_path = sys.argv[1]

with open(smk_path) as fh:
    src = fh.read()

if '-Xss512m' in src and '-Xms4000M -Xmx' in src:
    print("patch_se_preprocess_rule.py: already patched, skipping.")
    sys.exit(0)

patched = re.sub(r'-Xss[0-9]+[mMgGkK]', '-Xss512m', src)
patched = patched.replace('-Xms{resources.java_mem_mb}M -Xmx', '-Xms4000M -Xmx')

if patched == src:
    print(
        "patch_se_preprocess_rule.py: WARNING — patterns not found, file unchanged.",
        file=sys.stderr,
    )
    sys.exit(1)

with open(smk_path, 'w') as fh:
    fh.write(patched)

print(f"patch_se_preprocess_rule.py: patched {smk_path}.")
