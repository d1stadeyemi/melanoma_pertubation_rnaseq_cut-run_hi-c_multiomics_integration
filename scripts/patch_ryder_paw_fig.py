#!/usr/bin/env python
"""
Patch Ryder's paw.py so `_3_fgSignalConversion.pdf` stops overlapping titles.
Three surgical edits: enlarge the 2x2 figure + constrained_layout; shrink the
long equation title; drop tight_layout (conflicts with constrained_layout).
Asserts each target matches exactly once; safe to re-run (skips if already done).
Usage:  python scripts/ryder/patch_paw_fig.py /path/to/paw.py
"""
import sys
f = sys.argv[1]
s = open(f).read()

edits = [
    ("    fig, axs = pylab.subplots(2, 2, figsize=(4, 4))",
     "    fig, axs = pylab.subplots(2, 2, figsize=(10, 8), constrained_layout=True)"),
    ('''            "after normalization M~A PCC:%.3f\\nlog2(%s)=%.3flog2(%.3f%s)+%.3f"
            % (diff_scaled.corr(avg_scaled), tgtLabel, alpha, fg_sf, tgtLabel,
               beta))''',
     '''            "after normalization M~A PCC:%.3f\\nlog2(%s)=%.3flog2(%.3f%s)+%.3f"
            % (diff_scaled.corr(avg_scaled), tgtLabel, alpha, fg_sf, tgtLabel,
               beta), fontsize=8)'''),
    ('''            "after normalization M~A PCC:%.3f\\nlog2(%s)=%.3flog2(%.3f%s)%.3f" %
            (diff_scaled.corr(avg_scaled), tgtLabel, alpha, fg_sf, tgtLabel,
             beta))''',
     '''            "after normalization M~A PCC:%.3f\\nlog2(%s)=%.3flog2(%.3f%s)%.3f" %
            (diff_scaled.corr(avg_scaled), tgtLabel, alpha, fg_sf, tgtLabel,
             beta), fontsize=8)'''),
    ('''    pylab.tight_layout()
    pylab.savefig(out_prefix + "_3_fgSignalConversion.pdf")''',
     '''    pylab.savefig(out_prefix + "_3_fgSignalConversion.pdf")'''),
]
for old, new in edits:
    if new in s and old not in s:
        print(f"[skip] already patched: {old[:48]!r}..."); continue
    n = s.count(old)
    assert n == 1, f"expected 1 match, found {n} for:\n{old[:120]}"
    s = s.replace(old, new); print(f"[ok]   patched: {old[:48]!r}...")
open(f, "w").write(s); print("done.")