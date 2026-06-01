"""
=============================================================================
compartment_heatmap_and_structure.py
  Two compact, Q1-standard "multi-omics" visuals that COMPRESS several violin
  panels and anchor the structural story:

  FIGURE 1 — Marks x Switch-category HEATMAP
    Rows = histone marks, columns = switch categories, cell = MEAN parental
    signal, Z-SCORED PER MARK (across categories). Z-scoring per row is
    essential: raw mean signal is not comparable across marks (different
    dynamic ranges), so a raw heatmap would only show "which mark is brightest".
    Per-mark z-score shows, comparably, WHICH CATEGORY each mark is enriched/
    depleted in — the same information as the baseline violins, in one panel.
    This is a VISUAL SUMMARY; the inferential backing remains the violin +
    Kruskal-Wallis/Dunn/Cliff's-delta tests (heatmaps imply structure but do
    not test it).

  FIGURE 2 — PC1 compartment-structure plot
    Honest, fully-supported substitute for a saddle plot (which would need the
    Hi-C CONTACT MATRIX, not just the PC1 eigenvector). Shows the Parental vs
    Bottom10 PC1-value distributions and the directional shift, making the
    A/B reorganization visible as a recognizable Hi-C-signature figure from
    the data actually available.

INPUTS
  --classified  Switch_classified.txt from RF_SHAP.py (per-bin: marks + switch_cat)
                — provides the mark signal columns for the heatmap.
  --pc1         Hi-C PC1 / merged table (Parental_orig, Bottom10_orig, switch
                columns) — provides PC1 values for the structure plot.

OUTPUTS
  marks_by_switch_heatmap.pdf/png
  pc1_compartment_structure.pdf/png

USAGE
  python compartment_heatmap_and_structure.py \
      --classified results/exploratory/RF_SHAP/baseline/Switch_classified.txt \
      --pc1        data/processed/A375-CompartmentsAcrossConditions.txt \
      --outdir     results/exploratory/figures
"""

import argparse
import os
import sys

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    13,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "savefig.dpi":       600,
    "figure.facecolor":  "white",
})

SWITCH_LABELS = {"B_to_A": "B\u2192A\n(Activated)", "A_to_B": "A\u2192B\n(Silenced)",
                 "Stable_A": "Stable A", "Stable_B": "Stable B"}
CAT_ORDER = ["Stable_A", "A_to_B", "B_to_A", "Stable_B"]
CANDIDATE_MARKS = ["H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3", "CTCF"]


def classify_switch(df):
    """Same logic as RF_SHAP/Step1 (used only if switch_cat absent)."""
    cats = []
    for _, row in df.iterrows():
        cat = None
        v = str(row.get("P->B10_Both_original", "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if "toA" in v or "BtoA" in v: cat = "B_to_A"
            elif "toB" in v or "AtoB" in v: cat = "A_to_B"
        if cat is None:
            try:
                p, b = float(row.get("Parental_orig", np.nan)), float(row.get("Bottom10_orig", np.nan))
                if p > 0 and b > 0: cat = "Stable_A"
                elif p < 0 and b < 0: cat = "Stable_B"
            except (ValueError, TypeError):
                pass
        cats.append(cat)
    return cats


def parse_args():
    p = argparse.ArgumentParser(description="Marks-x-switch heatmap + PC1 structure plot.")
    p.add_argument("--classified", required=True,
                   help="Switch_classified.txt (per-bin marks + switch_cat).")
    p.add_argument("--pc1", required=True,
                   help="Hi-C PC1/merged table (Parental_orig, Bottom10_orig).")
    p.add_argument("--outdir", default="results/exploratory/figures")
    p.add_argument("--marks", nargs="+", default=None,
                   help="Marks to include (default: auto-detect among "
                        f"{CANDIDATE_MARKS}).")
    p.add_argument("--signal_suffix", default="_mean",
                   help="If mark columns are '<mark>_mean' style, the suffix to "
                        "match. Default '_mean'. Set '' for bare mark names.")
    return p.parse_args()


# =============================================================================
# FIGURE 1 — marks x switch-category heatmap (per-mark z-scored)
# =============================================================================
def build_heatmap(clf, marks, signal_suffix, outdir):
    if "switch_cat" not in clf.columns:
        sys.exit("[ERROR] 'switch_cat' not in classified file.")
    cats = [c for c in CAT_ORDER if c in clf["switch_cat"].unique()]

    # Resolve each mark to its signal column.
    def resolve(mark):
        for cand in (f"{mark}{signal_suffix}", mark):
            if cand in clf.columns:
                return cand
        # case-insensitive fallback
        lc = {c.lower(): c for c in clf.columns}
        for cand in (f"{mark}{signal_suffix}".lower(), mark.lower()):
            if cand in lc:
                return lc[cand]
        return None

    resolved = {m: resolve(m) for m in marks}
    resolved = {m: c for m, c in resolved.items() if c is not None}
    if len(resolved) < 2:
        sys.exit(f"[ERROR] <2 mark columns resolved from {marks}. "
                 f"Columns present: {list(clf.columns)}")
    print(f"[INFO] Heatmap marks -> columns: {resolved}")

    # Mean signal per (mark, category), then z-score PER MARK across categories.
    mat = np.zeros((len(resolved), len(cats)))
    mark_names = list(resolved.keys())
    for i, m in enumerate(mark_names):
        col = resolved[m]
        vals = pd.to_numeric(clf[col], errors="coerce")
        means = [np.nanmean(vals[clf["switch_cat"] == c]) for c in cats]
        means = np.array(means, dtype=float)
        mu, sd = np.nanmean(means), np.nanstd(means)
        mat[i] = (means - mu) / sd if sd > 0 else 0.0

    # Plot
    fig, ax = plt.subplots(figsize=(1.4 * len(cats) + 2.5, 0.7 * len(mark_names) + 2))
    norm = TwoSlopeNorm(vmin=min(-2, mat.min()), vcenter=0, vmax=max(2, mat.max()))
    im = ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([SWITCH_LABELS.get(c, c) for c in cats], fontsize=9)
    ax.set_yticks(range(len(mark_names)))
    ax.set_yticklabels(mark_names, fontsize=10)

    # Annotate each cell with its z-value.
    for i in range(len(mark_names)):
        for j in range(len(cats)):
            v = mat[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 1.3 else "#222222")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Parental signal\n(z-score per mark)", fontsize=9)
    ax.set_title("Parental chromatin signal across switch categories\n"
                 "(mean signal, z-scored within each mark)",
                 fontweight="bold", pad=10)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"marks_by_switch_heatmap.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)


# =============================================================================
# FIGURE 2 — PC1 compartment-structure plot (honest saddle-plot substitute)
# =============================================================================
def build_pc1_structure(pc1, outdir):
    pcol, bcol = "Parental_orig", "Bottom10_orig"
    if pcol not in pc1.columns or bcol not in pc1.columns:
        print(f"[WARN] {pcol}/{bcol} not in PC1 table; skipping structure plot.")
        return
    p = pd.to_numeric(pc1[pcol], errors="coerce")
    b = pd.to_numeric(pc1[bcol], errors="coerce")
    ok = p.notna() & b.notna()
    p, b = p[ok].values, b[ok].values

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.8))

    # Left: PC1 value distributions, Parental vs B10 (compartment-strength view).
    bins = np.linspace(min(p.min(), b.min()), max(p.max(), b.max()), 60)
    axL.hist(p, bins=bins, color="#85C1E9", alpha=0.6, label="Parental", density=True)
    axL.hist(b, bins=bins, color="#C0392B", alpha=0.5, label="Bottom10", density=True)
    axL.axvline(0, color="black", lw=0.8, linestyle="--")
    axL.set_xlabel("PC1 (compartment eigenvalue)\n\u2190 B compartment      A compartment \u2192")
    axL.set_ylabel("Density")
    axL.set_title("Compartment-strength distribution", fontweight="bold", pad=8)
    axL.legend(frameon=False, fontsize=9)
    for s in ("top", "right"): axL.spines[s].set_visible(False)

    # Right: Parental vs B10 PC1 scatter (the switch is the off-diagonal mass).
    # Quadrants: same-sign = stable; sign flip = switch.
    axR.axhline(0, color="grey", lw=0.6); axR.axvline(0, color="grey", lw=0.6)
    axR.scatter(p, b, s=3, alpha=0.15, color="#34495E", rasterized=True)
    lim = max(np.abs(np.concatenate([p, b]))) * 1.05
    axR.plot([-lim, lim], [-lim, lim], color="#C0392B", lw=0.8, linestyle="--",
             label="no change (y=x)")
    axR.set_xlim(-lim, lim); axR.set_ylim(-lim, lim)
    axR.set_xlabel("Parental PC1")
    axR.set_ylabel("Bottom10 PC1")
    axR.set_title("Per-bin PC1 shift\n(off-diagonal = reorganization; "
                  "sign flip = switch)", fontweight="bold", pad=8)
    axR.legend(frameon=False, fontsize=8, loc="upper left")
    for s in ("top", "right"): axR.spines[s].set_visible(False)

    # Annotate switch counts (sign flips) for context.
    n_BtoA = int(((p < 0) & (b > 0)).sum())   # was B, now A
    n_AtoB = int(((p > 0) & (b < 0)).sum())   # was A, now B
    axR.text(0.97, 0.03,
             f"B\u2192A: {n_BtoA:,}\nA\u2192B: {n_AtoB:,}",
             transform=axR.transAxes, ha="right", va="bottom", fontsize=8,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                       edgecolor="lightgrey", alpha=0.85))

    fig.suptitle("Hi-C compartment reorganization: Parental \u2192 Bottom10",
                 fontweight="bold", fontsize=12, y=1.02)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"pc1_compartment_structure.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Loading classified table (for heatmap)...")
    clf = pd.read_csv(args.classified, sep="\t", low_memory=False)
    if "switch_cat" not in clf.columns:
        print("[INFO] 'switch_cat' absent; classifying...")
        clf["switch_cat"] = classify_switch(clf)

    marks = args.marks or CANDIDATE_MARKS
    print("\n[INFO] === Figure 1: marks x switch-category heatmap ===")
    build_heatmap(clf, marks, args.signal_suffix, args.outdir)

    print("\n[INFO] === Figure 2: PC1 compartment-structure plot ===")
    pc1 = pd.read_csv(args.pc1, sep="\t", low_memory=False)
    build_pc1_structure(pc1, args.outdir)

    print("\n[INFO] Done. The heatmap is a VISUAL SUMMARY of the baseline "
          "violins (inferential backing stays with KW/Dunn/Cliff's delta). The "
          "PC1 structure plot is an honest saddle-plot substitute (a true saddle "
          "plot needs the Hi-C contact matrix, not just the eigenvector).")


if __name__ == "__main__":
    main()