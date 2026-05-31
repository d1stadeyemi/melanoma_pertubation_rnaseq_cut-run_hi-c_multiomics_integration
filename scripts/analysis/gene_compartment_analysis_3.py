"""
=============================================================================
gene_compartment_analysis.py  — STEP 3: expression change by switch direction
=============================================================================

GOAL
  Test the canonical structure->function relationship: do genes whose bin
  switched A->B (silenced) LOSE expression in Bottom10, and do genes whose bin
  switched B->A (activated) GAIN expression? This is the functional payoff of
  the compartment analysis — it shows switching is not just a structural
  rearrangement but tracks transcriptional output.

  Literature expectation (state in talk): compartment switches are associated
  with concordant expression changes — "genes which changed from A to B reduced
  their average expression levels, whereas genes which changed from B to A
  increased their average expression" (Dixon 2015; multiple Hi-C studies). The
  A/B compartment has a "contributory but not deterministic" role in expression,
  so expect a real but imperfect concordance, not a perfect one.

UNIT OF ANALYSIS = GENE (not bin)
  Steps 1-2 worked per-bin. Step 3 works per-GENE, because expression is a
  per-gene quantity. Each gene inherits the switch category of the bin its TSS
  falls in (via the gene list's 'coordinate'), then is joined to its DESeq2
  log2FC by gene symbol. Per-gene stats are valid here (expression is genuinely
  measured per gene — no pseudoreplication; that concern is specific to Step 4's
  bin-level mark values).

JOINS (two, each guarded)
  1. gene list  --coordinate-->  Step1 switch_cat   (each gene -> its bin's category)
  2. gene list  --gene_name-->   DESeq2 log2FC       (each gene -> its expression change)
  Gene-symbol joins are lossy (alias/case/version-suffix). The script reports
  the matched fraction and halts if implausibly low.

SIGN CONVENTION (verify on first run)
  DESeq2 log2FC must be oriented Bottom10-vs-Parental, the SAME direction as the
  switch calls, so that:
      A->B (silenced)  => expect log2FC < 0  (down in B10)
      B->A (activated) => expect log2FC > 0  (up in B10)
  If the medians come out with the OPPOSITE signs, your DESeq2 contrast was
  Parental-vs-B10; re-run with --flip_lfc to align. The script prints the
  observed median direction against the expectation so you can check.

INPUT
  --step1   Step1_bin_gene_density.txt   (bin coordinate + switch_cat)
  --gene_list  the hg19 gene list        (gene_name + coordinate)
  --deseq2  DESeq2 results               (gene_name + log2FC[/padj])

OUTPUTS (into --outdir)
  Step3_expression_by_switch.pdf/png   violin+box of log2FC by switch category
  Step3_expression_stats.txt           per-category n, median log2FC, MWU vs 0,
                                        and pairwise tests + Cliff's delta
  Step3_gene_expression_joined.txt     per-gene: gene, coordinate, switch_cat,
                                        log2FC, padj (auditable join result)

USAGE
  python gene_compartment_step3_expression.py \
      --step1     results/exploratory/gene_compartment/Step1_bin_gene_density.txt \
      --gene_list data/processed/GeneListTSShg19_cleaned.txt \
      --deseq2    data/processed/deseq2_results.txt \
      --outdir    results/exploratory/gene_compartment
"""

import argparse
import os
import sys
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal

try:
    import scikit_posthocs as sp
except ImportError:
    sp = None


# ── Aesthetics — consistent with Steps 1-2 and RF_SHAP ───────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    13,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "savefig.dpi":       600,
    "figure.facecolor":  "white",
})
SWITCH_COLORS = {"B_to_A": "#C0392B", "A_to_B": "#1A5276",
                 "Stable_A": "#85C1E9", "Stable_B": "#D5D8DC"}
SWITCH_LABELS = {"B_to_A": "B\u2192A (Activated)", "A_to_B": "A\u2192B (Silenced)",
                 "Stable_A": "Stable A", "Stable_B": "Stable B"}
CAT_ORDER = ["Stable_A", "A_to_B", "B_to_A", "Stable_B"]
# What each switching category should do to expression if structure tracks function.
EXPECTED_DIR = {"A_to_B": "negative (silenced in B10)",
                "B_to_A": "positive (activated in B10)"}


def cliffs_delta_fast(x, y):
    """Cliff's delta of x vs y. Positive => x tends to EXCEED y. (Correct sign.)"""
    x = np.asarray(x); y = np.asarray(y)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt = np.searchsorted(ys, x, side="left")
    le = np.searchsorted(ys, x, side="right")
    gt = ny - le
    d = (int(np.sum(lt)) - int(np.sum(gt))) / (nx * ny)
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


def _p_to_stars(p):
    if p < 1e-4: return "****"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"


def _bracket(ax, x1, x2, y, label, lw=1.0, fs=9):
    dy = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
    ax.plot([x1, x1, x2, x2], [y, y + dy, y + dy, y], lw=lw, c="black")
    ax.text((x1 + x2) / 2, y + dy * 1.2, label, ha="center", va="bottom", fontsize=fs)


def parse_args():
    p = argparse.ArgumentParser(
        description="Step 3: expression change (DESeq2) by switch direction.")
    p.add_argument("--step1", required=True,
                   help="Step1_bin_gene_density.txt (bin coordinate + switch_cat).")
    p.add_argument("--gene_list", required=True,
                   help="hg19 gene list (gene_name + coordinate).")
    p.add_argument("--deseq2", required=True,
                   help="DESeq2 results table (gene_name + log2FC).")
    p.add_argument("--outdir", default="results/exploratory/gene_compartment")
    p.add_argument("--coord_col", default="coordinate")
    p.add_argument("--gene_name_col", default="gene_name",
                   help="Gene-symbol column shared by gene_list and DESeq2.")
    p.add_argument("--lfc_col", default="log2FC",
                   help="DESeq2 log2 fold-change column [default %default].")
    p.add_argument("--padj_col", default="padj",
                   help="DESeq2 adjusted-p column (optional) [default %default].")
    p.add_argument("--flip_lfc", action="store_true",
                   help="Negate DESeq2 log2FC (use if its contrast is "
                        "Parental-vs-B10 instead of B10-vs-Parental).")
    p.add_argument("--min_match", type=float, default=0.70,
                   help="Min fraction of DESeq2 symbols that must match the gene "
                        "list, else halt (probable ID-format mismatch). Default 0.70.")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cc, gn = args.coord_col, args.gene_name_col

    # ── Load the three inputs ────────────────────────────────────────────────
    print("[INFO] Loading Step 1 (switch categories per bin)...")
    step1 = pd.read_csv(args.step1, sep="\t", low_memory=False)
    for c in (cc, "switch_cat"):
        if c not in step1.columns:
            sys.exit(f"[ERROR] '{c}' not in {args.step1}.")

    print("[INFO] Loading gene list...")
    genes = pd.read_csv(args.gene_list, sep="\t", low_memory=False)
    for c in (cc, gn):
        if c not in genes.columns:
            sys.exit(f"[ERROR] '{c}' not in gene list. Columns: {list(genes.columns)}")

    print("[INFO] Loading DESeq2 results...")
    de = pd.read_csv(args.deseq2, sep="\t", low_memory=False)
    if gn not in de.columns:
        sys.exit(f"[ERROR] '{gn}' not in DESeq2. Columns: {list(de.columns)}")
    if args.lfc_col not in de.columns:
        sys.exit(f"[ERROR] log2FC column '{args.lfc_col}' not in DESeq2. "
                 f"Columns: {list(de.columns)}")

    # ── Normalize join keys (symbols: strip + uppercase for matching only) ───
    genes[gn] = genes[gn].astype(str).str.strip()
    de[gn] = de[gn].astype(str).str.strip()
    genes["_key"] = genes[gn].str.upper()
    de["_key"] = de[gn].str.upper()

    # ── Join 1: gene -> its bin's switch category (via coordinate) ───────────
    genes[cc] = genes[cc].astype(str).str.strip()
    step1[cc] = step1[cc].astype(str).str.strip()
    g = genes.merge(step1[[cc, "switch_cat"]], on=cc, how="inner")
    g = g[g["switch_cat"].isin(CAT_ORDER)].copy()
    print(f"[INFO] Genes assigned to a classified bin: {len(g):,}")

    # ── Join 2: gene -> DESeq2 log2FC (via symbol), with match guard ─────────
    de_lfc = de[["_key", args.lfc_col]].copy()
    de_lfc[args.lfc_col] = pd.to_numeric(de_lfc[args.lfc_col], errors="coerce")
    if args.padj_col in de.columns:
        de_lfc[args.padj_col] = pd.to_numeric(de[args.padj_col], errors="coerce")
    de_lfc = de_lfc.dropna(subset=[args.lfc_col]).drop_duplicates("_key")

    de_keys = set(de_lfc["_key"])
    gene_keys = set(g["_key"])
    matched = de_keys & gene_keys
    frac = len(matched) / max(len(de_keys), 1)
    print(f"[INFO] DESeq2 symbols matching the gene list: "
          f"{len(matched):,}/{len(de_keys):,} ({frac:.1%})")
    if frac < args.min_match:
        sys.exit(
            f"[ERROR] Only {frac:.1%} of DESeq2 symbols matched the gene list "
            f"(below --min_match {args.min_match:.0%}). Likely an ID-format "
            f"mismatch (symbol vs Ensembl/Entrez, or version suffixes).\n"
            f"  Example DESeq2 keys: {list(de_keys)[:3]}\n"
            f"  Example gene keys:   {list(gene_keys)[:3]}\n"
            f"  Make both use the same ID type (gene symbols), then re-run."
        )

    g = g.merge(de_lfc, on="_key", how="inner")
    g["lfc"] = g[args.lfc_col]
    if args.flip_lfc:
        g["lfc"] = -g["lfc"]
        print("[INFO] --flip_lfc applied: log2FC negated to align with switch "
              "direction (B10 vs Parental).")
    g = g.dropna(subset=["lfc"])
    print(f"[INFO] Genes with both a switch category and expression: {len(g):,}")

    # ── Save the auditable joined per-gene table ─────────────────────────────
    keep = [gn, cc, "switch_cat", "lfc"]
    if args.padj_col in g.columns:
        keep.append(args.padj_col)
    keep = [c for c in keep if c in g.columns]
    g[keep].to_csv(os.path.join(args.outdir, "Step3_gene_expression_joined.txt"),
                   sep="\t", index=False)

    cats = [c for c in CAT_ORDER if c in g["switch_cat"].unique()]
    data = {c: g.loc[g["switch_cat"] == c, "lfc"].values for c in cats}

    # ── Per-category stats: median log2FC + Mann-Whitney vs 0 (sign test) ────
    # "Is this category's expression change significantly different from no
    # change?" Wilcoxon signed-rank vs 0 would also work; we report median +
    # a one-sample-style test via MWU against a zero-vector proxy is not ideal,
    # so use scipy's Wilcoxon signed-rank against 0 directly.
    from scipy.stats import wilcoxon
    rows = []
    for c in cats:
        vals = data[c]
        med = float(np.median(vals))
        # Wilcoxon signed-rank vs 0 (drop exact zeros as the test requires).
        nz = vals[vals != 0]
        if len(nz) > 10:
            try:
                w_stat, w_p = wilcoxon(nz)
            except ValueError:
                w_stat, w_p = np.nan, np.nan
        else:
            w_stat, w_p = np.nan, np.nan
        exp = EXPECTED_DIR.get(c, "n/a (stable)")
        rows.append({"switch_cat": c, "n_genes": len(vals),
                     "median_log2FC": med,
                     "wilcoxon_vs0_p": w_p, "expected_direction": exp})
    per_cat = pd.DataFrame(rows)
    print("\n[INFO] Expression change per category (median log2FC, B10/Parental):")
    print(per_cat.to_string(index=False))

    # ── Direction sanity check against literature expectation ────────────────
    def med(c): return float(np.median(data[c])) if c in data else np.nan
    print("\n[CHECK] Direction vs expectation:")
    if "A_to_B" in data:
        print(f"  A->B (silenced): median log2FC = {med('A_to_B'):+.3f} "
              f"(expect NEGATIVE)  -> "
              f"{'OK' if med('A_to_B') < 0 else 'UNEXPECTED (consider --flip_lfc)'}")
    if "B_to_A" in data:
        print(f"  B->A (activated): median log2FC = {med('B_to_A'):+.3f} "
              f"(expect POSITIVE)  -> "
              f"{'OK' if med('B_to_A') > 0 else 'UNEXPECTED (consider --flip_lfc)'}")

    # ── Pairwise: Kruskal-Wallis omnibus + Dunn + Cliff's delta ──────────────
    kw_H, kw_p = kruskal(*[data[c] for c in cats])
    pair_rows = []
    if sp is not None:
        dunn = sp.posthoc_dunn(g[["switch_cat", "lfc"]], val_col="lfc",
                               group_col="switch_cat", p_adjust="fdr_bh")
        for c1, c2 in combinations(cats, 2):
            p_pair = float(dunn.loc[c1, c2])
            d_pair, mag = cliffs_delta_fast(data[c1], data[c2])
            pair_rows.append({"group1": c1, "group2": c2, "dunn_p_BH": p_pair,
                              "cliffs_delta_g1_vs_g2": d_pair,
                              "effect_magnitude": mag, "stars": _p_to_stars(p_pair)})

    # ── Write combined stats file ────────────────────────────────────────────
    with open(os.path.join(args.outdir, "Step3_expression_stats.txt"), "w") as fh:
        fh.write("# Per-category expression change\n")
        per_cat.to_csv(fh, sep="\t", index=False)
        fh.write(f"\n# Kruskal-Wallis omnibus\nH\t{kw_H:.4f}\tdf\t{len(cats)-1}"
                 f"\tp\t{kw_p:.3e}\n")
        if pair_rows:
            fh.write("\n# Pairwise Dunn (BH) + Cliff's delta\n")
            pd.DataFrame(pair_rows).to_csv(fh, sep="\t", index=False)
    print(f"\n[INFO] Kruskal-Wallis: H={kw_H:.1f}, df={len(cats)-1}, p={kw_p:.3e}")
    if pair_rows:
        print(pd.DataFrame(pair_rows).to_string(index=False))

    # ── Figure: violin+box of log2FC by category, zero line for reference ────
    plot_data = [data[c] for c in cats]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.axhline(0, color="#888888", lw=0.8, linestyle="--", zorder=0)
    parts = ax.violinplot(plot_data, showmedians=False, showextrema=False,
                          widths=0.75)
    for pc, c in zip(parts["bodies"], cats):
        pc.set_facecolor(SWITCH_COLORS[c]); pc.set_alpha(0.75); pc.set_edgecolor("none")
    ax.boxplot(plot_data, widths=0.18, showfliers=False, patch_artist=True,
               medianprops=dict(color="black", linewidth=1.8),
               boxprops=dict(facecolor="white", alpha=0.9, linewidth=0.8),
               whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8))
    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels([SWITCH_LABELS[c] for c in cats], rotation=20, ha="right")
    ax.set_ylabel("Gene expression log2FC (Bottom10 / Parental)", labelpad=6)
    ax.set_title("Expression change across compartment-switch categories\n"
                 f"Kruskal-Wallis H = {kw_H:.1f}, "
                 + ("p < 1e-300" if kw_p == 0 or kw_p < 1e-300 else f"p = {kw_p:.1e}"),
                 fontweight="bold", pad=10)

    # Significance + effect-size brackets (same convention as Step 2).
    if pair_rows:
        cat_pos = {c: i + 1 for i, c in enumerate(cats)}
        sig = [r for r in pair_rows if r["dunn_p_BH"] < 0.05]
        sig.sort(key=lambda r: abs(cat_pos[r["group1"]] - cat_pos[r["group2"]]))
        top = max(np.max(d) for d in plot_data)
        tick = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.07
        base_y = top + tick * 0.6
        for k, r in enumerate(sig):
            label = f"{r['stars']}  \u03b4={r['cliffs_delta_g1_vs_g2']:+.2f} ({r['effect_magnitude']})"
            _bracket(ax, cat_pos[r["group1"]], cat_pos[r["group2"]],
                     base_y + k * tick, label)
        ax.set_ylim(ax.get_ylim()[0], base_y + (len(sig) + 1) * tick)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(args.outdir, f"Step3_expression_by_switch.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)

    print("\n[INFO] Step 3 complete. Verify the [CHECK] directions above match "
          "expectation (A->B down, B->A up). If they're reversed, re-run with "
          "--flip_lfc. Then proceed to Step 4 (histone-mark change overlay).")


if __name__ == "__main__":
    main()