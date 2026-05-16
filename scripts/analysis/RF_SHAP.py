"""
Epigenetic priming predicts compartment switching — Random Forest + SHAP

Usage:
    python scripts/RF_SHAP.py \
        data/processed/cut_run_signal_table.txt \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir results/exploratory/RF_SHAP \
        --include_ctcf            # optional: add parental CTCF as a 5th feature
        --shap_dependence_partners H3K27me3:H3K27Ac H3K4me1:H3K9me3  # optional: pairs of marks for SHAP dependence plots

Outputs:
    Fig1a_switch_distribution.pdf/png
    Fig1b_violin_marks_by_switch.pdf/png   (log1p, colored medians, Cliff's delta + MWU)
    Fig1c_RF_ROC_PR.pdf/png                (ROC + PR-AUC side-by-side, per-fold dot inset)
    Fig1d_SHAP_importance.pdf/png          (bar + beeswarm, aligned row order)
    Fig1d_dependence_<mark>.pdf/png        (SHAP dependence plots for top features)
    Fig1e_PCA_marks.pdf/png                (log1p, scatter only — no contours)
    Statistics.txt
    Switch_classified.txt
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_curve, auc,
    average_precision_score, precision_recall_curve,
)
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import shap

# ── Aesthetics ────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "Liberation Sans",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    13,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    1.0,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "savefig.dpi":       600,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
})

SWITCH_COLORS = {
    "B_to_A":   "#C0392B",
    "A_to_B":   "#1A5276",
    "Stable_A": "#85C1E9",
    "Stable_B": "#D5D8DC",
}
SWITCH_LABELS = {
    "B_to_A":   "B\u2192A (Activated)",
    "A_to_B":   "A\u2192B (Silenced)",
    "Stable_A": "Stable A",
    "Stable_B": "Stable B",
}


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Q1: Epigenetic priming predicts compartment switching"
    )
    p.add_argument("cutrun_file",  help="cut_run_signal_table.txt")
    p.add_argument("merged_file",  help="merged multi-omics table")
    p.add_argument("--outdir",     default="results/exploratory/RF_SHAP",
                   help="Output directory for figures and tables")
    p.add_argument("--include_ctcf", action="store_true",
                   help="Include parental_ctcf_mean as a 5th feature if present.")
    p.add_argument(
        "--shap_dependence_partners",
        nargs="+",
        default=None,
        metavar="MARK1:MARK2",
        help=(
            "Pairs of marks for SHAP dependence plots, formatted MARK1:MARK2. "
            "For each pair, two panels are generated side-by-side: "
            "MARK1 colored by MARK2, and MARK2 colored by MARK1. "
            "Provide multiple pairs separated by spaces. "
            "Example: --shap_dependence_partners H3K27me3:H3K27Ac H3K4me1:H3K9me3 "
            "If omitted, no dependence plots are generated."
        ),
    )    
    return p.parse_args()


def save_fig(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


# =============================================================================
# Helpers
# =============================================================================

def find_col(df, must, must_not=None):
    must_not = must_not or []
    cands = []
    for col in df.columns:
        c = col.lower()
        if all(k in c for k in must) and not any(k in c for k in must_not):
            cands.append(col)
    pref = [c for c in cands
            if not any(f"_r{i}" in c.lower() for i in range(1, 10))]
    return pref[0] if pref else (cands[0] if cands else None)


def cliffs_delta_fast(x, y):
    """
    Memory-efficient Cliff's delta via sort + searchsorted.
    Convention: delta(x, y) = P(x > y) - P(x < y) in {-1, +1}.
    Magnitude thresholds: Romano et al. 2006.
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt_per_x = np.searchsorted(ys, x, side="left")
    le_per_x = np.searchsorted(ys, x, side="right")
    gt_per_x = ny - le_per_x
    total_gt = int(np.sum(gt_per_x))
    total_lt = int(np.sum(lt_per_x))
    # x>y count per element = (ny - le); total_gt is P(x>y)*nx*ny
    d  = (total_gt - total_lt) / (nx * ny)
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


def classify_switch(merged_df):
    cats = []
    for _, row in merged_df.iterrows():
        cat = None
        v = str(row.get("P->B10_Both_original", "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if "toA" in v or "BtoA" in v:
                cat = "B_to_A"
            elif "toB" in v or "AtoB" in v:
                cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Switch_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if v == "BtoA":   cat = "B_to_A"
                elif v == "AtoB": cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Shift_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if "TowardA" in v:   cat = "B_to_A"
                elif "TowardB" in v: cat = "A_to_B"
        if cat is None:
            p_orig = row.get("Parental_orig", np.nan)
            b_orig = row.get("Bottom10_orig", np.nan)
            try:
                p_orig = float(p_orig); b_orig = float(b_orig)
                if   p_orig > 0 and b_orig > 0: cat = "Stable_A"
                elif p_orig < 0 and b_orig < 0: cat = "Stable_B"
            except (ValueError, TypeError):
                pass
        cats.append(cat)
    return cats


def load_data(cutrun_file, merged_file, include_ctcf=False):
    print("[INFO] Loading CUT&RUN signal table...")
    cutrun = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    print(f"  {len(cutrun):,} rows x {len(cutrun.columns)} columns")

    print("[INFO] Loading merged multi-omics table...")
    merged = pd.read_csv(merged_file, sep="\t", low_memory=False)
    print(f"  {len(merged):,} rows x {len(merged.columns)} columns")

    mark_spec = [
        ("H3K27me3", ["parental", "h3k27me3", "mean"]),
        ("H3K27Ac",  ["parental", "h3k27ac",  "mean"]),
        ("H3K4me1",  ["parental", "h3k4me1",  "mean"]),
        ("H3K9me3",  ["parental", "h3k9me3",  "mean"]),
    ]
    if include_ctcf:
        mark_spec.append(("CTCF", ["parental", "ctcf", "mean"]))

    mark_map = {}
    for mark, keywords in mark_spec:
        col = find_col(cutrun, keywords, ["minus", "unc", "ezh", "b10", "bottom"])
        if col:
            mark_map[mark] = col
            print(f"  Detected {mark}: {col}")
        else:
            if mark == "CTCF" and include_ctcf:
                print(f"  [WARN] --include_ctcf set but no parental CTCF column "
                      f"found; continuing without CTCF.")
            else:
                print(f"  [WARN] Could not detect {mark} column")

    if len(mark_map) < 2:
        sys.exit("[ERROR] Need at least 2 Parental mark columns.")

    feat_cols  = ["gene_name"] + list(mark_map.values())
    cutrun_sub = cutrun[feat_cols].dropna(subset=["gene_name"])
    cutrun_sub = cutrun_sub.drop_duplicates(subset="gene_name")
    for col in mark_map.values():
        cutrun_sub[col] = pd.to_numeric(cutrun_sub[col], errors="coerce")
    rename = {v: k for k, v in mark_map.items()}
    cutrun_sub = cutrun_sub.rename(columns=rename)

    print("[INFO] Classifying compartment switches...")
    switch_cols = [
        "gene_name", "Parental_orig", "Bottom10_orig",
        "P->B10_Switch_Orig", "P->B10_Shift_Orig",
        "P->B10_Both_original", "B10-P_Orig"
    ]
    avail = [c for c in switch_cols if c in merged.columns]
    merged_sw = merged[avail].drop_duplicates(subset="gene_name").copy()
    merged_sw["switch_cat"] = classify_switch(merged_sw)
    print(f"[INFO] Switch category counts:")
    print(merged_sw["switch_cat"].value_counts(dropna=False).to_string())

    df = cutrun_sub.merge(
        merged_sw[["gene_name", "switch_cat"] +
                  [c for c in avail if c != "gene_name"]],
        on="gene_name", how="inner"
    )
    df = df[df["switch_cat"].notna()].copy()
    df = df[df["switch_cat"].isin(["B_to_A", "A_to_B",
                                    "Stable_A", "Stable_B"])]

    mark_names = list(mark_map.keys())
    df[mark_names] = df[mark_names].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=mark_names)

    print(f"[INFO] Final dataset: {len(df):,} genes x {len(mark_names)} marks")
    return df, mark_names


# =============================================================================
# Figures
# =============================================================================

def plot_switch_distribution(df, outdir):
    """Fig 1a — bar chart of switch category counts (unchanged)."""
    cats   = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    cats   = [c for c in cats if c in df["switch_cat"].unique()]
    counts = [df["switch_cat"].value_counts().get(c, 0) for c in cats]
    colors = [SWITCH_COLORS[c] for c in cats]
    labels = [SWITCH_LABELS[c] for c in cats]

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bars = ax.bar(labels, counts, color=colors,
                  edgecolor="white", linewidth=0.5, width=0.6)
    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.012,
                f"{val:,}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Number of genes", labelpad=6)
    ax.set_title("Hi-C Compartment Switch Categories\nParental \u2192 Bottom10",
                 fontweight="bold", pad=10)
    ax.tick_params(axis="x", rotation=20)
    ax.set_ylim(0, max(counts) * 1.15)
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1a_switch_distribution")


def plot_violin_marks(df, mark_names, outdir):
    """
    Fig 1b — log1p violins per mark by switch category.
    Colored crossbar medians + Cliff's delta + MWU (Switchers vs Stable).
    """
    cat_order  = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    cat_order  = [c for c in cat_order if c in df["switch_cat"].unique()]
    n_marks    = len(mark_names)
    n_cols     = 3
    n_rows     = (n_marks + n_cols - 1) // n_cols   # ceil division

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    sw_mask = df["switch_cat"].isin(["B_to_A", "A_to_B"])
    st_mask = df["switch_cat"].isin(["Stable_A", "Stable_B"])

    for ax, mark in zip(axes, mark_names):
        data_per_cat = []
        keep_cats    = []
        for c in cat_order:
            vals = df.loc[df["switch_cat"] == c, mark].dropna().values
            vals = np.log1p(vals)
            if len(vals) > 1:
                data_per_cat.append(vals)
                keep_cats.append(c)
        if not data_per_cat:
            ax.set_visible(False); continue

        parts = ax.violinplot(data_per_cat, showmedians=False,
                              showextrema=False, widths=0.7)
        for pc, cat in zip(parts["bodies"], keep_cats):
            pc.set_facecolor(SWITCH_COLORS[cat])
            pc.set_alpha(0.75)
            pc.set_edgecolor("none")

        for i, (vals, cat) in enumerate(zip(data_per_cat, keep_cats), start=1):
            med = np.median(vals)
            ax.hlines(med, i - 0.22, i + 0.22,
                      colors="black", linewidth=2.2, zorder=3)
            ax.hlines(med, i - 0.18, i + 0.18,
                      colors=SWITCH_COLORS[cat], linewidth=1.4, zorder=4)

        sw = df.loc[sw_mask, mark].dropna().values
        st = df.loc[st_mask, mark].dropna().values
        if len(sw) > 0 and len(st) > 0:
            _, p_mwu = mannwhitneyu(sw, st, alternative="two-sided")
            d, mag   = cliffs_delta_fast(sw, st)
            p_str = "p < 1e-300" if (p_mwu == 0 or p_mwu < 1e-300) \
                    else f"p = {p_mwu:.1e}"
            ax.text(
                0.97, 0.97,
                f"Switchers vs Stable\nCliff's \u03b4 = {d:+.3f} ({mag})\n{p_str}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="black",
                bbox=dict(boxstyle="round,pad=0.25",
                          facecolor="white", edgecolor="lightgrey",
                          linewidth=0.5, alpha=0.85)
            )

        ax.set_xticks(range(1, len(keep_cats) + 1))
        ax.set_xticklabels([SWITCH_LABELS[c] for c in keep_cats],
                           fontsize=8, rotation=30, ha="right")
        ax.set_title(mark, fontweight="bold", fontsize=10)
        ax.set_ylabel("log1p(CPM signal)", labelpad=4)

    for ax in axes[len(mark_names):]:
        ax.set_visible(False)

    fig.suptitle(
        "Parental Epigenetic Marks at Compartment-Switching vs Stable Loci",
        fontweight="bold", fontsize=12, y=1.00
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_by_switch")


def run_rf_shap(df, mark_names, outdir, dependence_pairs=None):
    """
    Fig 1c — ROC + PR side-by-side, random/baseline reference lines,
             per-fold AUC dot-plot inset.
    Fig 1d — SHAP bar + beeswarm (aligned row ordering).
    Fig 1d dependence — paired SHAP dependence plots, ONLY if the user
                        provided --shap_dependence_partners pairs.
    """
    df_m = df[df["switch_cat"].isin(
        ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    )].copy()
    df_m["is_switcher"] = df_m["switch_cat"].isin(
        ["B_to_A", "A_to_B"]
    ).astype(int)

    X = df_m[mark_names].values
    y = df_m["is_switcher"].values

    if len(X) < 50 or y.sum() < 5:
        print("[WARN] Insufficient data for RF. Skipping.")
        return

    prevalence = y.mean()
    print(f"[INFO] RF dataset: {len(X)} genes, "
          f"{y.sum()} switchers, {(y==0).sum()} stable")
    print(f"[INFO] Positive-class prevalence: {prevalence:.3f} "
          f"(= random PR baseline)")

    rf = RandomForestClassifier(
        n_estimators=200, max_depth=6,
        class_weight="balanced", random_state=42, n_jobs=-1
    )
    cv    = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    tprs  = []
    aucs  = []
    prcs  = []   # (recall_grid, mean_precision_per_fold) for plotting
    aps   = []   # average precision per fold (PR-AUC)
    mfpr  = np.linspace(0, 1, 100)
    mrec  = np.linspace(0, 1, 100)

    # ── Fig 1c: ROC + PR side by side ─────────────────────────────────────────
    fig_perf, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 5))

    for train_idx, test_idx in cv.split(X, y):
        rf.fit(X[train_idx], y[train_idx])
        y_prob = rf.predict_proba(X[test_idx])[:, 1]

        # ROC
        fpr, tpr, _ = roc_curve(y[test_idx], y_prob)
        fold_auc = auc(fpr, tpr)
        aucs.append(fold_auc)
        tprs.append(np.interp(mfpr, fpr, tpr))
        ax_roc.plot(fpr, tpr, color="#85C1E9", alpha=0.15, lw=1)

        # PR
        precision, recall, _ = precision_recall_curve(y[test_idx], y_prob)
        ap = average_precision_score(y[test_idx], y_prob)
        aps.append(ap)
        # Interpolate precision onto a common recall grid (descending recall ->
        # reverse arrays before np.interp)
        prcs.append(np.interp(mrec, recall[::-1], precision[::-1]))
        ax_pr.plot(recall, precision, color="#85C1E9", alpha=0.15, lw=1)

    # ROC mean + band + random
    mean_tpr = np.mean(tprs, axis=0)
    mean_auc = np.mean(aucs); std_auc = np.std(aucs)
    ax_roc.plot(mfpr, mean_tpr, color="#C0392B", lw=2,
                label=f"Mean ROC (AUC = {mean_auc:.3f} \u00b1 {std_auc:.3f})")
    ax_roc.fill_between(
        mfpr, mean_tpr - np.std(tprs, axis=0),
        mean_tpr + np.std(tprs, axis=0),
        color="#E74C3C", alpha=0.12
    )
    ax_roc.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5,
                label="Random (AUC = 0.5)")
    ax_roc.set_xlabel("False Positive Rate", labelpad=6)
    ax_roc.set_ylabel("True Positive Rate", labelpad=6)
    ax_roc.set_title("ROC Curve (5-fold CV)", fontweight="bold", pad=8)
    ax_roc.legend(frameon=False, fontsize=9, loc="lower right")

    # Per-fold AUC dot inset
    ax_in = ax_roc.inset_axes([0.10, 0.62, 0.25, 0.22])
    ax_in.scatter(aucs, [1] * len(aucs),
                  color="#C0392B", s=35, zorder=3,
                  edgecolor="white", linewidth=0.6)
    ax_in.axvline(mean_auc, color="black", lw=1.0, linestyle="--", zorder=2)
    ax_in.set_title("Per-fold AUC", fontsize=8, pad=3)
    ax_in.tick_params(axis="x", labelsize=7)
    ax_in.set_yticks([])
    ax_in.set_xlim(min(aucs) - 0.005, max(aucs) + 0.005)
    ax_in.set_ylim(0.5, 1.5)
    for s in ("top", "right", "left"):
        ax_in.spines[s].set_visible(False)

    # PR mean + band + prevalence baseline
    mean_prec = np.mean(prcs, axis=0)
    mean_ap   = np.mean(aps); std_ap = np.std(aps)
    ax_pr.plot(mrec, mean_prec, color="#C0392B", lw=2,
               label=f"Mean PR (PR-AUC = {mean_ap:.3f} \u00b1 {std_ap:.3f})")
    ax_pr.fill_between(
        mrec, mean_prec - np.std(prcs, axis=0),
        mean_prec + np.std(prcs, axis=0),
        color="#E74C3C", alpha=0.12
    )
    ax_pr.axhline(prevalence, color="black", lw=0.8, linestyle="--", alpha=0.5,
                  label=f"Random (PR-AUC = {prevalence:.3f})")
    ax_pr.set_xlabel("Recall", labelpad=6)
    ax_pr.set_ylabel("Precision", labelpad=6)
    ax_pr.set_title("Precision-Recall Curve (5-fold CV)",
                    fontweight="bold", pad=8)
    ax_pr.set_xlim(0, 1); ax_pr.set_ylim(0, 1)
    ax_pr.legend(frameon=False, fontsize=9, loc="upper right")

    fig_perf.suptitle(
        "Random Forest: Parental Marks Predict Compartment Switching",
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig_perf, outdir, "Fig1c_RF_ROC_PR")
    print(f"[INFO] Cross-validated ROC-AUC: {mean_auc:.3f} \u00b1 {std_auc:.3f}")
    print(f"[INFO] Cross-validated PR-AUC : {mean_ap:.3f} \u00b1 {std_ap:.3f}")
    print(f"[INFO] PR-AUC lift over random baseline: "
          f"{mean_ap / prevalence:.2f}x")
    print(f"[INFO] Per-fold ROC-AUC: {[f'{a:.3f}' for a in aucs]}")
    print(f"[INFO] Per-fold PR-AUC : {[f'{a:.3f}' for a in aps]}")

    # ── SHAP on full dataset ──────────────────────────────────────────────────
    rf.fit(X, y)
    explainer = shap.TreeExplainer(rf)
    shap_vals = explainer.shap_values(X)
    if isinstance(shap_vals, list):
        sv = shap_vals[1]
    elif np.ndim(shap_vals) == 3:
        sv = shap_vals[:, :, 1]
    else:
        sv = shap_vals

    mean_shap = np.abs(sv).mean(axis=0).ravel()
    order     = np.argsort(mean_shap)[::-1]       # descending importance indices
    top_idx   = int(np.argmax(mean_shap))

    # Fig 1d: bar + beeswarm — both with most-important on TOP
    fig_shap, (ax_bar, ax_bee) = plt.subplots(
        1, 2, figsize=(11, 4.5),
        gridspec_kw={"width_ratios": [1.0, 1.2]}
    )

    # Bar: most-important on top => reverse for barh which plots first at bottom
    bar_colors = ["#C0392B" if int(i) == top_idx else "#1A5276" for i in order]
    bars = ax_bar.barh(
        [mark_names[i] for i in order[::-1]],
        [mean_shap[i] for i in order[::-1]],
        color=bar_colors[::-1], edgecolor="none", height=0.6
    )
    for bar, val in zip(bars, [mean_shap[i] for i in order[::-1]]):
        ax_bar.text(bar.get_width() + max(mean_shap) * 0.012,
                    bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}",
                    va="center", ha="left", fontsize=9)
    ax_bar.set_xlabel("Mean |SHAP value|\n(magnitude of contribution)", labelpad=6)
    ax_bar.set_title("Feature Importance (magnitude)",
                     fontweight="bold", pad=8)
    ax_bar.set_xlim(0, max(mean_shap) * 1.18)

    # Beeswarm: same visual ordering (top of plot = most important).
    # shap.summary_plot draws first column at the BOTTOM => reverse the
    # descending-importance order so the top mark appears at the top.
    feat_order = [mark_names[i] for i in order][::-1]
    X_df       = pd.DataFrame(X, columns=mark_names)[feat_order]
    sv_df      = pd.DataFrame(sv, columns=mark_names)[feat_order]

    plt.sca(ax_bee)
    shap.summary_plot(
        sv_df.values,
        X_df.values,
        feature_names=feat_order,
        plot_type="dot",
        show=False,
        color_bar=True,
        sort=False,
        plot_size=None,
    )
    ax_bee.set_xlabel("SHAP value (impact on switcher probability)",
                      labelpad=6)
    ax_bee.set_title("Feature Importance (direction)",
                     fontweight="bold", pad=8)

    fig_shap.suptitle(
        "Parental Marks \u2192 Compartment Switch Prediction (SHAP)",
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig_shap, outdir, "Fig1d_SHAP_importance")

    # ── SHAP dependence plots — only if user specified pairs ─────────────────
    # Each pair (MARK1, MARK2) produces a 1x2 figure with bidirectional views.
    if dependence_pairs:
        plot_shap_dependence(sv, X, mark_names, dependence_pairs, outdir)
    else:
        print("[INFO] No --shap_dependence_partners specified; skipping "
              "dependence plots.")


def plot_shap_dependence(sv, X, mark_names, pairs, outdir):
    """
    SHAP dependence plots for explicitly specified mark pairs.

    For each (MARK1, MARK2) pair, produce a 1x2 figure:
        Left panel  : SHAP(MARK1) vs log1p(MARK1), colored by log1p(MARK2)
        Right panel : SHAP(MARK2) vs log1p(MARK2), colored by log1p(MARK1)

    This lets the user inspect the bidirectional interaction between the two
    marks at the level of the model's per-gene contributions.

    Parameters
    ----------
    sv         : np.ndarray, shape (n_genes, n_features)
                 SHAP values for the positive (switcher) class.
    X          : np.ndarray, shape (n_genes, n_features)
                 Raw feature matrix (CPM signal, untransformed).
    mark_names : list of str
                 Column order of X (and of sv).
    pairs      : list of (str, str) tuples
                 Mark pairs to plot. Both names must be in mark_names.
    outdir     : str
                 Output directory.
    """
# Case-insensitive mark lookup: tolerates typos like 'H3K27ac' vs 'H3K27Ac'
    name_to_idx_ci = {m.lower(): i for i, m in enumerate(mark_names)}

    def resolve(name):
        """Return canonical mark name for a user-supplied string, or None."""
        idx = name_to_idx_ci.get(name.lower())
        if idx is None:
            return None
        return mark_names[idx]

    for pair in pairs:
        feat_a_in, feat_b_in = pair
        feat_a = resolve(feat_a_in)
        feat_b = resolve(feat_b_in)

        if feat_a is None:
            print(f"  [WARN] '{feat_a_in}' not in detected marks "
                  f"{mark_names}; skipping pair {feat_a_in}:{feat_b_in}")
            continue
        if feat_b is None:
            print(f"  [WARN] '{feat_b_in}' not in detected marks "
                  f"{mark_names}; skipping pair {feat_a_in}:{feat_b_in}")
            continue

        idx_a = name_to_idx_ci[feat_a.lower()]
        idx_b = name_to_idx_ci[feat_b.lower()]
        print(f"[INFO] SHAP dependence pair: {feat_a} <-> {feat_b}")

        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4.8))

        # ── Left panel: SHAP(feat_a) vs feat_a, colored by feat_b ────────────
        sc_l = ax_l.scatter(
            np.log1p(X[:, idx_a]), sv[:, idx_a],
            c=np.log1p(X[:, idx_b]),
            cmap="plasma", s=8, alpha=0.55,
            edgecolors="none", rasterized=True,
        )
        ax_l.axhline(0, color="black", lw=0.6, alpha=0.5)
        ax_l.set_xlabel(f"log1p({feat_a}) (Parental CPM signal)", labelpad=6)
        ax_l.set_ylabel(
            f"SHAP value for {feat_a}\n(impact on switcher probability)",
            labelpad=6,
        )
        ax_l.set_title(
            f"SHAP Dependence: {feat_a}\n(colored by {feat_b})",
            fontweight="bold", pad=8,
        )
        cbar_l = plt.colorbar(sc_l, ax=ax_l, pad=0.02)
        cbar_l.set_label(f"log1p({feat_b})", fontsize=9)
        cbar_l.ax.tick_params(labelsize=8)

        # ── Right panel: SHAP(feat_b) vs feat_b, colored by feat_a ───────────
        sc_r = ax_r.scatter(
            np.log1p(X[:, idx_b]), sv[:, idx_b],
            c=np.log1p(X[:, idx_a]),
            cmap="plasma", s=8, alpha=0.55,
            edgecolors="none", rasterized=True,
        )
        ax_r.axhline(0, color="black", lw=0.6, alpha=0.5)
        ax_r.set_xlabel(f"log1p({feat_b}) (Parental CPM signal)", labelpad=6)
        ax_r.set_ylabel(
            f"SHAP value for {feat_b}\n(impact on switcher probability)",
            labelpad=6,
        )
        ax_r.set_title(
            f"SHAP Dependence: {feat_b}\n(colored by {feat_a})",
            fontweight="bold", pad=8,
        )
        cbar_r = plt.colorbar(sc_r, ax=ax_r, pad=0.02)
        cbar_r.set_label(f"log1p({feat_a})", fontsize=9)
        cbar_r.ax.tick_params(labelsize=8)

        plt.tight_layout()
        # Safe filename: e.g. "Fig1d_dependence_H3K27me3_vs_H3K27Ac"
        stem = f"Fig1d_dependence_{feat_a}_vs_{feat_b}"
        save_fig(fig, outdir, stem)


def parse_dependence_pairs(raw_pairs):
    """
    Parse CLI strings like 'H3K27me3:H3K27Ac' into (str, str) tuples.
    Returns an empty list if raw_pairs is None or empty.
    Skips malformed entries with a warning.
    """
    if not raw_pairs:
        return []
    out = []
    for s in raw_pairs:
        if ":" not in s:
            print(f"[WARN] --shap_dependence_partners entry '{s}' missing ':' "
                  f"(expected MARK1:MARK2); skipping.")
            continue
        parts = s.split(":")
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            print(f"[WARN] --shap_dependence_partners entry '{s}' malformed; "
                  f"skipping.")
            continue
        out.append((parts[0].strip(), parts[1].strip()))
    return out


def plot_pca(df, mark_names, outdir):
    """
    Fig 1e — PCA of log1p(Parental marks), scatter only (no contours).
    Stable categories: tiny, very faint background.
    Switchers: visible foreground.
    """
    cat_order = ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]
    df_p = df[df["switch_cat"].isin(cat_order)].dropna(subset=mark_names)

    scaler = StandardScaler()
    Xsc    = scaler.fit_transform(np.log1p(df_p[mark_names].values))
    pca    = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xsc)
    var    = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(6.2, 5.6))

    for cat in ["Stable_A", "Stable_B"]:
        mask = df_p["switch_cat"].values == cat
        if mask.sum() == 0: continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=SWITCH_COLORS[cat], s=4, alpha=0.06,
            linewidths=0, rasterized=True, label=SWITCH_LABELS[cat]
        )

    for cat in ["A_to_B", "B_to_A"]:
        mask = df_p["switch_cat"].values == cat
        if mask.sum() == 0: continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=SWITCH_COLORS[cat], s=18, alpha=0.55,
            linewidths=0, rasterized=True, label=SWITCH_LABELS[cat]
        )

    xlo, xhi = np.percentile(coords[:, 0], [1, 99])
    ylo, yhi = np.percentile(coords[:, 1], [1, 99])
    pad_x = 0.05 * (xhi - xlo); pad_y = 0.05 * (yhi - ylo)
    ax.set_xlim(xlo - pad_x, xhi + pad_x)
    ax.set_ylim(ylo - pad_y, yhi + pad_y)

    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)", labelpad=6)
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)", labelpad=6)
    ax.set_title(
        "PCA: Parental Epigenetic Mark Space (log1p)\nby Compartment Switch Outcome",
        fontweight="bold", pad=10
    )

    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=SWITCH_COLORS[c],
                   markersize=8, label=SWITCH_LABELS[c])
        for c in ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
        if c in df_p["switch_cat"].values
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=9, loc="best")

    plt.tight_layout()
    save_fig(fig, outdir, "Fig1e_PCA_marks")


def save_statistics(df, mark_names, outdir):
    """Mann-Whitney U + Cliff's delta: switchers vs stable per mark."""
    results = []
    sw_mask = df["switch_cat"].isin(["B_to_A", "A_to_B"])
    st_mask = df["switch_cat"].isin(["Stable_A", "Stable_B"])
    for mark in mark_names:
        sw = df[sw_mask][mark].dropna().values
        st = df[st_mask][mark].dropna().values
        if len(sw) > 0 and len(st) > 0:
            stat, p = mannwhitneyu(sw, st, alternative="two-sided")
            d, mag  = cliffs_delta_fast(sw, st)
            results.append({
                "mark":             mark,
                "n_switchers":      len(sw),
                "n_stable":         len(st),
                "median_switchers": np.median(sw),
                "median_stable":    np.median(st),
                "mannwhitney_U":    stat,
                "p_value":          p,
                "cliffs_delta":     d,
                "effect_magnitude": mag,
            })
    res_df = pd.DataFrame(results)
    if len(res_df):
        res_df["padj_BH"] = multipletests(res_df["p_value"], method="fdr_bh")[1]
        out = os.path.join(outdir, "Statistics.txt")
        res_df.to_csv(out, sep="\t", index=False)
        print(f"[INFO] Statistics saved -> {out}")
        print(res_df[["mark", "median_switchers", "median_stable",
                       "cliffs_delta", "effect_magnitude",
                       "p_value", "padj_BH"]].to_string())


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df, mark_names = load_data(args.cutrun_file, args.merged_file,
                                include_ctcf=args.include_ctcf)

    out_genes = os.path.join(args.outdir, "Switch_classified.txt")
    df.to_csv(out_genes, sep="\t", index=False)
    print(f"[INFO] Classified gene table saved -> {out_genes}")

    print("\n[INFO] Generating Fig 1a...")
    plot_switch_distribution(df, args.outdir)

    print("[INFO] Generating Fig 1b...")
    plot_violin_marks(df, mark_names, args.outdir)

    pairs = parse_dependence_pairs(args.shap_dependence_partners)
    if pairs:
        print(f"[INFO] SHAP dependence pairs requested: {pairs}")
    print("[INFO] Random Forest + SHAP (Fig 1c ROC+PR; Fig 1d bar+beeswarm"
          + (" + dependence)..." if pairs else ")..."))
    run_rf_shap(df, mark_names, args.outdir, dependence_pairs=pairs)

    print("[INFO] Generating Fig 1e: PCA (log1p, scatter only)...")
    plot_pca(df, mark_names, args.outdir)

    print("[INFO] Computing statistics...")
    save_statistics(df, mark_names, args.outdir)

    print(f"\n[INFO] Q1 Python complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()