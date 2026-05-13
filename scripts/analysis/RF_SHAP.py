#!/usr/bin/env python3
"""
Epigenetic priming predicts compartment switching — Random Forest + SHAP

Usage:
    python scripts/RF_SHAP.py \
        outputs/tables/cut_run_signal_table.txt \
        outputs/tables/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir outputs/figures/Q1_RF_SHAP

Outputs:
    Fig1a_switch_distribution.pdf/png
    Fig1b_violin_marks_by_switch.pdf/png
    Fig1c_RF_ROC.pdf/png
    Fig1d_SHAP_importance.pdf/png
    Fig1e_PCA_marks.pdf/png
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
import matplotlib.patches as mpatches
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import shap

# ── Publication aesthetics ────────────────────────────────────────────────────
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

# ── Color palettes ────────────────────────────────────────────────────────────
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


def parse_args():
    p = argparse.ArgumentParser(
        description="Q1: Epigenetic priming predicts compartment switching"
    )
    p.add_argument("cutrun_file",  help="cut_run_signal_table.txt")
    p.add_argument("merged_file",  help="merged multi-omics table")
    p.add_argument("--outdir",     default="outputs/figures/Q1_RF_SHAP")
    return p.parse_args()


def save_fig(fig, outdir, stem):
    """Save as PDF (vector submission) and PNG (preview)."""
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


def find_col(df, must, must_not=None):
    """
    Find first column matching all keywords in must
    and none in must_not. Prefers group mean columns
    (no replicate suffix like _r1, _r2).
    """
    must_not = must_not or []
    cands = []
    for col in df.columns:
        c = col.lower()
        if all(k in c for k in must) and not any(k in c for k in must_not):
            cands.append(col)
    # Prefer group mean: no _r<digit> suffix
    pref = [c for c in cands
            if not any(f"_r{i}" in c.lower() for i in range(1, 10))]
    return pref[0] if pref else (cands[0] if cands else None)


def classify_switch(merged_df):
    """
    Assign compartment switch category per gene.
    Priority: Both_original > Switch_Orig > Shift_Orig > infer from PC1 > Unknown.
    Returns four categories: B_to_A, A_to_B, Stable_A, Stable_B.
    """
    cats = []
    for _, row in merged_df.iterrows():
        cat = None
        # Priority 1: Both (switch + shift)
        v = str(row.get("P->B10_Both_original", "")).strip()
        if v not in ["-", "", "nan", "None"]:
            if "toA" in v or "BtoA" in v:
                cat = "B_to_A"
            elif "toB" in v or "AtoB" in v:
                cat = "A_to_B"
        # Priority 2: Switch only
        if cat is None:
            v = str(row.get("P->B10_Switch_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if v == "BtoA":
                    cat = "B_to_A"
                elif v == "AtoB":
                    cat = "A_to_B"
        # Priority 3: Shift only
        if cat is None:
            v = str(row.get("P->B10_Shift_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if "TowardA" in v:
                    cat = "B_to_A"
                elif "TowardB" in v:
                    cat = "A_to_B"
        # Priority 4: Infer stable from PC1 sign
        if cat is None:
            p_orig = row.get("Parental_orig", np.nan)
            b_orig = row.get("Bottom10_orig", np.nan)
            try:
                p_orig = float(p_orig)
                b_orig = float(b_orig)
                if p_orig > 0 and b_orig > 0:
                    cat = "Stable_A"
                elif p_orig < 0 and b_orig < 0:
                    cat = "Stable_B"
            except (ValueError, TypeError):
                pass
        cats.append(cat)
    return cats


def load_data(cutrun_file, merged_file):
    """Load both files, detect mark columns, classify switches, merge."""
    print("[INFO] Loading CUT&RUN signal table...")
    cutrun = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    print(f"  {len(cutrun):,} rows x {len(cutrun.columns)} columns")

    print("[INFO] Loading merged multi-omics table...")
    merged = pd.read_csv(merged_file, sep="\t", low_memory=False)
    print(f"  {len(merged):,} rows x {len(merged.columns)} columns")

    # Detect Parental mark mean columns from CUT&RUN table
    mark_map = {}
    for mark, keywords in [
        ("H3K27me3", ["parental", "h3k27me3", "mean"]),
        ("H3K27Ac",  ["parental", "h3k27ac",  "mean"]),
        ("H3K4me1",  ["parental", "h3k4me1",  "mean"]),
        ("H3K9me3",  ["parental", "h3k9me3",  "mean"]),
    ]:
        col = find_col(cutrun, keywords, ["minus", "unc", "ezh", "b10", "bottom"])
        if col:
            mark_map[mark] = col
            print(f"  Detected {mark}: {col}")
        else:
            print(f"  [WARN] Could not detect {mark} column")

    if len(mark_map) < 2:
        sys.exit("[ERROR] Need at least 2 Parental mark columns.")

    # Build CUT&RUN features table
    feat_cols  = ["gene_name"] + list(mark_map.values())
    cutrun_sub = cutrun[feat_cols].dropna(subset=["gene_name"])
    cutrun_sub = cutrun_sub.drop_duplicates(subset="gene_name")

    # Convert to numeric
    for col in mark_map.values():
        cutrun_sub[col] = pd.to_numeric(cutrun_sub[col], errors="coerce")

    # Rename mark columns to clean names
    rename = {v: k for k, v in mark_map.items()}
    cutrun_sub = cutrun_sub.rename(columns=rename)

    # Classify switches from merged table
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

    # Merge features + switch categories
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

    print(f"[INFO] Final dataset: {len(df):,} genes with marks + switch call")
    return df, list(mark_map.keys())


def plot_switch_distribution(df, outdir):
    """
    Fig 1a: Bar chart of switch category counts.
    Establishes scale of genomic reorganization.
    """
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
    ax.set_title(
        "Hi-C Compartment Switch Categories\nParental \u2192 Bottom10",
        fontweight="bold", pad=10
    )
    ax.tick_params(axis="x", rotation=20)
    ax.set_ylim(0, max(counts) * 1.15)
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1a_switch_distribution")


def plot_violin_marks(df, mark_names, outdir):
    """
    Fig 1b: Violin plots of each Parental mark signal
    stratified by compartment switch category.
    Kruskal-Wallis p-value annotated per panel.
    """
    cat_order  = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    cat_order  = [c for c in cat_order if c in df["switch_cat"].unique()]
    short_labs = [SWITCH_LABELS[c] for c in cat_order]

    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()

    for ax, mark in zip(axes, mark_names):
        groups = [df[df["switch_cat"] == c][mark].dropna().values
                  for c in cat_order]
        valid_groups = [(lab, g) for lab, g in zip(short_labs, groups)
                        if len(g) > 1]
        if not valid_groups:
            ax.set_visible(False)
            continue
        labs_v, data_v = zip(*valid_groups)

        parts = ax.violinplot(
            data_v, showmedians=True, showextrema=False, widths=0.7
        )
        # Color each violin by its category
        for pc, cat in zip(parts["bodies"],
                           [c for c in cat_order
                            if len(df[df["switch_cat"] == c][mark].dropna()) > 1]):
            pc.set_facecolor(SWITCH_COLORS[cat])
            pc.set_alpha(0.75)
            pc.set_edgecolor("none")
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.5)

        # Kruskal-Wallis p-value
        if len(data_v) >= 2:
            _, p = kruskal(*data_v)
            ax.text(0.97, 0.97, f"p={p:.2e}",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=8, color="black")

        ax.set_xticks(range(1, len(labs_v) + 1))
        ax.set_xticklabels(labs_v, fontsize=8, rotation=30, ha="right")
        ax.set_title(mark, fontweight="bold", fontsize=10)
        ax.set_ylabel("CPM-normalized signal", labelpad=4)

    fig.suptitle(
        "Parental Epigenetic Marks at Compartment-Switching vs Stable Loci",
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_by_switch")


def run_rf_shap(df, mark_names, outdir):
    """
    Random Forest classification: Parental marks -> predict switcher vs stable.
    5-fold stratified CV, 50 repeats estimated via multiple seeds.
    SHAP TreeExplainer for feature importance.
    ROC curve with mean AUC ± SD.
    """
    # Binary: Switcher (B_to_A or A_to_B) vs Stable
    df_m = df[df["switch_cat"].isin(
        ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    )].copy()
    df_m["is_switcher"] = df_m["switch_cat"].isin(["B_to_A", "A_to_B"]).astype(int)

    X = df_m[mark_names].values
    y = df_m["is_switcher"].values

    if len(X) < 50 or y.sum() < 5:
        print("[WARN] Insufficient data for RF. Skipping.")
        return

    print(f"[INFO] RF dataset: {len(X)} genes, "
          f"{y.sum()} switchers, {(y==0).sum()} stable")

    rf = RandomForestClassifier(
        n_estimators=500, max_depth=6,
        class_weight="balanced", random_state=42, n_jobs=-1
    )
    cv    = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    tprs  = []
    aucs  = []
    mfpr  = np.linspace(0, 1, 100)

    # ── ROC plot ──────────────────────────────────────────────────────────────
    fig_roc, ax_roc = plt.subplots(figsize=(5, 5))

    for train_idx, test_idx in cv.split(X, y):
        rf.fit(X[train_idx], y[train_idx])
        y_prob = rf.predict_proba(X[test_idx])[:, 1]
        fpr, tpr, _ = roc_curve(y[test_idx], y_prob)
        fold_auc = auc(fpr, tpr)
        aucs.append(fold_auc)
        tprs.append(np.interp(mfpr, fpr, tpr))
        ax_roc.plot(fpr, tpr, color="#85C1E9", alpha=0.3, lw=1)

    mean_tpr = np.mean(tprs, axis=0)
    mean_auc = np.mean(aucs)
    std_auc  = np.std(aucs)

    ax_roc.plot(mfpr, mean_tpr, color="#C0392B", lw=2,
                label=f"Mean ROC (AUC = {mean_auc:.3f} \u00b1 {std_auc:.3f})")
    ax_roc.fill_between(
        mfpr,
        mean_tpr - np.std(tprs, axis=0),
        mean_tpr + np.std(tprs, axis=0),
        color="#E74C3C", alpha=0.12
    )
    ax_roc.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
    ax_roc.set_xlabel("False Positive Rate", labelpad=6)
    ax_roc.set_ylabel("True Positive Rate", labelpad=6)
    ax_roc.set_title(
        "Random Forest: Parental Marks\nPredict Compartment Switching (5-fold CV)",
        fontweight="bold", pad=10
    )
    ax_roc.legend(frameon=False, fontsize=9, loc="lower right")
    save_fig(fig_roc, outdir, "Fig1c_RF_ROC")
    print(f"[INFO] Cross-validated AUC: {mean_auc:.3f} ± {std_auc:.3f}")

    # ── SHAP on full dataset ───────────────────────────────────────────────────
    rf.fit(X, y)
    explainer = shap.TreeExplainer(rf)
    shap_vals = explainer.shap_values(X)
    # Handle all SHAP return shapes:
    #   - list of [class0, class1] arrays (older SHAP)         -> pick [1]
    #   - 3D array (n_samples, n_features, n_classes) (newer)  -> pick [..., 1]
    #   - 2D array (n_samples, n_features) (regression-like)   -> use as-is
    if isinstance(shap_vals, list):
        sv = shap_vals[1]
    elif np.ndim(shap_vals) == 3:
        sv = shap_vals[:, :, 1]
    else:
        sv = shap_vals

    mean_shap = np.abs(sv).mean(axis=0).ravel()   # force 1D
    order     = np.argsort(mean_shap)[::-1]       # 1D int array

    fig_shap, ax_shap = plt.subplots(figsize=(5.5, 4))
    top_idx    = int(np.argmax(mean_shap))
    bar_colors = ["#C0392B" if int(i) == top_idx else "#1A5276" for i in order]
    ax_shap.barh(
        [mark_names[i] for i in order[::-1]],
        [mean_shap[i] for i in order[::-1]],
        color=bar_colors[::-1], edgecolor="none", height=0.6
    )
    ax_shap.set_xlabel("Mean |SHAP value|\n(contribution to switcher prediction)",
                       labelpad=6)
    ax_shap.set_title(
        "Feature Importance: Parental Marks\n\u2192 Compartment Switch Prediction (SHAP)",
        fontweight="bold", pad=10
    )
    save_fig(fig_shap, outdir, "Fig1d_SHAP_importance")


def plot_pca(df, mark_names, outdir):
    """
    Fig 1e: PCA of Parental epigenetic mark space colored by switch outcome.
    Shows whether switchers separate from stable loci in feature space.
    """
    cat_order = ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]
    df_p = df[df["switch_cat"].isin(cat_order)].dropna(subset=mark_names)

    scaler = StandardScaler()
    Xsc    = scaler.fit_transform(df_p[mark_names].values)
    pca    = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xsc)
    var    = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(5.5, 5))
    for cat in cat_order:
        mask = df_p["switch_cat"].values == cat
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=SWITCH_COLORS[cat],
            s=10 if "Stable" in cat else 25,
            alpha=0.25 if "Stable" in cat else 0.80,
            linewidths=0, rasterized=True,
            label=SWITCH_LABELS[cat]
        )
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)", labelpad=6)
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)", labelpad=6)
    ax.set_title(
        "PCA: Parental Epigenetic Mark Space\nby Compartment Switch Outcome",
        fontweight="bold", pad=10
    )
    ax.legend(frameon=False, fontsize=9, markerscale=2, loc="best")
    save_fig(fig, outdir, "Fig1e_PCA_marks")


def save_statistics(df, mark_names, outdir):
    """Mann-Whitney U: switchers vs stable per mark. Save as supplementary table."""
    results = []
    switchers_mask = df["switch_cat"].isin(["B_to_A", "A_to_B"])
    stable_mask    = df["switch_cat"].isin(["Stable_A", "Stable_B"])
    for mark in mark_names:
        sw = df[switchers_mask][mark].dropna().values
        st = df[stable_mask][mark].dropna().values
        if len(sw) > 0 and len(st) > 0:
            stat, p = mannwhitneyu(sw, st, alternative="two-sided")
            results.append({
                "mark":             mark,
                "n_switchers":      len(sw),
                "n_stable":         len(st),
                "median_switchers": np.median(sw),
                "median_stable":    np.median(st),
                "mannwhitney_U":    stat,
                "p_value":          p,
            })
    res_df = pd.DataFrame(results)
    if len(res_df):
        res_df["padj_BH"] = multipletests(res_df["p_value"], method="fdr_bh")[1]
        out = os.path.join(outdir, "Statistics.txt")
        res_df.to_csv(out, sep="\t", index=False)
        print(f"[INFO] Statistics saved -> {out}")
        print(res_df[["mark", "median_switchers", "median_stable",
                       "p_value", "padj_BH"]].to_string())


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df, mark_names = load_data(args.cutrun_file, args.merged_file)

    # Save classified gene table for downstream R scripts
    out_genes = os.path.join(args.outdir, "Switch_classified.txt")
    df.to_csv(out_genes, sep="\t", index=False)
    print(f"[INFO] Classified gene table saved -> {out_genes}")

    print("\n[INFO] Generating Fig 1a...")
    plot_switch_distribution(df, args.outdir)

    print("[INFO] Generating Fig 1b...")
    plot_violin_marks(df, mark_names, args.outdir)

    print("[INFO] Running Random Forest + SHAP...")
    run_rf_shap(df, mark_names, args.outdir)

    print("[INFO] Generating Fig 1e: PCA...")
    plot_pca(df, mark_names, args.outdir)

    print("[INFO] Computing statistics...")
    save_statistics(df, mark_names, args.outdir)

    print(f"\n[INFO] Q1 Python complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()