"""
Epigenetic priming predicts compartment switching — Random Forest + SHAP

Usage:
    python scripts/RF_SHAP.py \
        data/processed/cut_run_signal_table.txt \
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \
        --outdir results/exploratory/RF_SHAP \
        --include_ctcf            # flag to include CTCF parental signal when available

Outputs:
    Fig1a_switch_distribution.pdf/png
    Fig1b_violin_marks_by_switch.pdf/png   (log1p Y, colored medians, Cliff's delta + MWU)
    Fig1c_RF_ROC.pdf/png                   (random baseline label + per-fold AUC inset)
    Fig1d_SHAP_importance.pdf/png          (bar + beeswarm)
    Fig1e_PCA_marks.pdf/png                (log1p inputs, KDE contours, faded stables)
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
import seaborn as sns
from scipy.stats import mannwhitneyu
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
    p.add_argument("--outdir",     default="results/exploratory/RF_SHAP",
                   help="Output directory for figures and tables")
    p.add_argument("--include_ctcf", action="store_true",
                   help="Include parental_ctcf_mean as a 5th feature if present.")
    return p.parse_args()


def save_fig(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


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


def cliffs_delta(x, y):
    """
    Cliff's delta: non-parametric effect size for ordinal/continuous data.
    Returns d in [-1, 1] and a qualitative magnitude string.

    Reference: Romano et al. 2006; standard in the genomics ML literature
    as the recommended companion to Mann-Whitney U.

    Magnitude thresholds (Romano et al. 2006):
        |d| < 0.147   negligible
        |d| < 0.33    small
        |d| < 0.474   medium
        else          large
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    # Vectorized count of x>y and x<y via sort + searchsorted
    y_sorted = np.sort(y)
    gt = np.sum(nx - np.searchsorted(y_sorted, x, side="right")) if False \
        else np.sum(x[:, None] > y[None, :])
    lt = np.sum(x[:, None] < y[None, :])
    d  = (gt - lt) / (nx * ny)
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


def cliffs_delta_fast(x, y):
    """
    Memory-efficient Cliff's delta via sort + searchsorted.
    For large n (15k stable vs 2.5k switchers, outer product blows up RAM).
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    # For each xi: count of y<xi (lt_xi) and y>xi (gt_xi)
    lt_per_x = np.searchsorted(ys, x, side="left")
    le_per_x = np.searchsorted(ys, x, side="right")
    gt_per_x = ny - le_per_x
    total_gt = int(np.sum(gt_per_x))
    total_lt = int(np.sum(lt_per_x))
    d  = (total_lt - total_gt) / (nx * ny)
    # NOTE: convention here is delta(x, y) = P(x>y) - P(x<y).
    # Using lt - gt above gives that sign because lt_per_x counts y<x.
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
                if v == "BtoA":
                    cat = "B_to_A"
                elif v == "AtoB":
                    cat = "A_to_B"
        if cat is None:
            v = str(row.get("P->B10_Shift_Orig", "")).strip()
            if v not in ["-", "", "nan", "None"]:
                if "TowardA" in v:
                    cat = "B_to_A"
                elif "TowardB" in v:
                    cat = "A_to_B"
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


def load_data(cutrun_file, merged_file, include_ctcf=False):
    """
    Load CUT&RUN and merged tables, detect Parental mark columns, classify
    compartment switches, merge into a per-gene feature+label table.

    If --include_ctcf is set and a parental CTCF mean column is detected,
    it will be added as a 5th feature. Otherwise the script proceeds with
    the four standard marks.
    """
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


def plot_switch_distribution(df, outdir):
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
    Fig 1b: log1p-transformed violins of each Parental mark by switch category.
    - Y-axis: log1p(CPM-normalized signal)  -> visible distributions, no outlier blowup
    - Median: COLORED crossbar per category (matches violin) instead of one black bar
    - Effect size annotation: Cliff's delta (Switchers vs Stable) with magnitude
    - Mann-Whitney U p-value (Switchers vs Stable) on same line
    """
    cat_order  = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    cat_order  = [c for c in cat_order if c in df["switch_cat"].unique()]
    short_labs = [SWITCH_LABELS[c] for c in cat_order]
    n_marks    = len(mark_names)
    n_rows     = (n_marks + 1) // 2

    fig, axes = plt.subplots(n_rows, 2, figsize=(8, 4 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    sw_mask = df["switch_cat"].isin(["B_to_A", "A_to_B"])
    st_mask = df["switch_cat"].isin(["Stable_A", "Stable_B"])

    for ax, mark in zip(axes, mark_names):
        # log1p transform for visualization
        data_per_cat = []
        keep_cats    = []
        for c in cat_order:
            vals = df.loc[df["switch_cat"] == c, mark].dropna().values
            vals = np.log1p(vals)
            if len(vals) > 1:
                data_per_cat.append(vals)
                keep_cats.append(c)
        if not data_per_cat:
            ax.set_visible(False)
            continue

        parts = ax.violinplot(
            data_per_cat, showmedians=False, showextrema=False, widths=0.7
        )
        for pc, cat in zip(parts["bodies"], keep_cats):
            pc.set_facecolor(SWITCH_COLORS[cat])
            pc.set_alpha(0.75)
            pc.set_edgecolor("none")

        # Colored median crossbars (one per violin, color-matched to category)
        for i, (vals, cat) in enumerate(zip(data_per_cat, keep_cats), start=1):
            med = np.median(vals)
            # Crossbar drawn as a short horizontal line in a contrasting shade
            ax.hlines(med, i - 0.22, i + 0.22,
                      colors="black", linewidth=2.2, zorder=3)
            ax.hlines(med, i - 0.18, i + 0.18,
                      colors=SWITCH_COLORS[cat], linewidth=1.4, zorder=4)

        # Switchers vs Stable: Mann-Whitney U + Cliff's delta
        sw = df.loc[sw_mask, mark].dropna().values
        st = df.loc[st_mask, mark].dropna().values
        if len(sw) > 0 and len(st) > 0:
            u_stat, p_mwu = mannwhitneyu(sw, st, alternative="two-sided")
            d, mag = cliffs_delta_fast(sw, st)
            # Format MWU p-value robustly (avoid e+00 underflow display)
            if p_mwu == 0 or p_mwu < 1e-300:
                p_str = "p < 1e-300"
            else:
                p_str = f"p = {p_mwu:.1e}"
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

    # Hide any leftover subplot if odd mark count
    for ax in axes[len(mark_names):]:
        ax.set_visible(False)

    fig.suptitle(
        "Parental Epigenetic Marks at Compartment-Switching vs Stable Loci",
        fontweight="bold", fontsize=12, y=1.00
    )
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_by_switch")


def run_rf_shap(df, mark_names, outdir):
    """
    Random Forest + SHAP.
    Fig 1c: ROC with random-chance label, per-fold curves at alpha 0.15,
            and inset histogram of per-fold AUCs.
    Fig 1d: SHAP bar (with numeric labels) AND SHAP beeswarm (directional).
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

    fig_roc, ax_roc = plt.subplots(figsize=(5.5, 5.2))

    for train_idx, test_idx in cv.split(X, y):
        rf.fit(X[train_idx], y[train_idx])
        y_prob = rf.predict_proba(X[test_idx])[:, 1]
        fpr, tpr, _ = roc_curve(y[test_idx], y_prob)
        fold_auc = auc(fpr, tpr)
        aucs.append(fold_auc)
        tprs.append(np.interp(mfpr, fpr, tpr))
        ax_roc.plot(fpr, tpr, color="#85C1E9", alpha=0.15, lw=1)

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
    ax_roc.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5,
                label="Random (AUC = 0.5)")
    ax_roc.set_xlabel("False Positive Rate", labelpad=6)
    ax_roc.set_ylabel("True Positive Rate", labelpad=6)
    ax_roc.set_title(
        "Random Forest: Parental Marks\nPredict Compartment Switching (5-fold CV)",
        fontweight="bold", pad=10
    )
    ax_roc.legend(frameon=False, fontsize=9, loc="lower right")

    # Inset histogram: per-fold AUCs
    ax_in = ax_roc.inset_axes([0.10, 0.62, 0.25, 0.22])
    ax_in.hist(aucs, bins=max(3, len(aucs)),
               color="#C0392B", edgecolor="white", linewidth=0.5)
    ax_in.set_title("Per-fold AUC", fontsize=8, pad=3)
    ax_in.tick_params(axis="both", labelsize=7)
    ax_in.set_xlim(min(aucs) - 0.01, max(aucs) + 0.01)
    for s in ("top", "right"):
        ax_in.spines[s].set_visible(False)

    save_fig(fig_roc, outdir, "Fig1c_RF_ROC")
    print(f"[INFO] Cross-validated AUC: {mean_auc:.3f} ± {std_auc:.3f}")
    print(f"[INFO] Per-fold AUCs: {[f'{a:.3f}' for a in aucs]}")

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
    order     = np.argsort(mean_shap)[::-1]
    top_idx   = int(np.argmax(mean_shap))

    # Fig 1d: bar + beeswarm side-by-side
    fig_shap, (ax_bar, ax_bee) = plt.subplots(
        1, 2, figsize=(11, 4.5),
        gridspec_kw={"width_ratios": [1.0, 1.2]}
    )

    # --- Bar (left) ---
    bar_colors = ["#C0392B" if int(i) == top_idx else "#1A5276" for i in order]
    bars = ax_bar.barh(
        [mark_names[i] for i in order[::-1]],
        [mean_shap[i] for i in order[::-1]],
        color=bar_colors[::-1], edgecolor="none", height=0.6
    )
    # Numeric labels at end of each bar
    for bar, val in zip(bars, [mean_shap[i] for i in order[::-1]]):
        ax_bar.text(bar.get_width() + max(mean_shap) * 0.012,
                    bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}",
                    va="center", ha="left", fontsize=9)
    ax_bar.set_xlabel("Mean |SHAP value|\n(magnitude of contribution)", labelpad=6)
    ax_bar.set_title("Feature Importance (magnitude)",
                     fontweight="bold", pad=8)
    ax_bar.set_xlim(0, max(mean_shap) * 1.18)

    # --- Beeswarm (right) ---
    # Order features the same way (top mark on top in beeswarm)
    feat_order = [mark_names[i] for i in order[::-1]]
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


def plot_pca(df, mark_names, outdir):
    """
    Fig 1e: PCA of log1p(Parental marks), colored by switch outcome.
    - Stable categories: tiny, very faint background (alpha ~0.06, s=4)
    - Switchers (B_to_A, A_to_B): visible foreground (alpha 0.55, s=18)
    - KDE contour overlays per category for visual structure cue
    """
    cat_order = ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]
    df_p = df[df["switch_cat"].isin(cat_order)].dropna(subset=mark_names)

    scaler = StandardScaler()
    Xsc    = scaler.fit_transform(np.log1p(df_p[mark_names].values))
    pca    = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xsc)
    var    = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(6.2, 5.6))

    # Background: stables (faint, small)
    for cat in ["Stable_A", "Stable_B"]:
        mask = df_p["switch_cat"].values == cat
        if mask.sum() == 0:
            continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=SWITCH_COLORS[cat], s=4, alpha=0.06,
            linewidths=0, rasterized=True, label=SWITCH_LABELS[cat]
        )

    # Foreground: switchers
    for cat in ["A_to_B", "B_to_A"]:
        mask = df_p["switch_cat"].values == cat
        if mask.sum() == 0:
            continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=SWITCH_COLORS[cat], s=18, alpha=0.55,
            linewidths=0, rasterized=True, label=SWITCH_LABELS[cat]
        )

    # KDE contour overlays per category
    for cat in cat_order:
        mask = df_p["switch_cat"].values == cat
        if mask.sum() < 50:
            continue
        try:
            sns.kdeplot(
                x=coords[mask, 0], y=coords[mask, 1],
                ax=ax,
                levels=4,
                color=SWITCH_COLORS[cat],
                linewidths=1.0,
                alpha=0.9,
            )
        except Exception as e:
            print(f"  [WARN] KDE failed for {cat}: {e}")

    # Zoom to bulk of data (1st-99th percentile of both axes)
    xlo, xhi = np.percentile(coords[:, 0], [1, 99])
    ylo, yhi = np.percentile(coords[:, 1], [1, 99])
    pad_x = 0.05 * (xhi - xlo)
    pad_y = 0.05 * (yhi - ylo)
    ax.set_xlim(xlo - pad_x, xhi + pad_x)
    ax.set_ylim(ylo - pad_y, yhi + pad_y)

    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)", labelpad=6)
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)", labelpad=6)
    ax.set_title(
        "PCA: Parental Epigenetic Mark Space (log1p)\nby Compartment Switch Outcome",
        fontweight="bold", pad=10
    )

    # Manual legend (forces visible markers regardless of plotted alpha)
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
    """
    Mann-Whitney U + Cliff's delta: switchers vs stable per mark.
    """
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

    print("[INFO] Running Random Forest + SHAP (Fig 1c + Fig 1d)...")
    run_rf_shap(df, mark_names, args.outdir)

    print("[INFO] Generating Fig 1e: PCA (log1p, faded stables, KDE)...")
    plot_pca(df, mark_names, args.outdir)

    print("[INFO] Computing statistics...")
    save_statistics(df, mark_names, args.outdir)

    print(f"\n[INFO] Q1 Python complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()