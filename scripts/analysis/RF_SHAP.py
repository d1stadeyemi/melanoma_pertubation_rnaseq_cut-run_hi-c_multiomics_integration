"""
Epigenetic priming predicts compartment switching — Random Forest / XGBoost + SHAP
(bin-keyed, multi-stat feature version)

This is the bin-based / 250kb-row adaptation of the original gene-based pipeline.
Each row of the prediction matrix is ONE 250kb genomic bin (one row per Hi-C
compartment label), and features are summary statistics of parental histone
marks within that bin (mean, max, std, coverage — any subset that the extract
script wrote).

WHY THIS DESIGN  ------------------------------------------------------------
  - One row per 250kb bin -> rows = labels -> standard StratifiedKFold is valid
    (no GroupKFold needed, no spatial-autocorrelation leakage between adjacent
    children of the same parent bin).
  - Multi-stat features (mean/max/std/coverage) preserve within-bin chromatin
    structure WITHOUT broadcasting one label across 50 5kb rows. This matches
    the Random Forest baseline of Zheng et al. (iScience 2024 "CoRNN"), which
    used mean+std per region and showed raw sub-bin features HURT the RF.

JOIN  -----------------------------------------------------------------------
  CUT&RUN and Hi-C are joined on the COORDINATE column (default: extra_1 in the
  CUT&RUN table, 'coordinate' in the Hi-C/merged file). Both must be at the
  same bin resolution (e.g. 250kb) or the join collapses. The script halts if
  the two key columns barely overlap, rather than silently producing all-NaN
  rows.

FEATURE COLUMNS  ------------------------------------------------------------
  Auto-detected from the CUT&RUN table. The detector picks any column of the
  form 'parental_<mark>_<stat>' for marks in --marks and stats in --stats.
  Defaults: marks = H3K27me3 H3K27Ac H3K4me1 H3K9me3 (+ CTCF if --include_ctcf);
  stats = mean max std coverage (filtered to whatever's present in the file).
  Treated-sample columns (containing 'unc'/'epz') and log2ratio columns are
  IGNORED — this ML uses only parental untreated baseline features.

OUTPUTS  --------------------------------------------------------------------
  All original figures/tables retained:
    Fig1a_switch_distribution_<mode>.pdf/png
    Fig1b_violin_marks_<mode>.pdf/png
    Fig1b_violin_pairwise_heatmap_<mode>.pdf/png  (all mode only)
    Fig1c_<model>_ROC_PR_<mode>.pdf/png
    Fig1d_SHAP_<model>_<mode>.pdf/png
    Fig1d_dependence_<pair>_<model>_<mode>.pdf/png
    Fig1d_dependence_<pair>_<model>_<mode>_zoomed.pdf/png
    Fig1e_PCA_marks_<mode>.pdf/png
    Statistics_<mode>.txt
    Statistics_dunn_<feature>_all.txt  (all mode only, per-feature)
    Switch_classified.txt              (always)
    SHAP_per_row_<model>_<mode>.txt    (was per_gene; now keyed on row_id)
  NEW:
    Feature_matrix_<mode>.txt          one file per comparison: row_id +
                                        all feature columns + switch label,
                                        for the exact rows the model scored.

USAGE -----------------------------------------------------------------------
    python scripts/RF_SHAP.py \\
        data/processed/cut_run_signal_table.txt \\
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \\
        --outdir results/exploratory/RF_SHAP \\
        --model rf \\
        --comparison both \\
        --include_ctcf

    # Restrict the feature set (e.g. just mean+std, matching CoRNN's RF baseline):
    python scripts/RF_SHAP.py ... --stats mean std

    # Override join keys if your columns are named differently:
    python scripts/RF_SHAP.py ... --cutrun_key extra_1 --hic_key coordinate
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
from matplotlib.colors import LinearSegmentedColormap
from itertools import combinations
from scipy.stats import mannwhitneyu, kruskal
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

try:
    from xgboost import XGBClassifier
    _HAS_XGB = True
except ImportError:
    XGBClassifier = None
    _HAS_XGB = False

try:
    import scikit_posthocs as sp
except ImportError:
    sp = None

# deep green = low feature value, deep yellow = high feature value
GREEN_YELLOW = LinearSegmentedColormap.from_list("green_yellow", ["#0AF16B", "#CEF3DD", "#F1EE09"])


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

COMPARISON_SPEC = {
    "A": {"neg": "Stable_A", "pos": "A_to_B",
          "label": "Stable_A vs A\u2192B", "short": "compA"},
    "B": {"neg": "Stable_B", "pos": "B_to_A",
          "label": "Stable_B vs B\u2192A", "short": "compB"},
}

DEFAULT_MARKS = ("H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3")
DEFAULT_STATS = ("mean", "max", "std", "coverage")
# Tokens that indicate a column belongs to a non-baseline track. Any column
# whose name contains one of these is rejected as a feature.
EXCLUDE_TOKENS = ("b10", "bottom", "unc", "ezh", "epz", "log2ratio", "minus")


# =============================================================================
# Model configuration
# =============================================================================
RF_PARAMS = dict(
    n_estimators=500,
    max_depth=12,
    class_weight="balanced_subsample",
    n_jobs=-1,
)
XGB_PARAMS = dict(
    n_estimators=500,
    max_depth=4,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_lambda=1.0,
    min_child_weight=1,
    tree_method="hist",
    n_jobs=-1,
    eval_metric="logloss",
)


def build_model(model_name, seed, scale_pos_weight=None):
    if model_name == "rf":
        return RandomForestClassifier(random_state=seed, **RF_PARAMS)
    elif model_name == "xgboost":
        if not _HAS_XGB:
            sys.exit("[ERROR] --model xgboost requested but xgboost is not installed.")
        params = dict(XGB_PARAMS)
        if scale_pos_weight is not None:
            params["scale_pos_weight"] = float(scale_pos_weight)
        return XGBClassifier(random_state=seed, **params)
    else:
        sys.exit(f"[ERROR] Unknown model '{model_name}'.")


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Epigenetic priming predicts compartment switching "
                    "(bin-keyed, multi-stat features).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("cutrun_file",  help="cut_run_signal_table.txt (250kb-bin)")
    p.add_argument("merged_file",  help="merged multi-omics table (Hi-C labels)")
    p.add_argument("--outdir",     default="results/exploratory/RF_SHAP",
                   help="Output directory")
    p.add_argument("--model", choices=["rf", "xgboost"], default="rf",
                   help="ML backend (default: rf)")

    # Join configuration ------------------------------------------------------
    p.add_argument("--cutrun_key", default=None,
                   help="Column in cutrun_file to join on. Default: auto-detect "
                        "first of [coordinate, coordinates, coord, coords] in "
                        "the file; if none, falls back to 'extra_1'.")
    p.add_argument("--hic_key", default=None,
                   help="Column in merged_file to join on. Default: auto-detect "
                        "first of [coordinate, coordinates, coord, coords].")
    p.add_argument("--row_id", default="region_id",
                   help="Column to use as the row identifier in outputs "
                        "(default: region_id). Falls back to the cutrun_key "
                        "value if region_id is absent.")
    p.add_argument("--min_overlap", type=float, default=0.10,
                   help="Minimum fraction of cutrun_key values that must be "
                        "found in hic_key for the join to proceed. Default "
                        "0.10 = 10%%. Below this the script halts with a "
                        "diagnostic (probable resolution mismatch).")

    # Feature configuration ---------------------------------------------------
    p.add_argument("--marks", nargs="+", default=list(DEFAULT_MARKS),
                   help="Histone marks to use as features (default: H3K27me3 "
                        "H3K27Ac H3K4me1 H3K9me3).")
    p.add_argument("--stats", nargs="+", default=None,
                   choices=list(DEFAULT_STATS),
                   help="Stat suffixes to use as features. Default: all stats "
                        "present in the file (intersection of {mean,max,std,"
                        "coverage} with what was actually extracted).")
    p.add_argument("--include_ctcf", action="store_true",
                   help="Include parental_ctcf_<stat> as additional features.")
    p.add_argument("--features", choices=["baseline", "dynamics", "all"],
                   default="baseline",
                   help=(
                       "Which feature SET to use for the ML:\n"
                       "  baseline (default) : parental_<mark>_<stat> columns "
                       "only. Answers 'does parental chromatin PRIME "
                       "switching?' — a chromatin-priming analysis.\n"
                       "  dynamics : log2ratio_p-b10_<mark>_<stat> columns "
                       "only (migration-axis log2 fold changes from the "
                       "extract script). Answers 'do migration-induced "
                       "CHANGES in marks predict switching?' — a chromatin-"
                       "dynamics analysis. Halts if no log2ratio columns "
                       "are found in the input.\n"
                       "  all : both baseline + dynamics features pooled. "
                       "Useful for assessing combined predictive power, "
                       "but interpret SHAP carefully — baseline and "
                       "dynamics features are correlated and will share "
                       "importance.\n"
                       "IMPORTANT FOR PUBLICATION: these answer DIFFERENT "
                       "biological questions. Report both 'baseline' and "
                       "'dynamics' analyses, not just whichever one gives a "
                       "more flattering result. Cherry-picking by feature "
                       "set is a form of analytical p-hacking."
                   ))

    # SHAP / comparison -------------------------------------------------------
    p.add_argument("--shap_dependence_partners", nargs="+", default=None,
                   metavar="FEAT1:FEAT2",
                   help="Pairs of FEATURE names for SHAP dependence plots, "
                        "e.g. H3K27me3_mean:H3K27Ac_mean. Use full "
                        "feature names as printed by the script.")
    p.add_argument("--comparison", choices=["A", "B", "both", "all"],
                   default="both",
                   help="A=Stable_A vs A->B, B=Stable_B vs B->A, both=A+B, "
                        "all=switchers-vs-stable (4-cat violin).")
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

def cliffs_delta_fast(x, y):
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "NA"
    ys = np.sort(y)
    lt_per_x = np.searchsorted(ys, x, side="left")
    le_per_x = np.searchsorted(ys, x, side="right")
    gt_per_x = ny - le_per_x
    d = (int(np.sum(lt_per_x)) - int(np.sum(gt_per_x))) / (nx * ny)
    ad = abs(d)
    if   ad < 0.147: mag = "negligible"
    elif ad < 0.330: mag = "small"
    elif ad < 0.474: mag = "medium"
    else:            mag = "large"
    return d, mag


def classify_switch(merged_df):
    """Classify each row into B_to_A / A_to_B / Stable_A / Stable_B / None.
    Reads switch columns in priority order, falls back to sign of original
    eigenvalues for stable categories. Identical logic to the gene-based
    version — operates row-wise, agnostic to the row identity."""
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


# -----------------------------------------------------------------------------
# Feature detection (baseline + dynamics)
# -----------------------------------------------------------------------------
def detect_baseline_features(cutrun_df, marks, stats):
    """
    Find PARENTAL BASELINE columns of the form 'parental_<mark>_<stat>'.
    These are the chromatin-priming features (does the baseline state predict
    switching?). Returns dict {display_name: source_column}, display names like
    'H3K27me3_mean'.
    """
    feat_map = {}
    cols_lc = {c: c.lower() for c in cutrun_df.columns}
    for mark in marks:
        for stat in stats:
            wanted = f"parental_{mark}_{stat}".lower()
            hits = []
            for col, lc in cols_lc.items():
                if lc == wanted and not any(tok in lc for tok in EXCLUDE_TOKENS
                                            if tok != "parental"):
                    # the EXCLUDE filter is moot here because we matched the
                    # exact baseline name, but kept for symmetry / safety.
                    hits.append(col)
            if hits:
                feat_map[f"{mark}_{stat}"] = hits[0]
                if len(hits) > 1:
                    print(f"  [WARN] Multiple columns match {wanted}; using "
                          f"{hits[0]}.")
    return feat_map


def detect_dynamics_features(cutrun_df, marks, stats):
    """
    Find MIGRATION LOG2 RATIO columns of the form 'log2ratio_p-b10_<mark>_<stat>'
    (matching the extract script's migration-axis naming). These are the
    chromatin-dynamics features (does migration-induced CHANGE in marks
    predict switching?).

    IMPORTANT: this filter is migration-specific. Drug-axis columns
    ('log2ratio_b10-b10_unc_*', 'log2ratio_p-p_unc_*') are EXCLUDED — those
    answer a different biological question (drug response, not migration)
    and would confound a switching predictor. Switching is a migration
    phenotype, so only migration log2 ratios make sense as features here.

    Display names use a 'd_' prefix to make them visually distinct from
    baseline features in shared plots (e.g. 'd_H3K27me3_mean' for the
    dynamics version vs 'H3K27me3_mean' for baseline). This prevents
    collisions if --features all combines both sets.
    """
    feat_map = {}
    cols_lc = {c: c.lower() for c in cutrun_df.columns}
    # Migration log2 naming from the extract script: log2ratio_p-b10_<mark>_<stat>
    # (parental abbreviated as 'p'; B10 verbatim). Match case-insensitively.
    for mark in marks:
        for stat in stats:
            wanted = f"log2ratio_p-b10_{mark}_{stat}".lower()
            hits = [col for col, lc in cols_lc.items() if lc == wanted]
            if hits:
                feat_map[f"d_{mark}_{stat}"] = hits[0]
                if len(hits) > 1:
                    print(f"  [WARN] Multiple columns match {wanted}; using "
                          f"{hits[0]}.")
    return feat_map


def auto_pick_join_key(df, override, file_label):
    """Return the join key column for one of the two files.
    If override is given and present in df, use it. Otherwise look for any of
    'coordinate','coordinates','coord','coords' (case-insensitive). If not found
    and this is the cutrun file, fall back to 'extra_1'. Otherwise error."""
    if override:
        if override in df.columns:
            return override
        sys.exit(f"[ERROR] {file_label}: --key '{override}' not in columns.")
    candidates = ("coordinate", "coordinates", "coord", "coords")
    cols_lc = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in cols_lc:
            return cols_lc[cand]
    if file_label == "cutrun_file" and "extra_1" in df.columns:
        return "extra_1"
    sys.exit(f"[ERROR] {file_label}: no coordinate-like column found and no "
             f"explicit --key given. Looked for: {candidates}.")


def load_data(cutrun_file, merged_file, marks, stats, include_ctcf,
              cutrun_key_override, hic_key_override, row_id_col,
              min_overlap_frac, feature_mode="baseline"):
    """
    Load both tables, auto-detect or use specified join keys, sanity-check
    overlap, build feature columns according to feature_mode, join, classify
    switches, return (df, feat_names, feat_map, row_id_col).

    feature_mode :
      'baseline' : parental_<mark>_<stat> columns only
      'dynamics' : log2ratio_p-b10_<mark>_<stat> columns only.
                   HALTS if no such columns are present (loud failure, per
                   user spec). This prevents silently running a 'dynamics'
                   analysis on a file that lacks log2 features.
      'all'      : both sets pooled. Dynamics features prefixed 'd_' to
                   avoid display collisions with baseline features.
    """
    print("[INFO] Loading CUT&RUN signal table...")
    cutrun = pd.read_csv(cutrun_file, sep="\t", low_memory=False)
    print(f"  {len(cutrun):,} rows x {len(cutrun.columns)} columns")

    print("[INFO] Loading merged multi-omics table...")
    merged = pd.read_csv(merged_file, sep="\t", low_memory=False)
    print(f"  {len(merged):,} rows x {len(merged.columns)} columns")

    # --- Join keys -----------------------------------------------------------
    ck = auto_pick_join_key(cutrun, cutrun_key_override, "cutrun_file")
    hk = auto_pick_join_key(merged, hic_key_override, "merged_file")
    print(f"[INFO] Join keys: cutrun.{ck}  <->  merged.{hk}")

    # --- Stats actually present in the file (intersect with what user asked) -
    detected_stats = []
    for s in stats:
        if any(c.lower().endswith(f"_{s.lower()}") for c in cutrun.columns):
            detected_stats.append(s)
    if not detected_stats:
        sys.exit(f"[ERROR] None of the requested stats {stats} present in "
                 f"{cutrun_file}. Re-check --stats and the extract output.")
    print(f"[INFO] Stats present in cutrun table (and requested): {detected_stats}")

    mark_list = list(marks)
    if include_ctcf:
        mark_list.append("CTCF")

    # --- Feature detection: branch on feature_mode --------------------------
    print(f"[INFO] Feature mode: {feature_mode}")
    feat_map = {}
    if feature_mode in ("baseline", "all"):
        baseline_map = detect_baseline_features(cutrun, mark_list, detected_stats)
        if feature_mode == "baseline" and not baseline_map:
            sys.exit("[ERROR] No baseline feature columns detected. Expected "
                     "names of the form 'parental_<mark>_<stat>'. Re-check "
                     "the extract script output and --marks/--stats.")
        feat_map.update(baseline_map)
        print(f"  [BASELINE] {len(baseline_map)} parental_<mark>_<stat> features")

    if feature_mode in ("dynamics", "all"):
        dynamics_map = detect_dynamics_features(cutrun, mark_list, detected_stats)
        if not dynamics_map:
            # Scream loudly per user spec — never silently fall through.
            sys.exit(
                "[ERROR] --features {} requested but NO migration log2 columns "
                "found in {}. Expected names of the form "
                "'log2ratio_p-b10_<mark>_<stat>'. Likely fix: re-run "
                "extract_bigwig_signal.py with --log2 to generate the "
                "migration-axis log2 ratios, then re-run this script.\n"
                "  Example of what's expected: log2ratio_p-b10_h3k27me3_mean".format(
                    feature_mode, cutrun_file)
            )
        feat_map.update(dynamics_map)
        print(f"  [DYNAMICS] {len(dynamics_map)} log2ratio_p-b10_<mark>_<stat> features")

    if not feat_map:
        sys.exit(f"[ERROR] No features detected for mode '{feature_mode}'.")
    print(f"[INFO] Total: {len(feat_map)} feature columns:")
    for disp, src in feat_map.items():
        print(f"  {disp:30s} <- {src}")

    # --- Overlap sanity check ------------------------------------------------
    ck_vals = set(cutrun[ck].dropna().astype(str))
    hk_vals = set(merged[hk].dropna().astype(str))
    common = ck_vals & hk_vals
    if len(ck_vals) == 0:
        sys.exit(f"[ERROR] cutrun_file column '{ck}' is empty.")
    overlap_frac = len(common) / len(ck_vals)
    print(f"[INFO] Key overlap: {len(common):,}/{len(ck_vals):,} cutrun rows "
          f"({overlap_frac:.1%}) match a merged row.")
    if overlap_frac < min_overlap_frac:
        sys.exit(
            f"[ERROR] Only {overlap_frac:.1%} of '{ck}' values are found in "
            f"'{hk}'. This is below --min_overlap {min_overlap_frac:.1%} and "
            f"almost certainly a resolution mismatch (e.g. 5kb cutrun vs 250kb "
            f"Hi-C). Re-extract CUT&RUN on the same bin grid as the Hi-C file.\n"
            f"  Example cutrun keys: {list(ck_vals)[:3]}\n"
            f"  Example hic keys:    {list(hk_vals)[:3]}"
        )

    # --- Build CUT&RUN subset (join key + features) --------------------------
    feat_src_cols = list(feat_map.values())
    keep_cols = [ck] + feat_src_cols
    if row_id_col in cutrun.columns and row_id_col != ck:
        keep_cols = [row_id_col] + keep_cols
    cutrun_sub = cutrun[keep_cols].copy()
    cutrun_sub = cutrun_sub.drop_duplicates(subset=ck)
    for col in feat_src_cols:
        cutrun_sub[col] = pd.to_numeric(cutrun_sub[col], errors="coerce")
    # Rename source columns to friendly display names.
    cutrun_sub = cutrun_sub.rename(columns={src: disp for disp, src in feat_map.items()})

    # --- Build merged subset (switch columns + join key) ---------------------
    switch_src_cols = [
        "Parental_orig", "Bottom10_orig",
        "P->B10_Switch_Orig", "P->B10_Shift_Orig",
        "P->B10_Both_original", "B10-P_Orig",
    ]
    avail_switch = [c for c in switch_src_cols if c in merged.columns]
    if not avail_switch:
        sys.exit(f"[ERROR] None of the switch columns {switch_src_cols} found "
                 f"in merged_file.")
    merged_sw = merged[[hk] + avail_switch].drop_duplicates(subset=hk).copy()
    merged_sw["switch_cat"] = classify_switch(merged_sw)
    print("[INFO] Switch category counts (full merged table):")
    print(merged_sw["switch_cat"].value_counts(dropna=False).to_string())

    # --- Join ---------------------------------------------------------------
    df = cutrun_sub.merge(merged_sw, left_on=ck, right_on=hk, how="inner")
    if ck != hk and hk in df.columns:
        df = df.drop(columns=[hk])

    df = df[df["switch_cat"].notna()].copy()
    df = df[df["switch_cat"].isin(["B_to_A", "A_to_B", "Stable_A", "Stable_B"])]

    feat_names = list(feat_map.keys())
    df[feat_names] = df[feat_names].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=feat_names)

    # Decide which row-id column to use for outputs.
    if row_id_col in df.columns:
        used_row_id = row_id_col
    else:
        print(f"  [WARN] --row_id '{row_id_col}' not found after merge; "
              f"falling back to join key '{ck}'.")
        used_row_id = ck

    print(f"[INFO] Final dataset: {len(df):,} bins x {len(feat_names)} features")
    return df, feat_names, feat_map, used_row_id


def subset_for_comparison(df, comp_key):
    spec = COMPARISON_SPEC[comp_key]
    sub = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])].copy()
    sub["is_switcher"] = (sub["switch_cat"] == spec["pos"]).astype(int)
    return sub, spec


# =============================================================================
# Fig 1a — switch distribution
# =============================================================================

def plot_switch_distribution(df, outdir, comparison):
    if comparison in ("both", "all"):
        cats = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]; title_sfx = ""
    elif comparison == "A":
        cats = ["Stable_A", "A_to_B"]; title_sfx = "  (A compartment only)"
    elif comparison == "B":
        cats = ["Stable_B", "B_to_A"]; title_sfx = "  (B compartment only)"
    else:
        cats = ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]; title_sfx = ""

    cats = [c for c in cats if c in df["switch_cat"].unique()]
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
    ax.set_ylabel("Number of 250kb bins", labelpad=6)
    ax.set_title("Hi-C Compartment Switch Categories\nParental \u2192 Bottom10"
                 + title_sfx,
                 fontweight="bold", pad=10)
    ax.tick_params(axis="x", rotation=20)
    ax.set_ylim(0, max(counts) * 1.15)
    plt.tight_layout()
    save_fig(fig, outdir, f"Fig1a_switch_distribution_{comparison}")


# =============================================================================
# Fig 1b — violins (operates on feature columns generically)
# =============================================================================

def _draw_violin_pair(ax, df, feat, neg_cat, pos_cat, title_extra=""):
    vals_neg = np.log1p(df.loc[df["switch_cat"] == neg_cat, feat].dropna().values)
    vals_pos = np.log1p(df.loc[df["switch_cat"] == pos_cat, feat].dropna().values)
    if len(vals_neg) < 2 or len(vals_pos) < 2:
        ax.set_visible(False); return

    parts = ax.violinplot([vals_neg, vals_pos],
                          showmedians=False, showextrema=False, widths=0.7)
    for pc, cat in zip(parts["bodies"], [neg_cat, pos_cat]):
        pc.set_facecolor(SWITCH_COLORS[cat]); pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    for i, (vals, cat) in enumerate(zip([vals_neg, vals_pos],
                                        [neg_cat, pos_cat]), start=1):
        med = np.median(vals)
        ax.hlines(med, i - 0.22, i + 0.22, colors="black", linewidth=2.2, zorder=3)
        ax.hlines(med, i - 0.18, i + 0.18, colors=SWITCH_COLORS[cat],
                  linewidth=1.4, zorder=4)

    _, p_mwu = mannwhitneyu(vals_pos, vals_neg, alternative="two-sided")
    d, mag = cliffs_delta_fast(vals_pos, vals_neg)
    p_str = "p < 1e-300" if (p_mwu == 0 or p_mwu < 1e-300) else f"p = {p_mwu:.1e}"
    ax.text(0.97, 0.97,
            f"{SWITCH_LABELS[pos_cat]} vs {SWITCH_LABELS[neg_cat]}\n"
            f"Cliff's \u03b4 = {d:+.3f} ({mag})\n{p_str}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="lightgrey", linewidth=0.5, alpha=0.85))
    ax.set_xticks([1, 2])
    ax.set_xticklabels([SWITCH_LABELS[neg_cat], SWITCH_LABELS[pos_cat]],
                       fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("log1p(signal)", labelpad=4)
    ax.set_title(feat + title_extra, fontweight="bold", fontsize=10)


def plot_violin_marks_AB(df, feat_names, outdir, comparison):
    if comparison in ("A", "B"):
        spec = COMPARISON_SPEC[comparison]
        n = len(feat_names)
        n_cols = min(4, n)
        n_rows = (n + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
        axes = np.atleast_2d(axes).flatten()
        for ax, feat in zip(axes, feat_names):
            _draw_violin_pair(ax, df, feat, spec["neg"], spec["pos"])
        for ax in axes[len(feat_names):]:
            ax.set_visible(False)
        fig.suptitle(f"Parental Features: {spec['label']}",
                     fontweight="bold", fontsize=12, y=1.00)
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1b_violin_marks_{comparison}")
        return

    # both mode: A and B side by side, one feature per row
    n = len(feat_names)
    fig, axes = plt.subplots(n, 2, figsize=(8, 4 * n))
    axes = np.atleast_2d(axes)
    for r_idx, feat in enumerate(feat_names):
        _draw_violin_pair(axes[r_idx, 0], df, feat,
                          COMPARISON_SPEC["A"]["neg"], COMPARISON_SPEC["A"]["pos"],
                          title_extra="\n(A compartment)")
        _draw_violin_pair(axes[r_idx, 1], df, feat,
                          COMPARISON_SPEC["B"]["neg"], COMPARISON_SPEC["B"]["pos"],
                          title_extra="\n(B compartment)")
    fig.suptitle("Parental Features: Compartment-Resolved Stable vs Switcher",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_both")


def _bracket(ax, x1, x2, y, label, lw=1.0, fs=8):
    dy = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
    ax.plot([x1, x1, x2, x2], [y, y + dy, y + dy, y], lw=lw, c="black")
    ax.text((x1 + x2) / 2, y + dy * 1.2, label,
            ha="center", va="bottom", fontsize=fs)


def _p_to_stars(p):
    if p < 1e-4:  return "****"
    if p < 1e-3:  return "***"
    if p < 1e-2:  return "**"
    if p < 5e-2:  return "*"
    return "ns"


def plot_violin_marks_all(df, feat_names, outdir):
    cat_order = [c for c in ["Stable_A", "Stable_B", "A_to_B", "B_to_A"]
                 if c in df["switch_cat"].unique()]
    n = len(feat_names)
    n_cols = min(4, n)
    n_rows = (n + n_cols - 1) // n_cols
    dunn_tables = {}

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 4.2 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    for ax, feat in zip(axes, feat_names):
        data_per_cat, keep_cats = [], []
        for c in cat_order:
            vals = df.loc[df["switch_cat"] == c, feat].dropna().values
            vals = np.log1p(vals)
            if len(vals) > 1:
                data_per_cat.append(vals); keep_cats.append(c)
        if not data_per_cat:
            ax.set_visible(False); continue

        parts = ax.violinplot(data_per_cat, showmedians=False,
                              showextrema=False, widths=0.7)
        for pc, cat in zip(parts["bodies"], keep_cats):
            pc.set_facecolor(SWITCH_COLORS[cat]); pc.set_alpha(0.75)
            pc.set_edgecolor("none")
        for i, (vals, cat) in enumerate(zip(data_per_cat, keep_cats), start=1):
            med = np.median(vals)
            ax.hlines(med, i - 0.22, i + 0.22, colors="black",
                      linewidth=2.2, zorder=3)
            ax.hlines(med, i - 0.18, i + 0.18, colors=SWITCH_COLORS[cat],
                      linewidth=1.4, zorder=4)

        if len(data_per_cat) >= 2:
            kw_stat, kw_p = kruskal(*data_per_cat)
        else:
            kw_stat, kw_p = np.nan, np.nan
        kw_str = ("p < 1e-300" if (kw_p == 0 or kw_p < 1e-300)
                  else f"p = {kw_p:.1e}")
        ax.text(0.03, 0.97,
                f"Kruskal-Wallis\nH = {kw_stat:.1f}\n{kw_str}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="lightgrey", linewidth=0.5, alpha=0.85))

        if sp is not None and len(data_per_cat) >= 2:
            long_df = pd.concat([pd.DataFrame({"val": d, "cat": c})
                                 for d, c in zip(data_per_cat, keep_cats)],
                                ignore_index=True)
            dunn = sp.posthoc_dunn(long_df, val_col="val", group_col="cat",
                                   p_adjust="fdr_bh")
            rows = []
            for c1, c2 in combinations(keep_cats, 2):
                p_pair = float(dunn.loc[c1, c2])
                d_pair, mag_pair = cliffs_delta_fast(
                    df.loc[df["switch_cat"] == c2, feat].dropna().values,
                    df.loc[df["switch_cat"] == c1, feat].dropna().values,
                )
                rows.append({"group1": c1, "group2": c2, "dunn_p_BH": p_pair,
                             "cliffs_delta_g2_vs_g1": d_pair,
                             "effect_magnitude": mag_pair,
                             "stars": _p_to_stars(p_pair)})
            dunn_tables[feat] = pd.DataFrame(rows)

            sig_pairs = [r for r in rows if r["dunn_p_BH"] < 0.05]
            cat_to_pos = {c: i + 1 for i, c in enumerate(keep_cats)}
            sig_pairs.sort(key=lambda r: abs(cat_to_pos[r["group1"]]
                                             - cat_to_pos[r["group2"]]))
            top = max(np.max(d) for d in data_per_cat)
            tick = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.07
            base_y = top + tick * 0.6
            for k, r in enumerate(sig_pairs):
                x1 = cat_to_pos[r["group1"]]; x2 = cat_to_pos[r["group2"]]
                _bracket(ax, x1, x2, base_y + k * tick, r["stars"])
            ax.set_ylim(ax.get_ylim()[0],
                        base_y + (len(sig_pairs) + 1) * tick)

        ax.set_xticks(range(1, len(keep_cats) + 1))
        ax.set_xticklabels([SWITCH_LABELS[c] for c in keep_cats],
                           fontsize=8, rotation=25, ha="right")
        ax.set_title(feat, fontweight="bold", fontsize=10)
        ax.set_ylabel("log1p(signal)", labelpad=4)

    for ax in axes[len(feat_names):]:
        ax.set_visible(False)

    fig.suptitle("Parental Features across Compartment-Switch Categories"
                 "\nKruskal-Wallis omnibus + Dunn's BH-adjusted post-hoc",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1b_violin_marks_all")

    for feat, tab in dunn_tables.items():
        out_tbl = os.path.join(outdir, f"Statistics_dunn_{feat}_all.txt")
        tab.to_csv(out_tbl, sep="\t", index=False)
        print(f"  Dunn table -> {out_tbl}")


# =============================================================================
# Fig 1c — ROC + PR
# =============================================================================

def _run_model_cv(X, y, model_name, seed=42):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs, aps, tprs, prcs = [], [], [], []
    mfpr = np.linspace(0, 1, 100); mrec = np.linspace(0, 1, 100)
    for tr_idx, te_idx in cv.split(X, y):
        y_tr = y[tr_idx]
        if model_name == "xgboost":
            n_pos = max(int(y_tr.sum()), 1); n_neg = int((y_tr == 0).sum())
            model = build_model(model_name, seed, scale_pos_weight=n_neg / n_pos)
        else:
            model = build_model(model_name, seed)
        model.fit(X[tr_idx], y_tr)
        y_prob = model.predict_proba(X[te_idx])[:, 1]
        fpr, tpr, _ = roc_curve(y[te_idx], y_prob)
        aucs.append(auc(fpr, tpr))
        tprs.append(np.interp(mfpr, fpr, tpr))
        precision, recall, _ = precision_recall_curve(y[te_idx], y_prob)
        aps.append(average_precision_score(y[te_idx], y_prob))
        prcs.append(np.interp(mrec, recall[::-1], precision[::-1]))
    if model_name == "xgboost":
        n_pos = max(int(y.sum()), 1); n_neg = int((y == 0).sum())
        model = build_model(model_name, seed, scale_pos_weight=n_neg / n_pos)
    else:
        model = build_model(model_name, seed)
    model.fit(X, y)
    return {
        "aucs": np.array(aucs), "aps": np.array(aps),
        "tprs": np.array(tprs), "prcs": np.array(prcs),
        "mfpr": mfpr, "mrec": mrec,
        "prevalence": float(y.mean()),
        "mean_auc": float(np.mean(aucs)), "std_auc": float(np.std(aucs)),
        "mean_ap":  float(np.mean(aps)),  "std_ap":  float(np.std(aps)),
        "model_refit": model, "model_name": model_name,
    }


def _draw_roc_panel(ax, res, title_extra=""):
    for tpr_fold in res["tprs"]:
        ax.plot(res["mfpr"], tpr_fold, color="#85C1E9", alpha=0.15, lw=1)
    mean_tpr = res["tprs"].mean(axis=0); std_tpr = res["tprs"].std(axis=0)
    ax.plot(res["mfpr"], mean_tpr, color="#C0392B", lw=2,
            label=f"Mean ROC (AUC = {res['mean_auc']:.3f} "
                  f"\u00b1 {res['std_auc']:.3f})")
    ax.fill_between(res["mfpr"], mean_tpr - std_tpr, mean_tpr + std_tpr,
                    color="#E74C3C", alpha=0.12)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="Random (AUC = 0.5)")
    ax.set_xlabel("False Positive Rate", labelpad=6)
    ax.set_ylabel("True Positive Rate", labelpad=6)
    ax.set_title("ROC Curve (5-fold CV)" + title_extra, fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax_in = ax.inset_axes([0.10, 0.62, 0.25, 0.22])
    ax_in.scatter(res["aucs"], [1] * len(res["aucs"]),
                  color="#C0392B", s=35, zorder=3,
                  edgecolor="white", linewidth=0.6)
    ax_in.axvline(res["mean_auc"], color="black", lw=1.0, linestyle="--", zorder=2)
    ax_in.set_title("Per-fold AUC", fontsize=8, pad=3)
    ax_in.tick_params(axis="x", labelsize=7); ax_in.set_yticks([])
    ax_in.set_xlim(res["aucs"].min() - 0.005, res["aucs"].max() + 0.005)
    ax_in.set_ylim(0.5, 1.5)
    for s in ("top", "right", "left"):
        ax_in.spines[s].set_visible(False)


def _draw_pr_panel(ax, res, title_extra=""):
    for prc_fold in res["prcs"]:
        ax.plot(res["mrec"], prc_fold, color="#85C1E9", alpha=0.15, lw=1)
    mean_prec = res["prcs"].mean(axis=0); std_prec = res["prcs"].std(axis=0)
    ax.plot(res["mrec"], mean_prec, color="#C0392B", lw=2,
            label=f"Mean PR (PR-AUC = {res['mean_ap']:.3f} "
                  f"\u00b1 {res['std_ap']:.3f})")
    ax.fill_between(res["mrec"], mean_prec - std_prec, mean_prec + std_prec,
                    color="#E74C3C", alpha=0.12)
    ax.axhline(res["prevalence"], color="black", lw=0.8, linestyle="--",
               alpha=0.5,
               label=f"Random (PR-AUC = {res['prevalence']:.3f})")
    ax.set_xlabel("Recall", labelpad=6); ax.set_ylabel("Precision", labelpad=6)
    ax.set_title("Precision-Recall Curve (5-fold CV)" + title_extra,
                 fontweight="bold", pad=8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def plot_roc_pr(results_by_mode, outdir, comparison, model_name):
    model_disp = "Random Forest" if model_name == "rf" else "XGBoost"
    if comparison in ("A", "B", "all"):
        mode_label = ("Switchers vs Stable" if comparison == "all"
                      else COMPARISON_SPEC[comparison]["label"])
        res = results_by_mode[comparison]
        fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 5))
        _draw_roc_panel(ax_roc, res); _draw_pr_panel(ax_pr, res)
        fig.suptitle(f"{model_disp}: Parental Features Predict "
                     f"Compartment Switching\n({mode_label})",
                     fontweight="bold", fontsize=12, y=1.02)
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1c_{model_name}_ROC_PR_{comparison}")
        return
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for r_idx, comp in enumerate(["A", "B"]):
        if comp not in results_by_mode: continue
        res = results_by_mode[comp]; spec = COMPARISON_SPEC[comp]
        _draw_roc_panel(axes[r_idx, 0], res, title_extra=f"\n({spec['label']})")
        _draw_pr_panel(axes[r_idx, 1], res, title_extra=f"\n({spec['label']})")
    fig.suptitle(f"{model_disp}: Parental Features Predict Compartment Switching"
                 "\nCompartment-resolved (A and B)",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    save_fig(fig, outdir, f"Fig1c_{model_name}_ROC_PR_both")


# =============================================================================
# Fig 1d — SHAP
# =============================================================================

def write_shap_per_row(sv, row_ids, feat_names, outdir, model_name, comp,
                       row_id_col):
    """Per-row SHAP values for one comparison. Columns: row_id_col,
    comparison, model, <feat>_SHAP for each feature."""
    comp_tag = COMPARISON_SPEC[comp]["short"] if comp in COMPARISON_SPEC else "all"
    out = pd.DataFrame(sv, columns=[f"{f}_SHAP" for f in feat_names])
    out.insert(0, row_id_col, row_ids)
    out.insert(1, "comparison", comp_tag)
    out.insert(2, "model", model_name)
    path = os.path.join(outdir, f"SHAP_per_row_{model_name}_{comp}.txt")
    out.to_csv(path, sep="\t", index=False)
    print(f"  Per-row SHAP -> {path}")


def _shap_values(model, X):
    explainer = shap.TreeExplainer(model)
    shap_vals = explainer.shap_values(X)
    if isinstance(shap_vals, list):
        sv = shap_vals[1] if len(shap_vals) > 1 else shap_vals[0]
    elif np.ndim(shap_vals) == 3:
        sv = shap_vals[:, :, 1]
    else:
        sv = shap_vals
    return np.asarray(sv)


def _draw_shap_panels(ax_bar, ax_bee, sv, X, feat_names, title_extra=""):
    mean_shap = np.abs(sv).mean(axis=0).ravel()
    order = np.argsort(mean_shap)[::-1]
    top_idx = int(np.argmax(mean_shap))
    bar_colors = ["#C0392B" if int(i) == top_idx else "#1A5276" for i in order]
    bars = ax_bar.barh([feat_names[i] for i in order[::-1]],
                       [mean_shap[i] for i in order[::-1]],
                       color=bar_colors[::-1], edgecolor="none", height=0.6)
    for bar, val in zip(bars, [mean_shap[i] for i in order[::-1]]):
        ax_bar.text(bar.get_width() + max(mean_shap) * 0.012,
                    bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}", va="center", ha="left", fontsize=8)
    ax_bar.set_xlabel("Mean |SHAP value|\n(magnitude)", labelpad=6)
    ax_bar.set_title("Feature Importance (magnitude)" + title_extra,
                     fontweight="bold", pad=8)
    ax_bar.set_xlim(0, max(mean_shap) * 1.18)

    feat_order = [feat_names[i] for i in order]
    X_df  = pd.DataFrame(X,  columns=feat_names)[feat_order]
    sv_df = pd.DataFrame(sv, columns=feat_names)[feat_order]
    plt.sca(ax_bee)
    shap.summary_plot(sv_df.values, X_df.values, feature_names=feat_order,
                      plot_type="dot", show=False, color_bar=True,
                      sort=False, plot_size=None, cmap=GREEN_YELLOW)
    ax_bee.set_xlabel("SHAP value (impact on switcher probability)", labelpad=6)
    ax_bee.set_title("Feature Importance (direction)" + title_extra,
                     fontweight="bold", pad=8)


def plot_shap(results_by_mode, feat_names, X_by_mode, outdir, comparison,
              model_name, row_id_col):
    model_disp = "Random Forest" if model_name == "rf" else "XGBoost"
    sv_store = {}
    if comparison in ("A", "B", "all"):
        res = results_by_mode[comparison]; X, row_ids = X_by_mode[comparison]
        sv = _shap_values(res["model_refit"], X)
        sv_store[comparison] = (sv, X)
        write_shap_per_row(sv, row_ids, feat_names, outdir, model_name,
                           comparison, row_id_col)
        fig, (ax_bar, ax_bee) = plt.subplots(
            1, 2, figsize=(11, max(4.5, 0.3 * len(feat_names) + 2)),
            gridspec_kw={"width_ratios": [1.0, 1.2]})
        mode_label = ("Switchers vs Stable" if comparison == "all"
                      else COMPARISON_SPEC[comparison]["label"])
        _draw_shap_panels(ax_bar, ax_bee, sv, X, feat_names)
        fig.suptitle(f"Parental Features \u2192 Compartment Switch Prediction "
                     f"(SHAP, {model_disp})\n({mode_label})",
                     fontweight="bold", fontsize=12, y=1.02)
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1d_SHAP_{model_name}_{comparison}")
        return sv_store

    fig, axes = plt.subplots(
        2, 2, figsize=(11, max(9, 0.6 * len(feat_names) + 4)),
        gridspec_kw={"width_ratios": [1.0, 1.2]})
    for r_idx, comp in enumerate(["A", "B"]):
        if comp not in results_by_mode: continue
        res = results_by_mode[comp]; X, row_ids = X_by_mode[comp]
        sv = _shap_values(res["model_refit"], X)
        sv_store[comp] = (sv, X)
        write_shap_per_row(sv, row_ids, feat_names, outdir, model_name, comp,
                           row_id_col)
        spec = COMPARISON_SPEC[comp]
        _draw_shap_panels(axes[r_idx, 0], axes[r_idx, 1], sv, X, feat_names,
                          title_extra=f"\n({spec['label']})")
    fig.suptitle(f"Parental Features \u2192 Compartment Switch Prediction "
                 f"(SHAP, {model_disp})"
                 "\nCompartment-resolved (A on top, B on bottom)",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    save_fig(fig, outdir, f"Fig1d_SHAP_{model_name}_both")
    return sv_store


# =============================================================================
# SHAP dependence plots
# =============================================================================

def _draw_dependence_panel(ax, x_vals, sv_vals, c_vals, feat_x, feat_c,
                            xlim=None, ylim=None):
    sc = ax.scatter(x_vals, sv_vals, c=c_vals, cmap="plasma", s=8, alpha=0.55,
                    edgecolors="none", rasterized=True)
    ax.axhline(0, color="black", lw=0.6, alpha=0.5)
    ax.set_xlabel(f"log1p({feat_x}) (Parental signal)", labelpad=6)
    ax.set_ylabel(f"SHAP value for {feat_x}\n(impact on switcher probability)",
                  labelpad=6)
    ax.set_title(f"SHAP Dependence: {feat_x}\n(colored by {feat_c})",
                 fontweight="bold", pad=8)
    if xlim is not None: ax.set_xlim(*xlim)
    if ylim is not None: ax.set_ylim(*ylim)
    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(f"log1p({feat_c})", fontsize=9)
    cbar.ax.tick_params(labelsize=8)


def _resolve_pair(pair, feat_names):
    name_to_idx_ci = {f.lower(): i for i, f in enumerate(feat_names)}
    a_in, b_in = pair
    a_idx = name_to_idx_ci.get(a_in.lower())
    b_idx = name_to_idx_ci.get(b_in.lower())
    if a_idx is None or b_idx is None:
        print(f"  [WARN] '{a_in}' or '{b_in}' not in feature list "
              f"{feat_names}; skipping.")
        return None
    return (feat_names[a_idx], feat_names[b_idx], a_idx, b_idx)


def plot_shap_dependence(sv_store, feat_names, pairs, outdir, comparison,
                         model_name, zoom_pct=99):
    if not pairs:
        print("[INFO] No --shap_dependence_partners; skipping dependence plots.")
        return

    def _draw_for_mode(ax_l, ax_r, sv, X, fa, fb, ia, ib,
                       xla=None, yla=None, xlb=None, ylb=None):
        _draw_dependence_panel(ax_l, np.log1p(X[:, ia]), sv[:, ia],
                               np.log1p(X[:, ib]), fa, fb, xla, yla)
        _draw_dependence_panel(ax_r, np.log1p(X[:, ib]), sv[:, ib],
                               np.log1p(X[:, ia]), fb, fa, xlb, ylb)

    def _compute_pct(arr, pct):
        lo = np.percentile(arr, 100 - pct); hi = np.percentile(arr, pct)
        pad = 0.04 * (hi - lo) if hi > lo else 0
        return (lo - pad, hi + pad)

    for pair in pairs:
        resolved = _resolve_pair(pair, feat_names)
        if resolved is None: continue
        fa, fb, ia, ib = resolved
        print(f"[INFO] SHAP dependence pair ({comparison}, {model_name}): "
              f"{fa} <-> {fb}")

        if comparison in ("A", "B", "all"):
            sv, X = sv_store[comparison]
            fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4.8))
            _draw_for_mode(ax_l, ax_r, sv, X, fa, fb, ia, ib)
            mode_label = ("Switchers vs Stable" if comparison == "all"
                          else COMPARISON_SPEC[comparison]["label"])
            fig.suptitle(f"({mode_label})", fontweight="bold", fontsize=10, y=1.02)
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{fa}_vs_{fb}_{model_name}_{comparison}")
            # zoomed version
            xla = _compute_pct(np.log1p(X[:, ia]), zoom_pct)
            yla = _compute_pct(sv[:, ia], zoom_pct)
            xlb = _compute_pct(np.log1p(X[:, ib]), zoom_pct)
            ylb = _compute_pct(sv[:, ib], zoom_pct)
            fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4.8))
            _draw_for_mode(ax_l, ax_r, sv, X, fa, fb, ia, ib, xla, yla, xlb, ylb)
            fig.suptitle(f"({mode_label})  —  zoomed to {zoom_pct}th percentile",
                         fontweight="bold", fontsize=10, y=1.02)
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{fa}_vs_{fb}_{model_name}_{comparison}_zoomed")
        else:
            fig, axes = plt.subplots(2, 2, figsize=(12, 9.6))
            for r_idx, comp in enumerate(["A", "B"]):
                if comp not in sv_store: continue
                sv, X = sv_store[comp]
                _draw_for_mode(axes[r_idx, 0], axes[r_idx, 1],
                               sv, X, fa, fb, ia, ib)
                spec = COMPARISON_SPEC[comp]
                axes[r_idx, 0].set_title(axes[r_idx, 0].get_title() +
                                         f"\n[{spec['label']}]",
                                         fontweight="bold", pad=8)
                axes[r_idx, 1].set_title(axes[r_idx, 1].get_title() +
                                         f"\n[{spec['label']}]",
                                         fontweight="bold", pad=8)
            fig.suptitle("SHAP Dependence — A on top, B on bottom",
                         fontweight="bold", fontsize=12, y=1.00)
            plt.tight_layout()
            save_fig(fig, outdir,
                     f"Fig1d_dependence_{fa}_vs_{fb}_{model_name}_both")


def parse_dependence_pairs(raw_pairs):
    if not raw_pairs: return []
    out = []
    for s in raw_pairs:
        if ":" not in s:
            print(f"[WARN] '{s}' missing ':'; skipping."); continue
        parts = s.split(":")
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            print(f"[WARN] '{s}' malformed; skipping."); continue
        out.append((parts[0].strip(), parts[1].strip()))
    return out


# =============================================================================
# Fig 1e — PCA
# =============================================================================

def _draw_pca_panel(ax, df_subset, feat_names, cats_fg, cats_bg):
    if df_subset.empty:
        ax.set_visible(False); return None
    Xsc = StandardScaler().fit_transform(np.log1p(df_subset[feat_names].values))
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xsc)
    var = pca.explained_variance_ratio_ * 100
    for cat in cats_bg:
        mask = df_subset["switch_cat"].values == cat
        if mask.sum() == 0: continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=SWITCH_COLORS[cat], s=4, alpha=0.06,
                   linewidths=0, rasterized=True, label=SWITCH_LABELS[cat])
    for cat in cats_fg:
        mask = df_subset["switch_cat"].values == cat
        if mask.sum() == 0: continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=SWITCH_COLORS[cat], s=18, alpha=0.55,
                   linewidths=0, rasterized=True, label=SWITCH_LABELS[cat])
    xlo, xhi = np.percentile(coords[:, 0], [1, 99])
    ylo, yhi = np.percentile(coords[:, 1], [1, 99])
    pad_x = 0.05 * (xhi - xlo); pad_y = 0.05 * (yhi - ylo)
    ax.set_xlim(xlo - pad_x, xhi + pad_x); ax.set_ylim(ylo - pad_y, yhi + pad_y)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)", labelpad=6)
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)", labelpad=6)
    legend_handles = [plt.Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor=SWITCH_COLORS[c],
                                  markersize=8, label=SWITCH_LABELS[c])
                      for c in (cats_fg + cats_bg)
                      if c in df_subset["switch_cat"].values]
    ax.legend(handles=legend_handles, frameon=False, fontsize=8.5, loc="best")
    return pca


def plot_pca(df, feat_names, outdir, comparison):
    if comparison == "all":
        df_sub = df[df["switch_cat"].isin(
            ["Stable_A", "Stable_B", "A_to_B", "B_to_A"])].dropna(subset=feat_names)
        fig, ax = plt.subplots(figsize=(6.2, 5.6))
        _draw_pca_panel(ax, df_sub, feat_names,
                        cats_fg=["A_to_B", "B_to_A"],
                        cats_bg=["Stable_A", "Stable_B"])
        ax.set_title("PCA: Parental Feature Space (log1p)"
                     "\nby Compartment Switch Outcome",
                     fontweight="bold", pad=10)
        plt.tight_layout()
        save_fig(fig, outdir, "Fig1e_PCA_marks_all")
        return
    if comparison in ("A", "B"):
        spec = COMPARISON_SPEC[comparison]
        df_sub = (df[df["switch_cat"].isin([spec["neg"], spec["pos"]])]
                  .dropna(subset=feat_names))
        fig, ax = plt.subplots(figsize=(6.2, 5.6))
        _draw_pca_panel(ax, df_sub, feat_names,
                        cats_fg=[spec["pos"]], cats_bg=[spec["neg"]])
        ax.set_title(f"PCA: Parental Feature Space (log1p)"
                     f"\n{spec['label']}",
                     fontweight="bold", pad=10)
        plt.tight_layout()
        save_fig(fig, outdir, f"Fig1e_PCA_marks_{comparison}")
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    for ax, comp in zip(axes, ["A", "B"]):
        spec = COMPARISON_SPEC[comp]
        df_sub = (df[df["switch_cat"].isin([spec["neg"], spec["pos"]])]
                  .dropna(subset=feat_names))
        _draw_pca_panel(ax, df_sub, feat_names,
                        cats_fg=[spec["pos"]], cats_bg=[spec["neg"]])
        ax.set_title(f"PCA: {spec['label']}", fontweight="bold", pad=10)
    fig.suptitle("PCA: Parental Feature Space (log1p)"
                 "\nCompartment-resolved (A on left, B on right)",
                 fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    save_fig(fig, outdir, "Fig1e_PCA_marks_both")


# =============================================================================
# Statistics
# =============================================================================

def save_statistics(df, feat_names, outdir, comparison):
    def _per_feat(neg_cats, pos_cats, label):
        out = []
        neg_mask = df["switch_cat"].isin(neg_cats)
        pos_mask = df["switch_cat"].isin(pos_cats)
        for feat in feat_names:
            neg_vals = df[neg_mask][feat].dropna().values
            pos_vals = df[pos_mask][feat].dropna().values
            if len(neg_vals) == 0 or len(pos_vals) == 0:
                continue
            stat, p = mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")
            d, mag = cliffs_delta_fast(pos_vals, neg_vals)
            out.append({"comparison": label, "feature": feat,
                        "n_positive": len(pos_vals),
                        "n_negative": len(neg_vals),
                        "median_positive": np.median(pos_vals),
                        "median_negative": np.median(neg_vals),
                        "mannwhitney_U": stat, "p_value": p,
                        "cliffs_delta": d, "effect_magnitude": mag})
        return out

    rows = []
    if comparison == "A":
        rows += _per_feat(["Stable_A"], ["A_to_B"], "Stable_A_vs_AtoB")
    elif comparison == "B":
        rows += _per_feat(["Stable_B"], ["B_to_A"], "Stable_B_vs_BtoA")
    elif comparison == "both":
        rows += _per_feat(["Stable_A"], ["A_to_B"], "Stable_A_vs_AtoB")
        rows += _per_feat(["Stable_B"], ["B_to_A"], "Stable_B_vs_BtoA")
    elif comparison == "all":
        rows += _per_feat(["Stable_A", "Stable_B"], ["B_to_A", "A_to_B"],
                          "Switchers_vs_Stable")
    if not rows:
        print("[WARN] No statistics rows generated."); return
    res_df = pd.DataFrame(rows)
    res_df["padj_BH"] = multipletests(res_df["p_value"], method="fdr_bh")[1]
    out = os.path.join(outdir, f"Statistics_{comparison}.txt")
    res_df.to_csv(out, sep="\t", index=False)
    print(f"[INFO] Statistics saved -> {out}")


# =============================================================================
# Feature matrix output (NEW — the exact prediction matrix)
# =============================================================================

def save_feature_matrix(sub, feat_names, row_id_col, outdir, comp_tag):
    """Write the exact rows/features/label that entered the model — auditable.
    Columns: row_id_col, <features...>, switch_cat, is_switcher (if present)."""
    cols = [row_id_col] + feat_names + ["switch_cat"]
    if "is_switcher" in sub.columns:
        cols.append("is_switcher")
    cols = [c for c in cols if c in sub.columns]
    out = sub[cols].copy()
    path = os.path.join(outdir, f"Feature_matrix_{comp_tag}.txt")
    out.to_csv(path, sep="\t", index=False)
    print(f"[INFO] Feature matrix -> {path} ({len(out):,} x {len(cols)})")


# =============================================================================
# Orchestrator
# =============================================================================

def run_ml_for_comparison(df, feat_names, comparison, model_name, outdir,
                           row_id_col):
    results, Xs = {}, {}
    if comparison == "all":
        sub = df[df["switch_cat"].isin(
            ["B_to_A", "A_to_B", "Stable_A", "Stable_B"])].copy()
        sub["is_switcher"] = sub["switch_cat"].isin(
            ["B_to_A", "A_to_B"]).astype(int)
        X = sub[feat_names].values; y = sub["is_switcher"].values
        if len(X) < 50 or y.sum() < 5:
            print("[WARN] Insufficient data for model (all). Skipping.")
            return results, Xs
        print(f"[INFO] {model_name} (all): N={len(X)}, switchers={y.sum()}, "
              f"stable={(y == 0).sum()}, prevalence={y.mean():.3f}")
        save_feature_matrix(sub, feat_names, row_id_col, outdir, "all")
        results["all"] = _run_model_cv(X, y, model_name)
        Xs["all"] = (X, sub[row_id_col].values)
        return results, Xs

    comps = ["A", "B"] if comparison == "both" else [comparison]
    for comp in comps:
        sub, spec = subset_for_comparison(df, comp)
        if sub["is_switcher"].sum() < 5 or len(sub) < 50:
            print(f"[WARN] Insufficient data for model ({comp}). Skipping.")
            continue
        X = sub[feat_names].values; y = sub["is_switcher"].values
        print(f"[INFO] {model_name} ({spec['label']}): N={len(X)}, "
              f"positive={y.sum()} ({spec['pos']}), "
              f"negative={(y == 0).sum()} ({spec['neg']}), "
              f"prevalence={y.mean():.3f}")
        save_feature_matrix(sub, feat_names, row_id_col, outdir, comp)
        results[comp] = _run_model_cv(X, y, model_name)
        Xs[comp] = (X, sub[row_id_col].values)
    return results, Xs


def main():
    args = parse_args()
    # Namespace outputs by feature_mode so 'baseline' and 'dynamics' runs land
    # in separate subdirectories and never overwrite each other. This matters
    # because the publication-defensible practice is to RUN BOTH and compare,
    # not pick one. Preserving both output sets makes that easy.
    args.outdir = os.path.join(args.outdir, args.features)
    os.makedirs(args.outdir, exist_ok=True)
    print(f"[INFO] Output directory (namespaced by feature mode): {args.outdir}")

    if args.model == "xgboost" and not _HAS_XGB:
        sys.exit("[ERROR] --model xgboost requested but xgboost not installed.")
    if args.comparison == "all" and sp is None:
        print("[WARN] 'all' mode requested but scikit-posthocs not installed. "
              "Dunn's post-hoc and brackets will be skipped.")
    print(f"[INFO] Model backend: {args.model}")

    stats = list(DEFAULT_STATS) if args.stats is None else args.stats
    df, feat_names, feat_map, row_id_col = load_data(
        args.cutrun_file, args.merged_file,
        marks=args.marks, stats=stats, include_ctcf=args.include_ctcf,
        cutrun_key_override=args.cutrun_key, hic_key_override=args.hic_key,
        row_id_col=args.row_id, min_overlap_frac=args.min_overlap,
        feature_mode=args.features,
    )

    out_genes = os.path.join(args.outdir, "Switch_classified.txt")
    df.to_csv(out_genes, sep="\t", index=False)
    print(f"[INFO] Classified bin table -> {out_genes}")

    print(f"\n[INFO] Comparison mode: {args.comparison}")
    print("\n[INFO] Fig 1a: switch distribution...")
    plot_switch_distribution(df, args.outdir, args.comparison)

    print("\n[INFO] Fig 1b: violins...")
    if args.comparison == "all":
        plot_violin_marks_all(df, feat_names, args.outdir)
    else:
        plot_violin_marks_AB(df, feat_names, args.outdir, args.comparison)

    print(f"\n[INFO] ML pipeline ({args.model} + CV)...")
    results_by_mode, X_by_mode = run_ml_for_comparison(
        df, feat_names, args.comparison, args.model, args.outdir, row_id_col)
    if not results_by_mode:
        print("[ERROR] No model results produced; aborting downstream plots.")
        sys.exit(1)

    print("\n[INFO] Fig 1c: ROC + PR...")
    plot_roc_pr(results_by_mode, args.outdir, args.comparison, args.model)

    print("\n[INFO] Fig 1d: SHAP bar + beeswarm...")
    sv_store = plot_shap(results_by_mode, feat_names, X_by_mode,
                         args.outdir, args.comparison, args.model, row_id_col)

    pairs = parse_dependence_pairs(args.shap_dependence_partners)
    if pairs:
        print(f"\n[INFO] Fig 1d: SHAP dependence pairs ({pairs})...")
        plot_shap_dependence(sv_store, feat_names, pairs,
                             args.outdir, args.comparison, args.model,
                             zoom_pct=99)

    print("\n[INFO] Fig 1e: PCA...")
    plot_pca(df, feat_names, args.outdir, args.comparison)

    print("\n[INFO] Statistics...")
    save_statistics(df, feat_names, args.outdir, args.comparison)

    print(f"\n[INFO] Q1 complete ({args.model}) "
          f"-> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()