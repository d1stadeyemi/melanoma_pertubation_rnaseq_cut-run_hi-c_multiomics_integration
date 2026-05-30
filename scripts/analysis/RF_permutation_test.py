"""
Permutation Test — empirical null distribution for RF / XGBoost ROC-AUC & PR-AUC.

PURPOSE
-------
Establishes whether the observed cross-validated AUC from RF_SHAP.py is higher
than what the same model achieves on the SAME features with the labels randomly
shuffled. Shuffling destroys any real feature->label relationship, so the
distribution of AUCs over many shuffles is the empirical null. The observed AUC
is then placed against that null to get a permutation p-value
(Phipson & Smyth 2010 (+1)-corrected, so p is never exactly zero).

WHY THIS SCRIPT CONSUMES RF_SHAP.py OUTPUT DIRECTLY
---------------------------------------------------
It reads the SAME `Switch_classified.txt` that RF_SHAP.py writes, and detects
the feature columns the SAME WAY RF_SHAP.py names them. Critically:

  1. FEATURES ARE AUTO-DETECTED, not hardcoded. RF_SHAP renames feature columns
     to display names that INCLUDE the stat suffix, e.g. 'H3K27me3_mean',
     'H3K9me3_std' (baseline) or 'd_H3K27me3_mean' (dynamics). An older
     hardcoded list of bare mark names ['H3K27me3', ...] matched NOTHING and
     the script exited. We now detect any column matching
     '<mark>_<stat>' or 'd_<mark>_<stat>' so the permutation test works on
     whatever feature set RF_SHAP actually produced.

  2. THE MODEL IS IDENTICAL to RF_SHAP.py. A permutation null is only a valid
     reference for the observed statistic if BOTH are produced by the exact
     same estimator. RF_SHAP uses RF(n_estimators=500, max_depth=12,
     class_weight='balanced_subsample'); XGBoost params likewise. We import
     those exact settings here. (The previous version used max_depth=6 /
     class_weight='balanced', producing a null from a DIFFERENT model than the
     one whose AUC we were testing — an invalid comparison.)

  3. FEATURE MODE IS INFERRED FROM THE FILE. Because RF_SHAP namespaces its
     output (<outdir>/baseline/, <outdir>/dynamics/), simply pointing this
     script at the Switch_classified.txt inside one of those folders tests that
     mode. The detected feature names are printed so the null is unambiguous.

CONTRACT WITH RF_SHAP.py
------------------------
  Input  : Switch_classified.txt   (the per-bin table RF_SHAP writes; contains
                                     the renamed feature columns + 'switch_cat')
  Match  : --comparison must match the RF_SHAP run you are validating
           (A / B / both / all).

USAGE
-----
    # validate the baseline run:
    python scripts/RF_permutation_test.py \
        results/exploratory/RF_SHAP/baseline/Switch_classified.txt \
        --outdir results/exploratory/RF_permutation_test/baseline \
        --comparison both --model rf --n_perm 1000 --n_jobs 8

    # validate the dynamics run (features auto-detected as d_<mark>_<stat>):
    python scripts/RF_permutation_test.py \
        results/exploratory/RF_SHAP/dynamics/Switch_classified.txt \
        --outdir results/exploratory/RF_permutation_test/dynamics \
        --comparison both --model rf --n_perm 1000

OUTPUTS (per mode)
------------------
    Fig1c_permutation_null_<mode>.pdf/png
    Permutation_results_<mode>.txt     (per-permutation null ROC & PR)
    Permutation_summary_<mode>.txt     (observed, null mean/sd/min/max, p)

NOTE ON n_perm
--------------
The smallest empirical p achievable is 1/(n_perm+1). With n_perm=100 the floor
is ~0.0099, which is fine for "clearly better than chance" but cannot resolve
small p-values. For a publication claim use n_perm>=1000 (floor ~0.001). The
default below is 1000; lower it only for quick checks.
"""

import argparse
import os
import re
import sys
import time
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score

# XGBoost is optional — only needed if --model xgboost is requested.
try:
    from xgboost import XGBClassifier
    _HAS_XGB = True
except ImportError:
    XGBClassifier = None
    _HAS_XGB = False

plt.rcParams.update({
    "font.family":       "Liberation Sans",
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

# -----------------------------------------------------------------------------
# MODEL CONFIG — MUST be identical to RF_SHAP.py so the null matches the model
# whose AUC we are testing. These are copied verbatim from RF_SHAP.py.
# If RF_SHAP.py's RF_PARAMS / XGB_PARAMS ever change, change them here too.
# -----------------------------------------------------------------------------
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

# Comparison spec — identical to RF_SHAP.py so subsets are defined the same way.
COMPARISON_SPEC = {
    "A": {"neg": "Stable_A", "pos": "A_to_B", "label": "Stable_A vs A\u2192B"},
    "B": {"neg": "Stable_B", "pos": "B_to_A", "label": "Stable_B vs B\u2192A"},
}

# Marks/stats vocabulary used to AUTO-DETECT feature columns by pattern.
# Matches RF_SHAP.py's DEFAULT_MARKS x DEFAULT_STATS, plus the 'd_' dynamics
# prefix and optional CTCF. We detect columns of the form:
#     <mark>_<stat>        (baseline, e.g. H3K27me3_mean)
#     d_<mark>_<stat>      (dynamics, e.g. d_H3K27me3_mean)
KNOWN_MARKS = ("H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3", "CTCF")
KNOWN_STATS = ("mean", "max", "std", "coverage")


def parse_args():
    p = argparse.ArgumentParser(
        description="RF/XGBoost permutation test (consumes RF_SHAP.py output).")
    p.add_argument("classified_file",
                   help="Switch_classified.txt from RF_SHAP.py "
                        "(use the one inside the baseline/ or dynamics/ "
                        "subdir to test that feature mode).")
    p.add_argument("--outdir",
                   default="results/exploratory/RF_permutation_test",
                   help="Output directory.")
    p.add_argument("--comparison",
                   choices=["A", "B", "both", "all"], default="both",
                   help="Match this to the RF_SHAP.py run being validated. "
                        "Default: both.")
    p.add_argument("--model", choices=["rf", "xgboost"], default="rf",
                   help="Must match the RF_SHAP.py --model used. Default: rf.")
    p.add_argument("--n_perm", type=int, default=1000,
                   help="Number of label permutations (default: 1000). "
                        "Smallest achievable p is 1/(n_perm+1).")
    p.add_argument("--n_jobs", type=int, default=-1,
                   help="Cores for the model (-1 = all; default -1).")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed (default 42).")
    # Escape hatch: let the user name feature columns explicitly if their file
    # uses non-standard names that auto-detection misses.
    p.add_argument("--features", nargs="+", default=None,
                   help="Explicit feature column names to use, overriding "
                        "auto-detection. Use the exact column names from "
                        "Switch_classified.txt (e.g. H3K27me3_mean "
                        "H3K9me3_std).")
    return p.parse_args()


def detect_feature_columns(df, explicit=None):
    """
    Auto-detect the feature columns RF_SHAP.py wrote, by matching the naming
    convention it uses: '<mark>_<stat>' or 'd_<mark>_<stat>' (case-insensitive
    on mark and stat). Returns the list of matching column names, in file order.

    This is the key compatibility fix: RF_SHAP renames features to include the
    stat suffix, so a bare-mark-name lookup finds nothing. We instead accept any
    column whose name parses as a known mark + known stat, optionally with the
    'd_' dynamics prefix.

    If `explicit` is provided, we use exactly those columns (after verifying
    they exist), bypassing detection entirely.
    """
    if explicit:
        missing = [c for c in explicit if c not in df.columns]
        if missing:
            sys.exit(f"[ERROR] --features columns not in file: {missing}\n"
                     f"  Available columns: {list(df.columns)}")
        print(f"[INFO] Using {len(explicit)} explicitly-specified feature "
              f"columns.")
        return list(explicit)

    marks_lc = {m.lower() for m in KNOWN_MARKS}
    stats_lc = {s.lower() for s in KNOWN_STATS}
    # Pattern: optional 'd_' prefix, then mark, then '_', then stat.
    # We parse rather than regex-per-combo so the message can explain misses.
    detected = []
    detected_baseline, detected_dynamics = [], []
    for col in df.columns:
        lc = col.lower()
        is_dyn = lc.startswith("d_")
        core = lc[2:] if is_dyn else lc
        # split into mark_part + stat_part on the LAST underscore
        if "_" not in core:
            continue
        mark_part, _, stat_part = core.rpartition("_")
        if mark_part in marks_lc and stat_part in stats_lc:
            detected.append(col)
            (detected_dynamics if is_dyn else detected_baseline).append(col)

    if detected_baseline and detected_dynamics:
        print(f"[INFO] Detected BOTH baseline ({len(detected_baseline)}) and "
              f"dynamics ({len(detected_dynamics)}) feature columns "
              f"(file came from --features all).")
    elif detected_dynamics:
        print(f"[INFO] Detected {len(detected_dynamics)} DYNAMICS feature "
              f"columns (d_<mark>_<stat>).")
    elif detected_baseline:
        print(f"[INFO] Detected {len(detected_baseline)} BASELINE feature "
              f"columns (<mark>_<stat>).")
    return detected


def cv_aucs(X, y, model_name, n_jobs, seed):
    """
    One 5-fold stratified-CV pass -> (mean ROC-AUC, mean PR-AUC).

    Uses the SAME model and the SAME CV scheme as RF_SHAP.py so the null is a
    valid reference for the observed statistic. For XGBoost we set
    scale_pos_weight per-fold from the TRAIN fold only (matching RF_SHAP.py's
    handling and avoiding leakage of test-fold prevalence into the model).
    """
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    roc_per_fold, pr_per_fold = [], []
    for tr_idx, te_idx in cv.split(X, y):
        y_tr = y[tr_idx]
        if model_name == "xgboost":
            n_pos = max(int(y_tr.sum()), 1)
            n_neg = int((y_tr == 0).sum())
            model = XGBClassifier(random_state=seed,
                                  scale_pos_weight=n_neg / n_pos,
                                  **XGB_PARAMS)
        else:
            model = RandomForestClassifier(random_state=seed, **RF_PARAMS)
        model.fit(X[tr_idx], y_tr)
        y_prob = model.predict_proba(X[te_idx])[:, 1]
        roc_per_fold.append(roc_auc_score(y[te_idx], y_prob))
        pr_per_fold.append(average_precision_score(y[te_idx], y_prob))
    return float(np.mean(roc_per_fold)), float(np.mean(pr_per_fold))


def empirical_p(observed, null_arr):
    """One-sided empirical p with Phipson & Smyth (2010) +1 correction.
    p = (#{null >= observed} + 1) / (n_perm + 1). Never returns exactly 0,
    which is the statistically correct behaviour for a finite permutation set."""
    n = len(null_arr)
    return (np.sum(null_arr >= observed) + 1) / (n + 1)


def build_subset(df, comparison, feat_names):
    """
    Build (X, y, label) tuples per comparison sub-mode, identically to
    RF_SHAP.py's subsetting. Returns a dict keyed by 'A'/'B'/'all'.

    'all'  -> switchers (B_to_A + A_to_B) vs stable (Stable_A + Stable_B)
    'A'    -> Stable_A vs A_to_B
    'B'    -> Stable_B vs B_to_A
    'both' -> entries 'A' and 'B'
    Rows with NaN in any feature are dropped (same as RF_SHAP).
    """
    out = {}
    if comparison == "all":
        sub = df[df["switch_cat"].isin(
            ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
        )].dropna(subset=feat_names).copy()
        sub["is_switcher"] = sub["switch_cat"].isin(
            ["B_to_A", "A_to_B"]).astype(int)
        out["all"] = (sub[feat_names].values,
                      sub["is_switcher"].values,
                      "Switchers vs Stable")
        return out

    comps = ["A", "B"] if comparison == "both" else [comparison]
    for comp in comps:
        spec = COMPARISON_SPEC[comp]
        sub = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])] \
                .dropna(subset=feat_names).copy()
        sub["is_switcher"] = (sub["switch_cat"] == spec["pos"]).astype(int)
        out[comp] = (sub[feat_names].values,
                     sub["is_switcher"].values,
                     spec["label"])
    return out


def run_permutation(X, y, model_name, n_perm, n_jobs, seed, label):
    """Observed CV AUCs + permutation null. Returns a dict of arrays + stats."""
    print(f"\n[INFO] === Comparison: {label} ===")
    print(f"[INFO] N={len(X)}, positive={int(y.sum())}, "
          f"negative={int((y == 0).sum())}, prevalence={y.mean():.3f}")

    print("[INFO] Computing OBSERVED 5-fold CV AUCs...")
    t0 = time.time()
    obs_roc, obs_pr = cv_aucs(X, y, model_name, n_jobs, seed=seed)
    print(f"[INFO] Observed ROC-AUC = {obs_roc:.4f}")
    print(f"[INFO] Observed PR-AUC  = {obs_pr:.4f}")
    print(f"[INFO] Observed run took {(time.time() - t0):.1f} s")

    print(f"[INFO] Running {n_perm} label permutations "
          f"(empirical p floor = {1.0/(n_perm+1):.4g})...")
    rng = np.random.default_rng(seed)
    null_roc = np.zeros(n_perm)
    null_pr  = np.zeros(n_perm)
    t_start = time.time()
    for i in range(n_perm):
        # Shuffle the labels only; features are untouched. Each permutation gets
        # its own CV seed so fold splits vary, exactly as an independent draw.
        y_shuffled = rng.permutation(y)
        roc_i, pr_i = cv_aucs(X, y_shuffled, model_name, n_jobs,
                              seed=seed + i + 1)
        null_roc[i] = roc_i
        null_pr[i]  = pr_i
        if (i + 1) % 10 == 0 or i == 0:
            elapsed  = time.time() - t_start
            per_perm = elapsed / (i + 1)
            eta = per_perm * (n_perm - (i + 1))
            print(f"  perm {i + 1:>4d}/{n_perm}  "
                  f"ROC={roc_i:.3f} PR={pr_i:.3f}  "
                  f"({elapsed/60:.1f} min elapsed; ETA {eta/60:.1f} min)")
    total_min = (time.time() - t_start) / 60
    print(f"[INFO] Permutations complete in {total_min:.1f} min")

    p_roc = empirical_p(obs_roc, null_roc)
    p_pr  = empirical_p(obs_pr,  null_pr)
    print(f"[INFO] Null ROC: mean={null_roc.mean():.4f} SD={null_roc.std():.4f}")
    print(f"[INFO] Null PR : mean={null_pr.mean():.4f}  SD={null_pr.std():.4f}")
    print(f"[INFO] empirical p (ROC) = {p_roc:.4g}")
    print(f"[INFO] empirical p (PR)  = {p_pr:.4g}")

    return {
        "label": label,
        "obs_roc": obs_roc, "obs_pr": obs_pr,
        "null_roc": null_roc, "null_pr": null_pr,
        "p_roc": p_roc, "p_pr": p_pr,
        "prevalence": float(y.mean()),
    }


def save_tables(res, mode_key, outdir, n_perm, model_name, feat_names):
    """Per-mode permutation results + summary tables (tab-separated)."""
    perm_df = pd.DataFrame({
        "permutation":  np.arange(1, n_perm + 1),
        "null_ROC_AUC": res["null_roc"],
        "null_PR_AUC":  res["null_pr"],
    })
    perm_path = os.path.join(outdir, f"Permutation_results_{mode_key}.txt")
    perm_df.to_csv(perm_path, sep="\t", index=False)
    print(f"  Saved -> {perm_path}")

    # Summary carries the model and feature set so the result is self-documenting.
    feat_str = ";".join(feat_names)
    summary = pd.DataFrame([
        {"metric": "ROC-AUC", "comparison": res["label"], "model": model_name,
         "n_features": len(feat_names), "features": feat_str,
         "observed": res["obs_roc"],
         "null_mean": res["null_roc"].mean(),
         "null_sd":   res["null_roc"].std(),
         "null_min":  res["null_roc"].min(),
         "null_max":  res["null_roc"].max(),
         "n_perm":    n_perm,
         "empirical_p_one_sided": res["p_roc"]},
        {"metric": "PR-AUC",  "comparison": res["label"], "model": model_name,
         "n_features": len(feat_names), "features": feat_str,
         "observed": res["obs_pr"],
         "null_mean": res["null_pr"].mean(),
         "null_sd":   res["null_pr"].std(),
         "null_min":  res["null_pr"].min(),
         "null_max":  res["null_pr"].max(),
         "n_perm":    n_perm,
         "empirical_p_one_sided": res["p_pr"]},
    ])
    summary_path = os.path.join(outdir, f"Permutation_summary_{mode_key}.txt")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"  Saved -> {summary_path}")


def _draw_null_panel(ax, null_arr, observed, prev_for_random, p_emp,
                     n_perm, metric_label, random_label):
    """One null-distribution histogram panel with observed/null-mean/random lines."""
    ax.hist(null_arr, bins=30, color="#85C1E9",
            edgecolor="white", linewidth=0.5, alpha=0.95)
    ax.axvline(observed, color="#C0392B", lw=2.0,
               label=f"Observed = {observed:.3f}")
    ax.axvline(null_arr.mean(), color="black", lw=1.0, linestyle="--",
               label=f"Null mean = {null_arr.mean():.3f}")
    ax.axvline(prev_for_random, color="grey", lw=0.8, linestyle=":",
               label=random_label)
    ax.set_xlabel(f"{metric_label} (5-fold CV)", labelpad=6)
    ax.set_ylabel("Permutations", labelpad=6)
    
    # Format p: at the permutation floor (observed beat all shuffles) the value
    # is censored, so report it as an inequality (Phipson & Smyth 2010; North
    # et al. 2002) rather than a spuriously precise decimal. Above the floor,
    # print the computed value.
    p_floor = 1.0 / (n_perm + 1)
    if p_emp <= p_floor + 1e-12:           # at the resolution limit
        # round the threshold to one significant figure for a clean label,
        # e.g. n_perm=1000 -> "p < 0.001"; n_perm=100 -> "p < 0.01"
        thresh = 10 ** np.floor(np.log10(p_floor) + 1)
        p_str = f"p < {thresh:g}"
    else:
        p_str = f"empirical p = {p_emp:.3g}"
    ax.set_title(
        f"{metric_label} Null Distribution\n"
        f"{p_str} (n = {n_perm})",
        fontweight="bold", pad=8,
    )

    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def plot_nulls(results_by_mode, outdir, comparison, n_perm):
    """A/B/all -> 1x2 (ROC + PR). both -> 2x2 (A top, B bottom)."""
    if comparison in ("A", "B", "all"):
        res = results_by_mode[comparison]
        fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 4.5))
        _draw_null_panel(ax_roc, res["null_roc"], res["obs_roc"], 0.5,
                         res["p_roc"], n_perm, "ROC-AUC", "Random (0.5)")
        _draw_null_panel(ax_pr, res["null_pr"],  res["obs_pr"],
                         res["prevalence"], res["p_pr"], n_perm,
                         "PR-AUC", f"Random ({res['prevalence']:.3f})")
        fig.suptitle(f"Permutation Test: {res['label']}",
                     fontweight="bold", fontsize=12, y=1.02)
        plt.tight_layout()
        for ext in ("pdf", "png"):
            out = os.path.join(outdir,
                               f"Fig1c_permutation_null_{comparison}.{ext}")
            fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
            print(f"  Saved -> {out}")
        plt.close(fig)
        return

    # both -> 2x2
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for r_idx, comp in enumerate(["A", "B"]):
        if comp not in results_by_mode:
            continue
        res = results_by_mode[comp]
        _draw_null_panel(axes[r_idx, 0], res["null_roc"], res["obs_roc"], 0.5,
                         res["p_roc"], n_perm, "ROC-AUC", "Random (0.5)")
        _draw_null_panel(axes[r_idx, 1], res["null_pr"],  res["obs_pr"],
                         res["prevalence"], res["p_pr"], n_perm,
                         "PR-AUC", f"Random ({res['prevalence']:.3f})")
        axes[r_idx, 0].set_title(
            axes[r_idx, 0].get_title() + f"\n[{res['label']}]",
            fontweight="bold", pad=8)
        axes[r_idx, 1].set_title(
            axes[r_idx, 1].get_title() + f"\n[{res['label']}]",
            fontweight="bold", pad=8)
    fig.suptitle(
        "Permutation Test: Parental Marks vs Compartment Switching"
        "\nCompartment-resolved (A on top, B on bottom)",
        fontweight="bold", fontsize=12, y=1.00)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Fig1c_permutation_null_both.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.model == "xgboost" and not _HAS_XGB:
        sys.exit("[ERROR] --model xgboost requested but xgboost not installed.")

    print(f"[INFO] Loading {args.classified_file} ...")
    df = pd.read_csv(args.classified_file, sep="\t", low_memory=False)

    if "switch_cat" not in df.columns:
        sys.exit("[ERROR] 'switch_cat' column not found. This file does not "
                 "look like RF_SHAP.py's Switch_classified.txt.")

    # Auto-detect (or accept explicit) feature columns the way RF_SHAP names them.
    feat_names = detect_feature_columns(df, explicit=args.features)
    if len(feat_names) < 2:
        sys.exit(
            "[ERROR] Fewer than 2 feature columns detected.\n"
            "  RF_SHAP.py writes features as '<mark>_<stat>' (baseline) or "
            "'d_<mark>_<stat>' (dynamics), e.g. 'H3K27me3_mean'.\n"
            f"  Columns present: {list(df.columns)}\n"
            "  If your columns use non-standard names, pass them with "
            "--features explicitly."
        )
    print(f"[INFO] Features ({len(feat_names)}): {feat_names}")
    print(f"[INFO] Model: {args.model}   Comparison: {args.comparison}   "
          f"n_perm: {args.n_perm}")

    # Coerce feature columns to numeric (defensive — they should already be).
    df[feat_names] = df[feat_names].apply(pd.to_numeric, errors="coerce")

    subsets = build_subset(df, args.comparison, feat_names)
    if not subsets:
        sys.exit("[ERROR] No subsets built; check --comparison.")

    results_by_mode = {}
    for key, (X, y, label) in subsets.items():
        # Same guards as RF_SHAP.py: need enough rows and enough positives.
        if len(X) < 50 or y.sum() < 5:
            print(f"[WARN] Insufficient data for permutation ({key}). Skipping.")
            continue
        res = run_permutation(X, y, args.model, args.n_perm, args.n_jobs,
                              args.seed, label)
        results_by_mode[key] = res
        save_tables(res, key, args.outdir, args.n_perm, args.model, feat_names)

    if not results_by_mode:
        sys.exit("[ERROR] No permutation results produced.")

    print("\n[INFO] Plotting null distributions...")
    plot_nulls(results_by_mode, args.outdir, args.comparison, args.n_perm)

    print(f"\n[INFO] Permutation test complete -> "
          f"{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()