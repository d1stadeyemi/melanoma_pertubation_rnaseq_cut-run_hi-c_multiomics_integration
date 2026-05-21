"""
Permutation Test — empirical null distribution for RF ROC-AUC and PR-AUC.

Mirrors the --comparison flag from RF_SHAP.py:
    A     : permutation on Stable_A vs A->B subset
    B     : permutation on Stable_B vs B->A subset
    both  : run A and B separately; 4-panel figure (DEFAULT)
    all   : original switcher-vs-stable framing on full dataset

Usage:
    python scripts/RF_permutation_test.py \
        results/exploratory/RF_SHAP/Switch_classified.txt \
        --outdir results/exploratory/RF_permutation_test \
        --comparison both \
        --n_perm 100 \
        --n_jobs 8

Outputs (per mode):
    Fig1c_permutation_null_<mode>.pdf/png
    Permutation_results_<mode>.txt
    Permutation_summary_<mode>.txt
"""

import argparse
import os
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

# Same as RF_SHAP.py
COMPARISON_SPEC = {
    "A": {"neg": "Stable_A", "pos": "A_to_B", "label": "Stable_A vs A\u2192B"},
    "B": {"neg": "Stable_B", "pos": "B_to_A", "label": "Stable_B vs B\u2192A"},
}


def parse_args():
    p = argparse.ArgumentParser(description="RF permutation test")
    p.add_argument("classified_file",
                   help="Path to Switch_classified.txt from RF_SHAP.py")
    p.add_argument("--outdir",
                   default="results/exploratory/RF_permutation_test",
                   help="Output directory")
    p.add_argument("--comparison",
                   choices=["A", "B", "both", "all"], default="both",
                   help="Match this to the RF_SHAP.py mode used. Default: both")
    p.add_argument("--n_perm", type=int, default=100,
                   help="Number of label permutations (default: 100)")
    p.add_argument("--n_jobs", type=int, default=-1,
                   help="Cores for RF (-1 = all; default -1)")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for reproducibility (default: 42)")
    return p.parse_args()


def cv_aucs(X, y, n_jobs, seed):
    """5-fold CV RF -> (mean ROC-AUC, mean PR-AUC)."""
    rf = RandomForestClassifier(
        n_estimators=500, max_depth=6,
        class_weight="balanced", random_state=seed, n_jobs=n_jobs,
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    roc_per_fold, pr_per_fold = [], []
    for tr_idx, te_idx in cv.split(X, y):
        rf.fit(X[tr_idx], y[tr_idx])
        y_prob = rf.predict_proba(X[te_idx])[:, 1]
        roc_per_fold.append(roc_auc_score(y[te_idx], y_prob))
        pr_per_fold.append(average_precision_score(y[te_idx], y_prob))
    return float(np.mean(roc_per_fold)), float(np.mean(pr_per_fold))


def empirical_p(observed, null_arr):
    """One-sided empirical p with Phipson & Smyth +1 correction."""
    n = len(null_arr)
    return (np.sum(null_arr >= observed) + 1) / (n + 1)


def build_subset(df, comparison, mark_names):
    """
    Build (X, y, label) tuples per comparison sub-mode. Returns a dict.
        all  -> 1 entry under key 'all'
        A    -> 1 entry under key 'A'
        B    -> 1 entry under key 'B'
        both -> 2 entries, 'A' and 'B'
    """
    out = {}
    if comparison == "all":
        sub = df[df["switch_cat"].isin(
            ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
        )].dropna(subset=mark_names).copy()
        sub["is_switcher"] = sub["switch_cat"].isin(
            ["B_to_A", "A_to_B"]
        ).astype(int)
        X = sub[mark_names].values
        y = sub["is_switcher"].values
        out["all"] = (X, y, "Switchers vs Stable")
        return out

    comps = ["A", "B"] if comparison == "both" else [comparison]
    for comp in comps:
        spec = COMPARISON_SPEC[comp]
        sub = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])] \
                .dropna(subset=mark_names).copy()
        sub["is_switcher"] = (sub["switch_cat"] == spec["pos"]).astype(int)
        X = sub[mark_names].values
        y = sub["is_switcher"].values
        out[comp] = (X, y, spec["label"])
    return out


def run_permutation(X, y, n_perm, n_jobs, seed, label):
    """Run observed CV + permutation null. Returns dict of arrays + stats."""
    print(f"\n[INFO] === Comparison: {label} ===")
    print(f"[INFO] N={len(X)}, positive={y.sum()}, negative={(y == 0).sum()}, "
          f"prevalence={y.mean():.3f}")

    print("[INFO] Computing OBSERVED 5-fold CV AUCs...")
    t0 = time.time()
    obs_roc, obs_pr = cv_aucs(X, y, n_jobs, seed=seed)
    print(f"[INFO] Observed ROC-AUC = {obs_roc:.4f}")
    print(f"[INFO] Observed PR-AUC  = {obs_pr:.4f}")
    print(f"[INFO] Observed run took {(time.time() - t0):.1f} s")

    print(f"[INFO] Running {n_perm} label permutations...")
    rng = np.random.default_rng(seed)
    null_roc = np.zeros(n_perm)
    null_pr  = np.zeros(n_perm)
    t_start = time.time()
    for i in range(n_perm):
        y_shuffled = rng.permutation(y)
        roc_i, pr_i = cv_aucs(X, y_shuffled, n_jobs,
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


def save_tables(res, mode_key, outdir, n_perm):
    """Per-mode permutation results + summary tables."""
    perm_df = pd.DataFrame({
        "permutation":  np.arange(1, n_perm + 1),
        "null_ROC_AUC": res["null_roc"],
        "null_PR_AUC":  res["null_pr"],
    })
    perm_path = os.path.join(outdir, f"Permutation_results_{mode_key}.txt")
    perm_df.to_csv(perm_path, sep="\t", index=False)
    print(f"  Saved -> {perm_path}")

    summary = pd.DataFrame([
        {"metric": "ROC-AUC", "comparison": res["label"],
         "observed": res["obs_roc"],
         "null_mean": res["null_roc"].mean(),
         "null_sd":   res["null_roc"].std(),
         "null_min":  res["null_roc"].min(),
         "null_max":  res["null_roc"].max(),
         "n_perm":    n_perm,
         "empirical_p_one_sided": res["p_roc"]},
        {"metric": "PR-AUC",  "comparison": res["label"],
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
    """Helper: one null-distribution histogram panel."""
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
    ax.set_title(
        f"{metric_label} Null Distribution\n"
        f"empirical p = {p_emp:.3g} (n = {n_perm})",
        fontweight="bold", pad=8,
    )
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def plot_nulls(results_by_mode, outdir, comparison, n_perm):
    """
    Mode A/B/all: 1x2 figure (ROC null + PR null).
    Mode both:    2x2 figure -> top row A, bottom row B.
    """
    if comparison in ("A", "B", "all"):
        res = results_by_mode[comparison]
        fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 4.5))
        _draw_null_panel(ax_roc, res["null_roc"], res["obs_roc"], 0.5,
                          res["p_roc"], n_perm, "ROC-AUC", "Random (0.5)")
        _draw_null_panel(ax_pr, res["null_pr"],  res["obs_pr"],
                          res["prevalence"], res["p_pr"], n_perm,
                          "PR-AUC",
                          f"Random ({res['prevalence']:.3f})")
        fig.suptitle(
            f"Permutation Test: {res['label']}",
            fontweight="bold", fontsize=12, y=1.02,
        )
        plt.tight_layout()
        for ext in ("pdf", "png"):
            out = os.path.join(
                outdir, f"Fig1c_permutation_null_{comparison}.{ext}"
            )
            fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
            print(f"  Saved -> {out}")
        plt.close(fig)
        return

    # both -> 2x2
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for r_idx, comp in enumerate(["A", "B"]):
        res = results_by_mode[comp]
        _draw_null_panel(axes[r_idx, 0], res["null_roc"], res["obs_roc"], 0.5,
                          res["p_roc"], n_perm, "ROC-AUC", "Random (0.5)")
        _draw_null_panel(axes[r_idx, 1], res["null_pr"],  res["obs_pr"],
                          res["prevalence"], res["p_pr"], n_perm,
                          "PR-AUC",
                          f"Random ({res['prevalence']:.3f})")
        axes[r_idx, 0].set_title(
            axes[r_idx, 0].get_title() + f"\n[{res['label']}]",
            fontweight="bold", pad=8,
        )
        axes[r_idx, 1].set_title(
            axes[r_idx, 1].get_title() + f"\n[{res['label']}]",
            fontweight="bold", pad=8,
        )
    fig.suptitle(
        "Permutation Test: Parental Marks vs Compartment Switching"
        "\nCompartment-resolved (A on top, B on bottom)",
        fontweight="bold", fontsize=12, y=1.00,
    )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"Fig1c_permutation_null_both.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"[INFO] Loading {args.classified_file} ...")
    df = pd.read_csv(args.classified_file, sep="\t", low_memory=False)

    candidate_marks = ["H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3", "CTCF"]
    mark_names = [m for m in candidate_marks if m in df.columns]
    if len(mark_names) < 2:
        sys.exit("[ERROR] Need at least 2 mark columns in classified file.")
    print(f"[INFO] Marks: {mark_names}")
    print(f"[INFO] Comparison mode: {args.comparison}")

    subsets = build_subset(df, args.comparison, mark_names)
    if not subsets:
        sys.exit("[ERROR] No subsets built; check --comparison.")

    results_by_mode = {}
    for key, (X, y, label) in subsets.items():
        if len(X) < 50 or y.sum() < 5:
            print(f"[WARN] Insufficient data for permutation ({key}). "
                  f"Skipping.")
            continue
        res = run_permutation(X, y, args.n_perm, args.n_jobs,
                               args.seed, label)
        results_by_mode[key] = res
        save_tables(res, key, args.outdir, args.n_perm)

    if not results_by_mode:
        sys.exit("[ERROR] No permutation results produced.")

    print("\n[INFO] Plotting null distributions...")
    plot_nulls(results_by_mode, args.outdir, args.comparison, args.n_perm)

    print(f"\n[INFO] Permutation test complete -> "
          f"{os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()