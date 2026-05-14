"""
Permutation Test — empirical null distribution for RF ROC-AUC and PR-AUC.

Tests the null hypothesis that Parental epigenetic marks have NO relationship
to compartment switching. Shuffles labels n_perm times, refits the RF on each
shuffle, and compares the real cross-validated AUCs against the shuffled null.

Reuses the gene table produced by RF_SHAP.py (Switch_classified.txt) so the
exact same genes, same labels, and same features go into both analyses.

Usage:
    python scripts/RF_permutation_test.py \
        results/exploratory/RF_SHAP/Switch_classified.txt \
        --outdir results/exploratory/RF_permutation_test \
        --n_perm 100 \
        --n_jobs 8

Outputs (in --outdir):
    Fig1c_permutation_null.pdf/png   ROC-AUC + PR-AUC null vs observed
    Permutation_results.txt          per-permutation AUCs + empirical p-values
    Permutation_summary.txt          observed AUC, null mean/SD, empirical p

Runtime guidance (5-fold CV, 500-tree RF, 4 marks, ~21k genes):
    n_perm = 100   ~1.0-1.5 hr on 8 cores  (rotation report)
    n_perm = 1000  ~10-15 hr on 8 cores    (overnight, manuscript)

SLURM submission template:
    #SBATCH --nodes=1
    #SBATCH --time=24:00:00
    #SBATCH --ntasks-per-node=8
    #SBATCH --mem=32G
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


def parse_args():
    p = argparse.ArgumentParser(description="RF permutation test")
    p.add_argument("classified_file", default="results/exploratory/RF_SHAP/Switch_classified.txt",
                   help="Path to Switch_classified.txt from RF_SHAP.py")
    p.add_argument("--outdir", default="results/exploratory/RF_permutation_test",
                   help="Output directory")
    p.add_argument("--n_perm", type=int, default=100,
                   help="Number of label permutations (default: 100)")
    p.add_argument("--n_jobs", type=int, default=-1,
                   help="Cores for RF (-1 = all; default -1)")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for reproducibility (default: 42)")
    return p.parse_args()


def cv_aucs(X, y, n_jobs, seed):
    """
    Run a single 5-fold CV pass of the RF and return (mean ROC-AUC, mean PR-AUC).
    """
    rf = RandomForestClassifier(
        n_estimators=500, max_depth=6,
        class_weight="balanced", random_state=seed, n_jobs=n_jobs
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    roc_per_fold = []
    pr_per_fold  = []
    for train_idx, test_idx in cv.split(X, y):
        rf.fit(X[train_idx], y[train_idx])
        y_prob = rf.predict_proba(X[test_idx])[:, 1]
        roc_per_fold.append(roc_auc_score(y[test_idx], y_prob))
        pr_per_fold.append(average_precision_score(y[test_idx], y_prob))
    return float(np.mean(roc_per_fold)), float(np.mean(pr_per_fold))


def empirical_p(observed, null_arr):
    """
    One-sided empirical p-value: P(null >= observed).
    With +1 correction (Phipson & Smyth 2010) to avoid p = 0.
    """
    n = len(null_arr)
    return (np.sum(null_arr >= observed) + 1) / (n + 1)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load classified gene table from RF_SHAP.py ───────────────────────────
    print(f"[INFO] Loading {args.classified_file} ...")
    df = pd.read_csv(args.classified_file, sep="\t", low_memory=False)

    # Detect mark columns (4 standard + optional CTCF)
    candidate_marks = ["H3K27me3", "H3K27Ac", "H3K4me1", "H3K9me3", "CTCF"]
    mark_names = [m for m in candidate_marks if m in df.columns]
    if len(mark_names) < 2:
        sys.exit("[ERROR] Need at least 2 mark columns in classified file.")
    print(f"[INFO] Marks: {mark_names}")

    df = df[df["switch_cat"].isin(
        ["B_to_A", "A_to_B", "Stable_A", "Stable_B"]
    )].dropna(subset=mark_names).copy()
    df["is_switcher"] = df["switch_cat"].isin(["B_to_A", "A_to_B"]).astype(int)

    X = df[mark_names].values
    y = df["is_switcher"].values
    print(f"[INFO] N genes: {len(X):,}  Switchers: {y.sum():,}  "
          f"Stable: {(y == 0).sum():,}")
    print(f"[INFO] Positive prevalence: {y.mean():.3f}")

    # ── 1. Observed AUCs ─────────────────────────────────────────────────────
    print("\n[INFO] Computing OBSERVED 5-fold CV AUCs (no permutation)...")
    t0 = time.time()
    obs_roc, obs_pr = cv_aucs(X, y, args.n_jobs, seed=args.seed)
    print(f"[INFO] Observed ROC-AUC = {obs_roc:.4f}")
    print(f"[INFO] Observed PR-AUC  = {obs_pr:.4f}")
    print(f"[INFO] Observed run took {(time.time() - t0):.1f} s")

    # ── 2. Permutation null ──────────────────────────────────────────────────
    print(f"\n[INFO] Running {args.n_perm} label permutations...")
    rng = np.random.default_rng(args.seed)
    null_roc = np.zeros(args.n_perm)
    null_pr  = np.zeros(args.n_perm)

    t_start = time.time()
    for i in range(args.n_perm):
        y_shuffled = rng.permutation(y)
        # Use a different RF seed per permutation so the null reflects both
        # label-shuffling noise and model stochasticity.
        roc_i, pr_i = cv_aucs(X, y_shuffled, args.n_jobs,
                              seed=args.seed + i + 1)
        null_roc[i] = roc_i
        null_pr[i]  = pr_i

        # Progress every 10
        if (i + 1) % 10 == 0 or i == 0:
            elapsed = time.time() - t_start
            per_perm = elapsed / (i + 1)
            eta = per_perm * (args.n_perm - (i + 1))
            print(f"  perm {i+1:>4d}/{args.n_perm}  "
                  f"ROC={roc_i:.3f}  PR={pr_i:.3f}  "
                  f"({elapsed/60:.1f} min elapsed; ETA {eta/60:.1f} min)")

    total_min = (time.time() - t_start) / 60
    print(f"[INFO] Permutations complete in {total_min:.1f} min")

    # ── 3. Empirical p-values ────────────────────────────────────────────────
    p_roc = empirical_p(obs_roc, null_roc)
    p_pr  = empirical_p(obs_pr,  null_pr)

    print(f"\n[INFO] Null ROC-AUC: mean = {null_roc.mean():.4f}  "
          f"SD = {null_roc.std():.4f}")
    print(f"[INFO] Null PR-AUC : mean = {null_pr.mean():.4f}  "
          f"SD = {null_pr.std():.4f}")
    print(f"[INFO] Observed ROC-AUC = {obs_roc:.4f}  "
          f"empirical p = {p_roc:.4g}")
    print(f"[INFO] Observed PR-AUC  = {obs_pr:.4f}  "
          f"empirical p = {p_pr:.4g}")

    # ── 4. Save per-permutation table ────────────────────────────────────────
    perm_df = pd.DataFrame({
        "permutation": np.arange(1, args.n_perm + 1),
        "null_ROC_AUC": null_roc,
        "null_PR_AUC":  null_pr,
    })
    perm_path = os.path.join(args.outdir, "Permutation_results.txt")
    perm_df.to_csv(perm_path, sep="\t", index=False)
    print(f"[INFO] Per-permutation results saved -> {perm_path}")

    # ── 5. Save summary ──────────────────────────────────────────────────────
    summary = pd.DataFrame([
        {"metric": "ROC-AUC",
         "observed": obs_roc,
         "null_mean": null_roc.mean(),
         "null_sd":   null_roc.std(),
         "null_min":  null_roc.min(),
         "null_max":  null_roc.max(),
         "n_perm":    args.n_perm,
         "empirical_p_one_sided": p_roc},
        {"metric": "PR-AUC",
         "observed": obs_pr,
         "null_mean": null_pr.mean(),
         "null_sd":   null_pr.std(),
         "null_min":  null_pr.min(),
         "null_max":  null_pr.max(),
         "n_perm":    args.n_perm,
         "empirical_p_one_sided": p_pr},
    ])
    summary_path = os.path.join(args.outdir, "Permutation_summary.txt")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"[INFO] Summary saved -> {summary_path}")
    print(summary.to_string())

    # ── 6. Figure: null distributions + observed lines ───────────────────────
    fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 4.5))

    # ROC null
    ax_roc.hist(null_roc, bins=30, color="#85C1E9",
                edgecolor="white", linewidth=0.5, alpha=0.95)
    ax_roc.axvline(obs_roc, color="#C0392B", lw=2.0,
                   label=f"Observed = {obs_roc:.3f}")
    ax_roc.axvline(null_roc.mean(), color="black",
                   lw=1.0, linestyle="--",
                   label=f"Null mean = {null_roc.mean():.3f}")
    ax_roc.axvline(0.5, color="grey", lw=0.8, linestyle=":",
                   label="Random (0.5)")
    ax_roc.set_xlabel("ROC-AUC (5-fold CV)", labelpad=6)
    ax_roc.set_ylabel("Permutations", labelpad=6)
    ax_roc.set_title(
        f"ROC-AUC Null Distribution\nempirical p = {p_roc:.3g} "
        f"(n = {args.n_perm})",
        fontweight="bold", pad=8
    )
    ax_roc.legend(frameon=False, fontsize=9, loc="upper right")

    # PR null
    ax_pr.hist(null_pr, bins=30, color="#85C1E9",
               edgecolor="white", linewidth=0.5, alpha=0.95)
    ax_pr.axvline(obs_pr, color="#C0392B", lw=2.0,
                  label=f"Observed = {obs_pr:.3f}")
    ax_pr.axvline(null_pr.mean(), color="black",
                  lw=1.0, linestyle="--",
                  label=f"Null mean = {null_pr.mean():.3f}")
    ax_pr.axvline(y.mean(), color="grey", lw=0.8, linestyle=":",
                  label=f"Random ({y.mean():.3f})")
    ax_pr.set_xlabel("PR-AUC (5-fold CV)", labelpad=6)
    ax_pr.set_ylabel("Permutations", labelpad=6)
    ax_pr.set_title(
        f"PR-AUC Null Distribution\nempirical p = {p_pr:.3g} "
        f"(n = {args.n_perm})",
        fontweight="bold", pad=8
    )
    ax_pr.legend(frameon=False, fontsize=9, loc="upper right")

    fig.suptitle(
        "Permutation Test: Parental Marks vs Compartment Switching",
        fontweight="bold", fontsize=12, y=1.02
    )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(args.outdir, f"Fig1c_permutation_null.{ext}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05, dpi=600)
        print(f"  Saved -> {out}")
    plt.close(fig)

    print(f"\n[INFO] Permutation test complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()