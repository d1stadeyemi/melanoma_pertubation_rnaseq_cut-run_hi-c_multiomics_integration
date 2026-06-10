"""
confirm_pc1_concordance.py
Check that Parental_orig ~ Parental_DMSO and Bottom10_orig ~ Bottom10_DMSO
are concordant measurements of the same compartment state (correlated, same
sign/compartment call, no large systematic bias) — without expecting identical
values.

Usage:
  python confirm_pc1_concordance.py --pc1 data/processed/A375-CompartmentsAcrossConditions.txt
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


def concordance(a, b, name_a, name_b):
    """Run the standard agreement checks on two continuous vectors."""
    s = pd.DataFrame({"a": pd.to_numeric(a, errors="coerce"),
                      "b": pd.to_numeric(b, errors="coerce")}).dropna()
    a, b = s["a"].values, s["b"].values
    n = len(a)

    # 1. correlation (track together?)
    pear = pearsonr(a, b)[0]
    spear = spearmanr(a, b)[0]

    # 2. sign / compartment-call agreement (call the same A/B?)
    same_sign = np.mean(np.sign(a) == np.sign(b))
    # cross-tab of the disagreements
    a_pos_b_neg = int(np.sum((a > 0) & (b < 0)))
    a_neg_b_pos = int(np.sum((a < 0) & (b > 0)))

    # 3. bias / systematic offset (Bland-Altman + regression)
    diff = a - b
    mean_diff = float(np.mean(diff))          # systematic bias (should be ~0)
    sd_diff = float(np.std(diff, ddof=1))
    loa_lo, loa_hi = mean_diff - 1.96*sd_diff, mean_diff + 1.96*sd_diff
    slope, intercept = np.polyfit(b, a, 1)     # a ~ slope*b + intercept; want ~1, ~0

    print(f"\n=== {name_a}  vs  {name_b}   (n={n:,}) ===")
    print(f"  Pearson r        : {pear:.4f}   (track together; want high)")
    print(f"  Spearman rho     : {spear:.4f}")
    print(f"  Same-sign (A/B)  : {same_sign:.1%}   (call same compartment; want >~95%)")
    print(f"     disagreements : {name_a}>0 & {name_b}<0 : {a_pos_b_neg:,};  "
          f"{name_a}<0 & {name_b}>0 : {a_neg_b_pos:,}")
    print(f"  Mean difference  : {mean_diff:+.5f}  (systematic bias; want ~0)")
    print(f"  95% limits agree : [{loa_lo:+.4f}, {loa_hi:+.4f}]  (Bland-Altman LoA)")
    print(f"  Regression slope : {slope:.3f}   intercept: {intercept:+.5f}  "
          f"(want slope~1, intercept~0)")
    return dict(name=f"{name_a}_vs_{name_b}", n=n, pearson=pear, spearman=spear,
                same_sign=same_sign, mean_diff=mean_diff,
                slope=slope, intercept=intercept,
                a=a, b=b, diff=diff, mean_d=mean_diff, loa=(loa_lo, loa_hi),
                la=name_a, lb=name_b)

def make_plots(results, outpath):
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for col, r in enumerate(results):
        a, b, diff = r["a"], r["b"], r["diff"]
        # scatter with y=x
        ax = axes[0, col]
        lim = np.nanpercentile(np.abs(np.concatenate([a, b])), 99.5)
        ax.scatter(b, a, s=6, alpha=0.3, color="#1A5276", edgecolors="none")
        ax.plot([-lim, lim], [-lim, lim], "r--", lw=1, label="y = x")
        ax.axhline(0, color="#aaa", lw=0.6); ax.axvline(0, color="#aaa", lw=0.6)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel(r["lb"]); ax.set_ylabel(r["la"])
        ax.set_title(f"{r['la']} vs {r['lb']}\nPearson r={r['pearson']:.3f}, "
                     f"same-sign={r['same_sign']:.1%}", fontsize=10, fontweight="bold")
        ax.legend(frameon=False, fontsize=8)
        # Bland-Altman
        ax2 = axes[1, col]
        mean_ab = (a + b) / 2
        ax2.scatter(mean_ab, diff, s=6, alpha=0.3, color="#27AE60", edgecolors="none")
        ax2.axhline(r["mean_d"], color="black", lw=1, label=f"mean diff={r['mean_d']:+.4f}")
        ax2.axhline(r["loa"][0], color="#888", ls="--", lw=0.8, label="95% LoA")
        ax2.axhline(r["loa"][1], color="#888", ls="--", lw=0.8)
        ax2.set_xlabel(f"mean of the two"); ax2.set_ylabel(f"{r['la']} − {r['lb']}")
        ax2.set_title("Bland–Altman (agreement)", fontsize=10, fontweight="bold")
        ax2.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"\n[INFO] Plots -> {outpath}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pc1", required=True)
    ap.add_argument("--out_plot", default="pc1_concordance.png")
    ap.add_argument("--parental_orig", default="Parental_orig")
    ap.add_argument("--parental_dmso", default="Parental_DMSO")
    ap.add_argument("--b10_orig", default="Bottom10_orig")
    ap.add_argument("--b10_dmso", default="Bottom10-DMSO_fixed_inverted")
    args = ap.parse_args()

    df = pd.read_csv(args.pc1, sep="\t", low_memory=False)
    print(f"[INFO] {len(df):,} bins x {len(df.columns)} columns")
    for c in (args.parental_orig, args.parental_dmso, args.b10_orig, args.b10_dmso):
        if c not in df.columns:
            raise SystemExit(f"[ERROR] '{c}' not in columns:\n  " + "\n  ".join(df.columns))

    results = [
        concordance(df[args.parental_orig], df[args.parental_dmso],
                    args.parental_orig, args.parental_dmso),
        concordance(df[args.b10_orig], df[args.b10_dmso],
                    args.b10_orig, args.b10_dmso),
    ]
    make_plots(results, args.out_plot)

    print("\n[VERDICT] Concordant if: Pearson high (>~0.8), same-sign high "
          "(>~95%), mean difference near 0, slope near 1. High correlation with "
          "low sign-agreement would mean they track but call different "
          "compartments — investigate. A large mean difference means a "
          "systematic offset (possible scaling/inversion issue).")


if __name__ == "__main__":
    main()