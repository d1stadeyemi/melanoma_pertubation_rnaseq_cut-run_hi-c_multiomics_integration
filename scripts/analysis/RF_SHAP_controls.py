"""
RF_SHAP_controls.py — adversarial controls for the chromatin-priming claim.

PURPOSE ---------------------------------------------------------------------
  RF_SHAP.py makes the POSITIVE claim ("parental histone marks predict
  compartment switching"). This script's only job is to TRY TO BREAK that
  claim, by ruling out the two alternative explanations a reviewer is most
  likely to raise:

  CONTROL #1 — "It's not the marks, it's just the starting compartment."
    Marks correlate with baseline compartment (A regions are H3K4me1-high,
    etc.), and switching is constrained by where a bin starts. So a model
    could score well by merely reading out baseline compartment strength
    (PC1), not chromatin content. We test this with four feature variants
    evaluated on the SAME folds (so differences are paired/comparable):

        marks_only      : parental marks (the headline feature set)
        pc1_only        : baseline PC1 magnitude alone (Parental_orig)
        marks_plus_pc1  : marks + baseline PC1
        scrambled       : marks shuffled across bins (must collapse to ~0.5)

    Key reads:
      - marks_only >> pc1_only            -> marks carry real, independent signal
      - marks_plus_pc1 ~ pc1_only         -> marks add nothing beyond compartment
      - marks_plus_pc1 >> pc1_only        -> marks add information beyond compartment
      - scrambled ~ 0.5                   -> the model uses chromatin content,
                                             not bin structure / class frequency
    NOTE: the A/B compartment-resolved design (Stable_A vs A->B; Stable_B vs
    B->A) ALREADY holds the *starting compartment sign* constant within each
    comparison. pc1_only here tests the residual effect of compartment
    *strength* within that stratum — the quantitative version of #1.

  CONTROL #2 — "Your CV is optimistic because adjacent genomic bins leak."
    Neighbouring 250kb bins share compartment and marks, so random 5-fold CV
    can place a bin's neighbour in training and the bin itself in test. We
    re-evaluate every variant under CHROMOSOME-HELD-OUT CV (whole chromosomes
    never span train/test), via StratifiedGroupKFold grouped on chromosome.
      - standard CV  = optimistic upper bound
      - chromosome CV = honest cross-genome generalization
    Report BOTH. A large drop standard->chromosome flags leakage.

  (Permutation/label-shuffle null, control #1, is NOT here — already done at
   n=1000 with a clean result.)

DESIGN NOTES ----------------------------------------------------------------
  - Reuses RF_SHAP.load_data so the join / feature detection logic is identical
    to the main analysis (no divergence, no duplicated bugs).
  - Baseline PC1 = parental eigenvalue column (default 'Parental_orig'),
    already carried through load_data's merged subset.
  - All variants in a comparison are restricted to rows with a valid PC1 value
    so every variant scores the EXACT same rows (paired ΔAUROC is meaningful).
  - Folds are precomputed ONCE per (comparison, cv_scheme) and reused across
    all four variants -> per-fold deltas are paired.

USAGE -----------------------------------------------------------------------
    python scripts/RF_SHAP_controls.py \\
        data/processed/cut_run_signal_table.txt \\
        data/processed/merged_hic_rnaseq_cutrun_parent_bottom10.txt \\
        --outdir results/exploratory/RF_SHAP_controls \\
        --model rf \\
        --comparison both \\
        --include_ctcf

    # override the baseline-PC1 column or chromosome source if needed:
    python scripts/RF_SHAP_controls.py ... --pc1_col Parental_orig --chrom_col region_id

OUTPUTS (into --outdir) -----------------------------------------------------
    Controls_summary_<comp>.txt    variant x cv_scheme grid: mean/std AUROC,
                                    per-fold AUROCs, n rows, prevalence
    Controls_delta_<comp>.txt      paired ΔAUROC for the key contrasts
    Controls_pc1_diagnostic.txt    point-biserial corr(baseline PC1, switcher)
                                    within each comparison (the #1 quick check)
    Controls_AUROC_<comp>.pdf/png  grouped-bar figure (variant x cv_scheme)
    Controls_summary_all.txt       everything concatenated
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import pointbiserialr

# --- Reuse the main pipeline's loader + model builder (DRY) ------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
try:
    from RF_SHAP import (
        load_data, build_model, COMPARISON_SPEC, DEFAULT_MARKS, DEFAULT_STATS,
    )
except Exception as e:  # pragma: no cover
    sys.exit(f"[ERROR] Could not import from RF_SHAP.py (must sit beside this "
             f"script): {e}")

VARIANTS = ("marks_only", "pc1_only", "marks_plus_pc1", "scrambled")
PC1_FEATNAME = "baseline_pc1"


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Adversarial controls (#1 baseline-compartment, #2 "
                    "chromosome-held-out CV) for the priming claim.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("cutrun_file")
    p.add_argument("merged_file")
    p.add_argument("--outdir", default="results/exploratory/RF_SHAP_controls")
    p.add_argument("--model", choices=["rf", "xgboost"], default="rf")
    p.add_argument("--comparison", choices=["A", "B", "both"], default="both",
                   help="A=Stable_A vs A->B, B=Stable_B vs B->A, both=A+B. "
                        "('all' is intentionally unsupported here: the #1 "
                        "control is only interpretable within a single "
                        "starting compartment.)")

    # join / feature passthrough to load_data
    p.add_argument("--cutrun_key", default=None)
    p.add_argument("--hic_key", default=None)
    p.add_argument("--row_id", default="region_id")
    p.add_argument("--min_overlap", type=float, default=0.50)
    p.add_argument("--marks", nargs="+", default=list(DEFAULT_MARKS))
    p.add_argument("--stats", nargs="+", default=None,
                   choices=list(DEFAULT_STATS))
    p.add_argument("--include_ctcf", action="store_true")

    # control-specific
    p.add_argument("--pc1_col", default="Parental_orig",
                   help="Column holding the baseline (parental) PC1/eigenvalue "
                        "used as the compartment-strength feature for the #1 "
                        "control [default: Parental_orig].")
    p.add_argument("--chrom_col", default=None,
                   help="Column to parse chromosome from for chromosome-held-out "
                        "CV. Default: auto-detect a coordinate/region column "
                        "whose values look like 'chrN...'.")
    p.add_argument("--n_splits", type=int, default=5,
                   help="CV folds for both schemes [default 5].")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# =============================================================================
# Helpers
# =============================================================================
def _parse_chrom(series):
    """Extract a chromosome token (chr1, chrX, ...) from a coordinate/region
    string column. Returns a pandas Series of chromosome labels (or NaN)."""
    s = series.astype(str).str.strip()
    # leading 'chr' + digits/letters, tolerant of separators ':', '_', '-', tab
    chrom = s.str.extract(r"(?i)^(chr[0-9xymt]+)", expand=False)
    # fallback: bare leading token like '1', 'X' -> prefix chr
    bare = s.str.extract(r"^([0-9]{1,2}|[XYMxym])[\:\_\-\t]", expand=False)
    chrom = chrom.fillna(bare.where(bare.isna(), "chr" + bare.astype(str)))
    return chrom


def find_chrom_source(df, preferred):
    """Pick a column to parse chromosome from. Try `preferred` first, then any
    object column whose values mostly look like 'chr...'. Returns column name
    or None."""
    cand_order = []
    if preferred and preferred in df.columns:
        cand_order.append(preferred)
    # common coordinate-like names
    for c in ("region_id", "coordinate", "coordinates", "coord", "coords",
              "extra_1"):
        if c in df.columns and c not in cand_order:
            cand_order.append(c)
    # any other object column as last resort
    for c in df.columns:
        if df[c].dtype == object and c not in cand_order:
            cand_order.append(c)

    for c in cand_order:
        chrom = _parse_chrom(df[c])
        frac = chrom.notna().mean()
        if frac >= 0.80:
            return c, chrom
    return None, None


def make_folds(splitter, X_ref, y, groups):
    """Precompute fold (train, test) index pairs ONCE so all feature variants
    are scored on identical splits (paired ΔAUROC)."""
    if groups is None:
        return list(splitter.split(X_ref, y))
    return list(splitter.split(X_ref, y, groups))


def eval_variant(X, y, folds, model_name, seed):
    """Per-fold AUROC for one feature matrix on fixed folds."""
    aucs = []
    for tr, te in folds:
        y_tr = y[tr]
        # skip degenerate folds (a chromosome split can yield one-class test)
        if len(np.unique(y_tr)) < 2 or len(np.unique(y[te])) < 2:
            aucs.append(np.nan)
            continue
        if model_name == "xgboost":
            n_pos = max(int(y_tr.sum()), 1)
            n_neg = int((y_tr == 0).sum())
            model = build_model(model_name, seed, scale_pos_weight=n_neg / n_pos)
        else:
            model = build_model(model_name, seed)
        model.fit(X[tr], y_tr)
        prob = model.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], prob))
    return np.array(aucs, dtype=float)


def build_variant_matrices(sub, feat_names, pc1_name, seed):
    """Return {variant: X_matrix} for one comparison subset, all on the same
    rows. `scrambled` permutes each mark column independently."""
    rng = np.random.default_rng(seed)
    marks_X = sub[feat_names].values.astype(float)
    pc1_X = sub[[pc1_name]].values.astype(float)

    scrambled = marks_X.copy()
    for j in range(scrambled.shape[1]):
        scrambled[:, j] = rng.permutation(scrambled[:, j])

    return {
        "marks_only":     marks_X,
        "pc1_only":       pc1_X,
        "marks_plus_pc1": np.hstack([marks_X, pc1_X]),
        "scrambled":      scrambled,
    }


# =============================================================================
# Core per-comparison routine
# =============================================================================
def run_comparison(df, feat_names, comp, args, chrom_series):
    spec = COMPARISON_SPEC[comp]
    sub = df[df["switch_cat"].isin([spec["neg"], spec["pos"]])].copy()
    sub["is_switcher"] = (sub["switch_cat"] == spec["pos"]).astype(int)

    # baseline PC1 feature; restrict to rows with a valid value so every
    # variant scores identical rows (paired comparison requirement).
    if args.pc1_col not in sub.columns:
        sys.exit(f"[ERROR] --pc1_col '{args.pc1_col}' not present after load. "
                 f"Available eigen-like cols: "
                 f"{[c for c in sub.columns if 'orig' in c.lower()]}")
    sub[PC1_FEATNAME] = pd.to_numeric(sub[args.pc1_col], errors="coerce")
    n_before = len(sub)
    sub = sub[sub[PC1_FEATNAME].notna()].copy()
    # align chromosome labels to the surviving rows
    chrom = chrom_series.loc[sub.index] if chrom_series is not None else None
    n_after = len(sub)
    if n_after < n_before:
        print(f"  [INFO] {comp}: dropped {n_before - n_after} rows lacking a "
              f"valid {args.pc1_col}; {n_after} retained (paired across variants).")

    y = sub["is_switcher"].values
    if len(sub) < 50 or y.sum() < 5 or (y == 0).sum() < 5:
        print(f"  [WARN] {comp}: insufficient data "
              f"(N={len(sub)}, pos={int(y.sum())}). Skipping.")
        return None

    prevalence = float(y.mean())
    print(f"  [INFO] {comp} ({spec['label']}): N={len(sub)}, "
          f"pos={int(y.sum())}, neg={int((y==0).sum())}, "
          f"prevalence={prevalence:.3f}")

    # --- #1 quick diagnostic: does baseline PC1 alone correlate w/ switching? -
    try:
        r_pb, p_pb = pointbiserialr(y, sub[PC1_FEATNAME].values)
    except Exception:
        r_pb, p_pb = np.nan, np.nan
    print(f"  [#1 diag] {comp}: point-biserial corr(baseline PC1, switcher) "
          f"r={r_pb:+.3f}, p={p_pb:.2e}")

    Xvar = build_variant_matrices(sub, feat_names, PC1_FEATNAME, args.seed)

    # --- CV schemes ----------------------------------------------------------
    schemes = {}
    schemes["standard"] = (
        StratifiedKFold(n_splits=args.n_splits, shuffle=True,
                        random_state=args.seed),
        None,
    )
    groups = None
    if chrom is not None:
        groups = chrom.values
        n_chrom = pd.Series(groups).nunique()
        if n_chrom < args.n_splits:
            print(f"  [WARN] {comp}: only {n_chrom} chromosomes < n_splits "
                  f"{args.n_splits}; chromosome CV may be unstable.")
        schemes["chromosome"] = (
            StratifiedGroupKFold(n_splits=args.n_splits, shuffle=True,
                                 random_state=args.seed),
            groups,
        )
    else:
        print(f"  [WARN] {comp}: no chromosome source -> skipping chromosome CV.")

    # reference matrix for fold generation (rows are identical across variants)
    X_ref = Xvar["marks_only"]

    rows = []   # for summary table
    auc_store = {}  # auc_store[scheme][variant] = per-fold array
    for scheme_name, (splitter, grp) in schemes.items():
        folds = make_folds(splitter, X_ref, y, grp)
        auc_store[scheme_name] = {}
        for v in VARIANTS:
            aucs = eval_variant(Xvar[v], y, folds, args.model, args.seed)
            auc_store[scheme_name][v] = aucs
            rows.append({
                "comparison": spec["short"],
                "comparison_label": spec["label"],
                "cv_scheme": scheme_name,
                "variant": v,
                "n_rows": len(sub),
                "prevalence": round(prevalence, 4),
                "n_features": Xvar[v].shape[1],
                "mean_auroc": round(np.nanmean(aucs), 4),
                "std_auroc": round(np.nanstd(aucs), 4),
                "per_fold_auroc": ";".join(
                    f"{a:.4f}" if np.isfinite(a) else "NA" for a in aucs),
                "n_valid_folds": int(np.isfinite(aucs).sum()),
            })
            print(f"    [{scheme_name:10s}] {v:15s} "
                  f"AUROC={np.nanmean(aucs):.3f} ± {np.nanstd(aucs):.3f} "
                  f"(folds={int(np.isfinite(aucs).sum())})")

    summary_df = pd.DataFrame(rows)

    # --- paired ΔAUROC contrasts (same folds within a scheme) ---------------
    delta_rows = []
    contrasts = [
        ("marks_only", "pc1_only",
         "marks carry signal beyond compartment STRENGTH"),
        ("marks_plus_pc1", "pc1_only",
         "marks ADD information on top of compartment STRENGTH"),
        ("marks_plus_pc1", "marks_only",
         "compartment STRENGTH adds on top of marks"),
        ("marks_only", "scrambled",
         "real chromatin content vs shuffled (artifact check)"),
    ]
    for scheme_name in auc_store:
        for a, b, meaning in contrasts:
            da = auc_store[scheme_name][a]
            db = auc_store[scheme_name][b]
            paired = da - db                      # per-fold paired deltas
            valid = np.isfinite(paired)
            mean_d = np.nanmean(paired)
            n_a_gt_b = int(np.sum((paired > 0) & valid))
            delta_rows.append({
                "comparison": spec["short"],
                "cv_scheme": scheme_name,
                "contrast": f"{a} - {b}",
                "meaning": meaning,
                "mean_delta_auroc": round(mean_d, 4),
                "per_fold_delta": ";".join(
                    f"{d:+.4f}" if np.isfinite(d) else "NA" for d in paired),
                "folds_A_gt_B": f"{n_a_gt_b}/{int(valid.sum())}",
            })
    delta_df = pd.DataFrame(delta_rows)

    # --- figure --------------------------------------------------------------
    _plot_grid(summary_df, spec, args.outdir, args.model)

    return {
        "summary": summary_df,
        "delta": delta_df,
        "pc1_diag": {"comparison": spec["short"], "label": spec["label"],
                     "pointbiserial_r": round(r_pb, 4) if np.isfinite(r_pb) else np.nan,
                     "pointbiserial_p": p_pb, "n": len(sub),
                     "prevalence": round(prevalence, 4)},
    }


# =============================================================================
# Figure
# =============================================================================
def _plot_grid(summary_df, spec, outdir, model_name):
    os.makedirs(outdir, exist_ok=True)
    schemes = list(dict.fromkeys(summary_df["cv_scheme"]))  # preserve order
    colors = {"standard": "#1A5276", "chromosome": "#C0392B"}
    width = 0.38

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    xbase = np.arange(len(VARIANTS))
    for k, scheme in enumerate(schemes):
        means, stds = [], []
        for v in VARIANTS:
            r = summary_df[(summary_df["cv_scheme"] == scheme) &
                           (summary_df["variant"] == v)]
            means.append(float(r["mean_auroc"].iloc[0]) if len(r) else np.nan)
            stds.append(float(r["std_auroc"].iloc[0]) if len(r) else 0.0)
        off = (k - (len(schemes) - 1) / 2) * width
        bars = ax.bar(xbase + off, means, width, yerr=stds, capsize=3,
                      color=colors.get(scheme, "#888"), alpha=0.9,
                      edgecolor="white", linewidth=0.6,
                      label=f"{scheme} CV")
        for bar, m in zip(bars, means):
            if np.isfinite(m):
                ax.text(bar.get_x() + bar.get_width() / 2, m + 0.012,
                        f"{m:.2f}", ha="center", va="bottom", fontsize=8)

    ax.axhline(0.5, color="black", lw=0.8, ls="--", alpha=0.6)
    ax.text(len(VARIANTS) - 0.5, 0.505, "chance", fontsize=8, color="black",
            ha="right", va="bottom")
    ax.set_xticks(xbase)
    ax.set_xticklabels(["marks\nonly", "PC1\nonly", "marks\n+PC1",
                        "scrambled"], fontsize=9)
    ax.set_ylabel("AUROC (mean ± SD across folds)")
    ax.set_ylim(0.4, 1.0)
    ax.set_title(f"Control battery — {spec['label']}\n"
                 f"({'Random Forest' if model_name=='rf' else 'XGBoost'})",
                 fontweight="bold", fontsize=11, pad=8)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    plt.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"Controls_AUROC_{spec['short']}.{ext}")
        fig.savefig(path, bbox_inches="tight", dpi=600)
        print(f"  Saved -> {path}")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    print(f"[INFO] Controls output dir: {args.outdir}")
    print(f"[INFO] Model backend: {args.model}")

    stats = list(DEFAULT_STATS) if args.stats is None else args.stats

    # Load with the SAME loader as the main analysis; baseline features only.
    df, feat_names, feat_map, used_row_id = load_data(
        args.cutrun_file, args.merged_file,
        marks=args.marks, stats=stats, include_ctcf=args.include_ctcf,
        cutrun_key_override=args.cutrun_key, hic_key_override=args.hic_key,
        row_id_col=args.row_id, min_overlap_frac=args.min_overlap,
        feature_mode="baseline",
    )

    # Chromosome source for control #2.
    chrom_col, chrom_series = find_chrom_source(
        df, preferred=args.chrom_col or used_row_id)
    if chrom_col is not None:
        n_chrom = chrom_series.nunique()
        print(f"[INFO] Chromosome source for held-out CV: '{chrom_col}' "
              f"({n_chrom} chromosomes detected).")
    else:
        print("[WARN] No chromosome-like column found; chromosome-held-out CV "
              "(#2) will be skipped. Pass --chrom_col explicitly to enable.")

    comps = ["A", "B"] if args.comparison == "both" else [args.comparison]

    all_summ, all_delta, pc1_diag = [], [], []
    for comp in comps:
        print(f"\n[INFO] === Comparison {comp} ===")
        res = run_comparison(df, feat_names, comp, args, chrom_series)
        if res is None:
            continue
        # write per-comparison tables
        s_path = os.path.join(args.outdir, f"Controls_summary_{comp}.txt")
        d_path = os.path.join(args.outdir, f"Controls_delta_{comp}.txt")
        res["summary"].to_csv(s_path, sep="\t", index=False)
        res["delta"].to_csv(d_path, sep="\t", index=False)
        print(f"  Summary -> {s_path}")
        print(f"  Deltas  -> {d_path}")
        all_summ.append(res["summary"])
        all_delta.append(res["delta"])
        pc1_diag.append(res["pc1_diag"])

    if not all_summ:
        sys.exit("[ERROR] No comparisons produced results.")

    pd.concat(all_summ, ignore_index=True).to_csv(
        os.path.join(args.outdir, "Controls_summary_all.txt"),
        sep="\t", index=False)
    pd.concat(all_delta, ignore_index=True).to_csv(
        os.path.join(args.outdir, "Controls_delta_all.txt"),
        sep="\t", index=False)
    pd.DataFrame(pc1_diag).to_csv(
        os.path.join(args.outdir, "Controls_pc1_diagnostic.txt"),
        sep="\t", index=False)

    # --- console verdict ----------------------------------------------------
    print("\n" + "=" * 70)
    print("INTERPRETATION GUIDE (read against Controls_delta_all.txt):")
    print("  marks_only - pc1_only      > 0 and consistent  -> marks carry")
    print("                                 signal beyond compartment strength")
    print("  marks_plus_pc1 - pc1_only  >> 0                 -> marks ADD info")
    print("  marks_plus_pc1 - pc1_only  ~ 0                  -> marks redundant")
    print("                                 with compartment (RETHINK claim)")
    print("  scrambled                  ~ 0.5                -> no structural")
    print("                                 artifact; model uses chromatin")
    print("  standard -> chromosome AUROC drop large         -> spatial leakage;")
    print("                                 report chromosome CV as the honest #")
    print("=" * 70)
    print(f"\n[INFO] Controls complete -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()