#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unc_selective_erasure.py
===========================
Does UNC1999 (EZH2/EZH1 inhibitor) SELECTIVELY erase the H3K27me3 that
constricted migration installed?

Scientific core unchanged from prior versions. This revision keeps everything as
before and ADDS a COMPARTMENT-MATCHED reference analysis on top of the existing
outputs.

WHY compartment matching (read this):
    Diagnostic D5 shows dUNC tracks compartment (PC1): in the Parental arm dUNC is a
    steep function of PC1, in B10 it is flat. The switcher / reverter / gainer sets
    are compartment-PURE, but the generic "stable" reference is compartment-MIXED, so
    the naive target-vs-stable delta partly re-reads the dUNC~PC1 gradient rather than
    a set-specific effect (hence Parental over-shooting the stable median instead of
    hugging it). The fix is to contrast each target ONLY against stable bins of the
    SAME baseline compartment IN THE ARM BEING MEASURED. A switcher changes compartment
    between cell lines, so the matching flips per arm:
        A_to_B locus            = A in Parental , B in B10
        B_to_A locus            = B in Parental , A in B10
        reverter rev_atb (UNC A->B, starts A)   = A in B10 baseline
        reverter rev_bta (UNC B->A, starts B)   = B in B10 baseline
    This is a SELECTIVITY test and is valid under CPM (it is within-arm & within-
    compartment); it does NOT license cross-cell-line magnitude claims, which still
    require spike-in / greenlist normalization.

HOW TO READ THE DELTA (sign convention -- read before interpreting any panel):
    delta = Cliff's delta of (target dUNC) vs (same-compartment stable dUNC), same arm.
      delta < 0  -> target sits BELOW its stable reference -> loses MORE mark than its
                    same-compartment neighbours  ("extra loss" / preferentially erased)
      delta > 0  -> target sits ABOVE its stable reference -> loses LESS mark than its
                    same-compartment neighbours  ("spared" / relatively resistant)
    delta is RANK-based and each normalization applies ONE global factor per sample
    (i.e. one additive shift to every dUNC in that arm), so a global shift moves target
    and reference together and CANCELS. Consequence: these deltas are essentially
    normalization-invariant (cpm vs greenlist vs spike-in give ~the same delta). What
    normalization DOES move is the absolute violin HEIGHT: under cpm the fixed
    per-million total makes a globally-depleted arm read as apparent "gain" (positive
    dUNC) -- an artifact. Use --norm greenlist to place the heights on a scale where
    global loss reads as loss. NEVER negate dUNC to "fix" a sign: that reverses every
    rank in the arm and flips ALL deltas in that arm at once (it would turn the gainer
    depletion into an apparent gain), i.e. it destroys the result rather than correcting it.

This revision also retains the anomaly diagnostics suite (D1-D5) and the selectable
baseline-control mode --baseline-match {decile,fine,residual}.

New outputs (in --outdir):
    compartment_matched_results.tsv          per-contrast Cliff's delta / CI / perm-p
                                             (switchers + reverters; H3K27me3 gainers are
                                              added only with --cmatched-gainers), each vs
                                              same-compartment stable in its arm
    forest_effect_sizes_cmatched.png         forest of the compartment-matched contrasts
    violin_cmatched_b10_stableA.png          B10 arm, A-compartment sets vs Stable A
    violin_cmatched_b10_stableB.png          B10 arm, B-compartment sets vs Stable B
    violin_cmatched_parental_stableA.png     Parental arm, A-compartment set vs Stable A
    violin_cmatched_parental_stableB.png     Parental arm, B-compartment set vs Stable B
    violin_cmatched_combined_stableA.png     both arms, A-compartment (per-arm stable ref)
    violin_cmatched_combined_stableB.png     both arms, B-compartment (per-arm stable ref)
    violin_cmatched_arm_parental.png         PARENTAL only: A-compartment group + B-compartment
                                             group side by side, each vs its own stable
    violin_cmatched_arm_b10.png              B10 only: same layout (incl. reverters)

Environment: conda env multiomics-py (deeptools, numpy, scipy, pandas, matplotlib)

USAGE EXAMPLES
--------------
# Standard run (CPM), diagnostics + compartment-matched analysis both on by default:
python unc_selective_erasure.py \
    --sample-sheet sample_sheet.csv \
    --compartments compartments_250kb_hg19.tsv \
    --counts-table raw_bin_counts.tsv \
    --project-root . \
    --outdir results/exploratory/unc_selective_erasure \
    --plot-gainers

# Publication-grade: greenlist puts the violin HEIGHTS on a scale where the global
# H3K27me3 loss reads as loss (the deltas themselves barely move -- they are rank-based):
python unc_selective_erasure.py \
    --sample-sheet sample_sheet.csv --compartments comp_250kb.tsv \
    --counts-table raw_bin_counts.tsv \
    --norm greenlist --size-factors greenlist_size_factors.tsv \
    --baseline-match residual --n-perm 1000 --n-boot 1000 \
    --outdir results/final/unc_selective_erasure

# Only the 4 single-arm compartment-matched violins, skip the 2 combined panels:
python unc_selective_erasure.py ... --no-cmatched-combined

# Skip the 2 per-arm (both-compartments-on-one-axis) panels:
python unc_selective_erasure.py ... --no-cmatched-per-arm

# Add the gainer contrasts (OPT-IN: off unless this flag is given):
python unc_selective_erasure.py ... --cmatched-gainers

# Skip the compartment-matched block entirely (old behaviour):
python unc_selective_erasure.py ... --no-compartment-matched

# Validate the statistics core without any data:
python unc_selective_erasure.py --self-test
"""
import argparse, os, subprocess, sys, shutil, json, datetime as _dt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL_COORD="coordinate"; COL_CHR="chr"; COL_START="start"; COL_END="end"
PC1_PAR_DMSO="Parental_DMSO"; PC1_PAR_UNC="Parental_UNC"
PC1_B10_DMSO="Bottom10-DMSO"; PC1_B10_UNC="Bottom10-UNC"
PC1_PAR_ORIG="Parental_orig"; PC1_B10_ORIG="Bottom10_orig"
MIG_BOTH_DMSO="P->B10_Both_DMSO"; MIG_SWITCH_DMSO="P->B10_Switch_DMSO"; MIG_SHIFT_DMSO="P->B10_Shift_DMSO"
MIG_BOTH_ORIG="P->B10_Both_original"; MIG_SWITCH_ORIG="P->B10_Switch_Orig"; MIG_SHIFT_ORIG="P->B10_Shift_Orig"
UNC_BOTH_B10="B10_DMSO->UNC_Both"; UNC_SWITCH_B10="B10_DMSO->UNC_switch"; UNC_SHIFT_B10="B10_DMSO->UNC_shift"
NULL_TOKENS={"-","","nan","none","na","."}


def _assign_roles(cond, trt):
    role=[]
    for c,t in zip(cond.str.strip().str.lower(), trt.str.strip().str.lower()):
        cc=("parental" if c.startswith(("parent","ctrl","control"))
            else ("b10" if c in {"b10","bottom","bottom10"} else c))
        tt=("unc" if t.startswith("unc") else ("none" if t in {"none","dmso"} else t))
        role.append(f"{cc}_{tt}")
    return role


def load_samples(size_factor_path, sample_sheet_path, project_root, factor="H3K27me3", drop_samples=()):
    sheet=pd.read_csv(sample_sheet_path)
    for col in ("SampleID","bamReads"):
        if col not in sheet.columns: sys.exit(f"[ERROR] sample sheet needs a '{col}' column")
    sheet_small=sheet[["SampleID","bamReads"]].drop_duplicates("SampleID")
    if size_factor_path:
        sf=pd.read_csv(size_factor_path, sep="\t")
        need={"SampleID","Condition","Treatment","deseq2_size_factor"}
        if need-set(sf.columns): sys.exit(f"[ERROR] size-factor file missing columns: {need-set(sf.columns)}")
        df=sf.merge(sheet_small, on="SampleID", how="left")
    else:
        need={"SampleID","Condition","Treatment","Factor"}
        if need-set(sheet.columns): sys.exit(f"[ERROR] sample sheet missing columns: {need-set(sheet.columns)}")
        df=sheet[sheet["Factor"].str.strip().str.lower()==factor.strip().lower()].drop_duplicates("SampleID").copy()
        df["deseq2_size_factor"]=np.nan
    drop={s.strip() for s in drop_samples if s.strip()}
    if drop: df=df[~df["SampleID"].isin(drop)].copy()
    if df["bamReads"].isna().any():
        sys.exit(f"[ERROR] no BAM in sample sheet for: {df.loc[df['bamReads'].isna(),'SampleID'].tolist()}")
    df["role"]=_assign_roles(df["Condition"], df["Treatment"])
    df["bam"]=df["bamReads"].apply(lambda p: os.path.normpath(os.path.join(project_root, str(p))))
    expected={"parental_none","parental_unc","b10_none","b10_unc"}
    if set(df["role"])!=expected: sys.exit(f"[ERROR] roles {sorted(set(df['role']))} != {sorted(expected)}")
    if any(v!=2 for v in df["role"].value_counts().to_dict().values()):
        sys.exit(f"[ERROR] expected 2 replicates per role, got {df['role'].value_counts().to_dict()}")
    return df.reset_index(drop=True)


def load_bins(compartments_path, keep_chroms):
    comp=pd.read_csv(compartments_path, sep="\t", dtype={COL_CHR:str})
    for c in (COL_COORD,COL_CHR,COL_START,COL_END):
        if c not in comp.columns: sys.exit(f"[ERROR] compartment file missing column '{c}'")
    comp=comp[comp[COL_CHR].isin(keep_chroms)].copy()
    comp[COL_START]=comp[COL_START].astype(int); comp[COL_END]=comp[COL_END].astype(int)
    comp["bed_start"]=comp[COL_START]-1; comp["bed_end"]=comp[COL_END]
    return comp.sort_values([COL_CHR,"bed_start"]).reset_index(drop=True)


def write_bins_bed(comp, path):
    comp[[COL_CHR,"bed_start","bed_end",COL_COORD]].to_csv(path, sep="\t", header=False, index=False)
    return path


def compute_blacklist_fraction(bins, blacklist_path):
    bl=pd.read_csv(blacklist_path, sep="\t", header=None, usecols=[0,1,2], names=["chr","start","end"])
    bl["chr"]=bl["chr"].astype(str)
    bl_by={c:g[["start","end"]].values for c,g in bl.groupby("chr")}
    chrs=bins["chr"].astype(str).values; starts=bins["bed_start"].values; ends=bins["bed_end"].values
    frac=np.zeros(len(bins), float)
    for i in range(len(bins)):
        regions=bl_by.get(chrs[i])
        if regions is None: continue
        s,e=starts[i],ends[i]
        ov=np.clip(np.minimum(e,regions[:,1])-np.maximum(s,regions[:,0]),0,None).sum()
        frac[i]=ov/max(e-s,1)
    return frac


def count_reads(bins_bed, samples, out_npz, out_raw, blacklist, threads, min_mapq, extend_reads, use_deeptools_blacklist=False):
    if shutil.which("multiBamSummary") is None: sys.exit("[ERROR] multiBamSummary not found; activate env 'multiomics-py'")
    for b in samples["bam"]:
        if not os.path.exists(b): sys.exit(f"[ERROR] BAM not found: {b}")
        if not (os.path.exists(b+".bai") or os.path.exists(b[:-4]+".bai")): sys.exit(f"[ERROR] index (.bai) missing for: {b}")
    cmd=["multiBamSummary","BED-file","--BED",bins_bed,"--bamfiles",*list(samples["bam"]),
         "--labels",*list(samples["SampleID"]),"-o",out_npz,"--outRawCounts",out_raw,
         "-p",str(threads),"--minMappingQuality",str(min_mapq)]
    if blacklist and use_deeptools_blacklist: cmd+=["--blackListFileName",blacklist]
    if extend_reads: cmd+=["--extendReads"]
    print("[run]"," ".join(cmd), flush=True); subprocess.run(cmd, check=True)
    return out_raw


def load_raw_counts(raw_path, samples):
    raw=pd.read_csv(raw_path, sep="\t")
    raw.columns=[c.strip().strip("#").strip().strip("'\"") for c in raw.columns]
    ren={}
    for c in raw.columns:
        cl=c.lower()
        if cl in ("chr","chrom","#chr","'chr'"): ren[c]="chr"
        elif cl=="start": ren[c]="bed_start"
        elif cl=="end": ren[c]="bed_end"
    raw=raw.rename(columns=ren)
    raw["chr"]=raw["chr"].astype(str); raw["bed_start"]=raw["bed_start"].astype(int); raw["bed_end"]=raw["bed_end"].astype(int)
    if len([s for s in samples["SampleID"] if s in raw.columns])!=len(samples):
        sys.exit("[ERROR] counts matrix columns do not cover all samples")
    return raw[["chr","bed_start","bed_end"]+list(samples["SampleID"])]


def compute_divisors(raw, samples, norm):
    if norm=="greenlist":
        if samples["deseq2_size_factor"].isna().any(): sys.exit("[ERROR] --norm greenlist but size factors are missing")
        return {sid:float(sf) for sid,sf in zip(samples["SampleID"], samples["deseq2_size_factor"])}, "size_factor"
    div={}
    for sid in samples["SampleID"]:
        lib=float(raw[sid].sum())
        if lib<=0: sys.exit(f"[ERROR] sample {sid} has zero library size")
        div[sid]=lib/1e6
    return div, "libsize/1e6"


def normalize_and_aggregate(raw, samples, pseudocount, norm):
    divisor,_=compute_divisors(raw, samples, norm)
    mean_div=float(np.mean(list(divisor.values())))
    role_of=dict(zip(samples["SampleID"], samples["role"]))
    out=raw[["chr","bed_start","bed_end"]].copy()
    by_role={r:[] for r in ("parental_none","parental_unc","b10_none","b10_unc")}
    for sid in samples["SampleID"]:
        prior_raw=pseudocount*divisor[sid]/mean_div
        by_role[role_of[sid]].append((raw[sid].astype(float)+prior_raw)/divisor[sid])
    for r,cols in by_role.items():
        out[r]=np.mean(np.vstack([c.values for c in cols]), axis=0)
    return out


def categorize_direction(both, switch, shift, pc1_par, pc1_b10):
    def clean(v): return str(v).strip()
    v=clean(both)
    if v.lower() not in NULL_TOKENS:
        if "toA" in v or "BtoA" in v: return "B_to_A"
        if "toB" in v or "AtoB" in v: return "A_to_B"
    v=clean(switch)
    if v.lower() not in NULL_TOKENS:
        if v=="BtoA": return "B_to_A"
        if v=="AtoB": return "A_to_B"
    v=clean(shift)
    if v.lower() not in NULL_TOKENS:
        if "TowardA" in v: return "B_to_A"
        if "TowardB" in v: return "A_to_B"
    try:
        p,b=float(pc1_par),float(pc1_b10)
        if p>0 and b>0: return "Stable_A"
        if p<0 and b<0: return "Stable_B"
    except (ValueError,TypeError): pass
    return None


def categorize_movement(both, switch, shift):
    s=str(both).strip()
    if s.lower() not in NULL_TOKENS:
        if "toA" in s or "BtoA" in s: return "B_to_A"
        if "toB" in s or "AtoB" in s: return "A_to_B"
    s=str(switch).strip()
    if s.lower() not in NULL_TOKENS:
        if s=="BtoA": return "B_to_A"
        if s=="AtoB": return "A_to_B"
    s=str(shift).strip()
    if s.lower() not in NULL_TOKENS:
        if "TowardA" in s: return "B_to_A"
        if "TowardB" in s: return "A_to_B"
    return None


def build_master(comp, norm, migration_axis, top_gainer_pct):
    m=comp.merge(norm, on=["chr","bed_start","bed_end"], how="inner")
    m["dUNC_b10"]=np.log2(m["b10_unc"]/m["b10_none"])
    m["dUNC_parental"]=np.log2(m["parental_unc"]/m["parental_none"])
    m["migration_dH3K27me3"]=np.log2(m["b10_none"]/m["parental_none"])
    m["baseline_b10"]=np.log2(m["b10_none"]); m["baseline_parental"]=np.log2(m["parental_none"])
    if migration_axis=="dmso":
        both,sw,sh=MIG_BOTH_DMSO,MIG_SWITCH_DMSO,MIG_SHIFT_DMSO; pc1_par,pc1_b10=PC1_PAR_DMSO,PC1_B10_DMSO
    else:
        both,sw,sh=MIG_BOTH_ORIG,MIG_SWITCH_ORIG,MIG_SHIFT_ORIG; pc1_par,pc1_b10=PC1_PAR_ORIG,PC1_B10_ORIG
    m["mig_dir"]=[categorize_direction(r.get(both),r.get(sw),r.get(sh),r.get(pc1_par),r.get(pc1_b10)) for _,r in m.iterrows()]
    m["is_mig_switcher"]=m["mig_dir"].isin(["A_to_B","B_to_A"])
    m["unc_dir_b10"]=[categorize_movement(r.get(UNC_BOTH_B10),r.get(UNC_SWITCH_B10),r.get(UNC_SHIFT_B10)) for _,r in m.iterrows()]
    mig=m["mig_dir"]; unc=m["unc_dir_b10"]
    moved=mig.isin(["A_to_B","B_to_A"]) & unc.isin(["A_to_B","B_to_A"])
    m["is_reverter_b10"]=moved & (mig!=unc)
    def _hard(v):
        s=str(v).strip(); return s if s in ("AtoB","BtoA") else None
    mig_hard=m[MIG_SWITCH_DMSO if migration_axis=="dmso" else MIG_SWITCH_ORIG].map(_hard).map({"BtoA":"B_to_A","AtoB":"A_to_B"})
    unc_hard=m[UNC_SWITCH_B10].map(_hard).map({"BtoA":"B_to_A","AtoB":"A_to_B"})
    m["is_reverter_b10_strict"]=(mig_hard.isin(["A_to_B","B_to_A"]) & unc_hard.isin(["A_to_B","B_to_A"]) & (mig_hard!=unc_hard))
    thr=np.nanpercentile(m["migration_dH3K27me3"], 100-top_gainer_pct)
    m["is_gainer"]=m["migration_dH3K27me3"]>=thr
    m["is_stable"]=m["mig_dir"].isin(["Stable_A","Stable_B"]) & (~m["is_gainer"])
    return m


def cliffs_delta(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    x=x[np.isfinite(x)]; y=y[np.isfinite(y)]
    n1,n2=len(x),len(y)
    if n1==0 or n2==0: return np.nan
    U1=mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2.0*U1/(n1*n2)-1.0


def baseline_residualize(values, baseline, is_ref, n_ref=50, min_per_bin=10):
    """Continuously remove the baseline-dependent TREND from `values`, using the
    reference (stable) bins. dUNC=log2(unc)-baseline carries a structural slope on
    baseline; decile matching removes it in 10 steps, so a target sitting high/low
    WITHIN its deciles leaks the slope as a fake directional effect. Fit the stable
    median of `values` vs baseline (n_ref quantile knots) and subtract from every
    bin: a target tracking the trend then has residual ~0; one that deviates keeps
    a non-zero residual."""
    values=np.asarray(values,float); baseline=np.asarray(baseline,float); is_ref=np.asarray(is_ref,bool)
    rb,rv=baseline[is_ref],values[is_ref]
    keep=np.isfinite(rb)&np.isfinite(rv); rb,rv=rb[keep],rv[keep]
    if len(rb)<2*min_per_bin: return values-(np.nanmedian(rv) if len(rv) else 0.0)
    k=max(5, min(n_ref, len(rb)//min_per_bin))
    edges=np.unique(np.quantile(rb, np.linspace(0,1,k+1)))
    centers,meds=[],[]
    for j in range(len(edges)-1):
        m=(rb>=edges[j])&(rb<=edges[j+1])
        if m.sum()>=min_per_bin: centers.append(np.median(rb[m])); meds.append(np.median(rv[m]))
    if len(centers)<2: return values-np.nanmedian(rv)
    trend=np.interp(baseline, np.array(centers), np.array(meds))
    return values-trend


def stratified_matched_test(values, is_target, baseline, n_strata=10, n_perm=1000, seed=0, n_boot=0):
    values=np.asarray(values,float); is_target=np.asarray(is_target,bool); baseline=np.asarray(baseline,float)
    keep=np.isfinite(values)&np.isfinite(baseline)
    v,t,b=values[keep],is_target[keep],baseline[keep]
    if n_strata<=1:
        strata=np.zeros(len(v),int); n_strata=1
    else:
        edges=np.quantile(b, np.linspace(0,1,n_strata+1)); edges[0]-=1e-9; edges[-1]+=1e-9
        strata=np.clip(np.digitize(b, edges[1:-1]), 0, n_strata-1)
    def pooled_delta(target_labels):
        num=den=0.0
        for s in range(n_strata):
            idx=strata==s
            xs=v[idx & target_labels]; ys=v[idx & ~target_labels]
            if len(xs)==0 or len(ys)==0: continue
            d=cliffs_delta(xs,ys)
            if not np.isfinite(d): continue
            w=len(xs)*len(ys); num+=d*w; den+=w
        return num/den if den>0 else np.nan
    obs=pooled_delta(t)
    rng=np.random.default_rng(seed)
    perm=np.empty(n_perm)
    for i in range(n_perm):
        lab=np.zeros_like(t)
        for s in range(n_strata):
            idx=np.where(strata==s)[0]; k=int(t[idx].sum())
            if k==0 or k==len(idx): lab[idx]=t[idx]; continue
            lab[rng.choice(idx, size=k, replace=False)]=True
        perm[i]=pooled_delta(lab)
    perm=perm[np.isfinite(perm)]
    p_two=(np.sum(np.abs(perm)>=abs(obs))+1)/(len(perm)+1) if len(perm) else np.nan
    ci_low=ci_high=np.nan
    if n_boot>0:
        tgt_idx=[np.where((strata==s)&t)[0] for s in range(n_strata)]
        stb_idx=[np.where((strata==s)&~t)[0] for s in range(n_strata)]
        boots=np.empty(n_boot)
        for j in range(n_boot):
            num=den=0.0
            for s in range(n_strata):
                ti,si=tgt_idx[s],stb_idx[s]
                if len(ti)==0 or len(si)==0: continue
                xs=v[rng.choice(ti,size=len(ti),replace=True)]; ys=v[rng.choice(si,size=len(si),replace=True)]
                d=cliffs_delta(xs,ys)
                if not np.isfinite(d): continue
                w=len(ti)*len(si); num+=d*w; den+=w
            boots[j]=num/den if den>0 else np.nan
        boots=boots[np.isfinite(boots)]
        if len(boots): ci_low,ci_high=np.percentile(boots,[2.5,97.5])
    return {"cliffs_delta":obs,"perm_p_two_sided":p_two,"ci_low":ci_low,"ci_high":ci_high,
            "n_target":int(t.sum()),"n_stable":int((~t).sum()),"n_perm_valid":int(len(perm)),
            "median_target":float(np.nanmedian(v[t])) if t.any() else np.nan,
            "median_stable":float(np.nanmedian(v[~t])) if (~t).any() else np.nan}


def run_contrast(master, target_mask, stable_mask, value_col, baseline_col, n_strata, n_perm, seed,
                 n_boot=0, baseline_match="decile", n_ref_bins=50):
    target_mask=np.asarray(target_mask,bool); stable_mask=np.asarray(stable_mask,bool)
    sel=target_mask|stable_mask; sub=master.loc[sel]
    vals=sub[value_col].values.astype(float); base=sub[baseline_col].values.astype(float); is_t=target_mask[sel]
    if baseline_match=="residual":
        vals=baseline_residualize(vals, base, ~is_t, n_ref=n_ref_bins); strata=1
    elif baseline_match=="fine":
        strata=max(n_strata,50)
    else:
        strata=n_strata
    return stratified_matched_test(vals, is_t, base, n_strata=strata, n_perm=n_perm, seed=seed, n_boot=n_boot)


def interpret(delta, p):
    if not np.isfinite(delta) or not np.isfinite(p): return "undetermined"
    if p>=0.05: return "no enrichment"
    return "DEPLETION enriched at target (hypothesis-consistent)" if delta<0 else "depletion DE-enriched at target (opposite)"


# --- colour palettes -------------------------------------------------------------
# B10 palette (default) and an alternate Parental palette, so that when both arms are
# drawn on the SAME axis the same category never appears in the same colour twice.
# 'stable' stays red in both arms on purpose: it is the shared reference.
CAT_COLORS={"stable":"#b2182b","switcher":"#0e12e0","reverter":"#011b1d","gainer":"#0be96b"}
CAT_COLORS_PARENTAL={"stable":"#b2182b","switcher":"#59a14f","reverter":"#8c564b","gainer":"#f28e2b"}


def cat_color(cat, style=None):
    """Category -> fill colour. style='Parental' selects the alternate palette. Colour
    is resolved ONCE, at panel construction, so plotting functions never have to guess
    a panel's category back from its colour (that reverse-lookup is fragile)."""
    if style=="Parental": return CAT_COLORS_PARENTAL.get(cat, CAT_COLORS.get(cat,"#777777"))
    return CAT_COLORS.get(cat,"#777777")


def cliff_band(delta):
    a=abs(delta)
    if not np.isfinite(delta): return "na"
    if a<0.147: return "negligible"
    if a<0.330: return "small"
    if a<0.474: return "medium"
    return "large"


def cliff_tag(delta, p=np.nan):
    """Plain-language reading of the delta, so a reader does not have to decode the sign
    convention from memory. delta<0 = target loses MORE than its same-compartment stable
    ('extra loss'); delta>0 = loses LESS ('spared'). Non-significant or negligible
    effects get no direction word, because their sign is not interpretable."""
    if not np.isfinite(delta): return "na"
    band=cliff_band(delta)
    if np.isfinite(p) and p>=0.05: return f"{band} · n.s."
    if band=="negligible": return band
    return f"{band} · {'spared' if delta>0 else 'extra loss'}"


def p_stars(p):
    if not np.isfinite(p): return ""
    if p<0.001: return "***"
    if p<0.01: return "**"
    if p<0.05: return "*"
    return "ns"


def _panel_xlabel(p, n_fallback):
    """x tick label: name / [group] / (n=...) each on its OWN line."""
    parts=[p["label"]]
    g=p.get("group")
    if g: parts.append(f"[{g}]")
    parts.append(f"(n={p.get('n', n_fallback)})")
    return "\n".join(parts)


def plot_ecdf(panels, value_desc, outpath):
    try:
        fig,ax=plt.subplots(figsize=(6.5,5), constrained_layout=True)
        for p in panels:
            x=np.sort(np.asarray(p["values"],float)); x=x[np.isfinite(x)]
            if len(x)==0: continue
            lw=2.6 if p.get("ref") else 1.8
            ax.plot(x, np.linspace(0,1,len(x)), lw=lw, color=p["color"], label=f"{p['label']} (n={len(x)})")
        ax.axvline(0, color="grey", ls=":", lw=1)
        ax.set_xlabel(value_desc, fontsize=12); ax.set_ylabel("ECDF", fontsize=12)
        ax.tick_params(axis='both', labelsize=12)
        ax.set_title("UNC H3K27me3 effect: reprogrammed sets vs stable\n(left shift = more depletion)")
        ax.legend(fontsize=12, frameon=False, loc="lower right")
        fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] ECDF plot failed: {e}")


def plot_violin(panels, value_desc, outpath, show_box=True, cell_tag=""):
    try:
        usable=[p for p in panels if len(np.asarray(p["values"],float)[np.isfinite(np.asarray(p["values"],float))])>=2]
        if len(usable)<2: print("[warn] violin: <2 usable panels, skipping"); return
        clean=[np.asarray(p["values"],float) for p in usable]; clean=[v[np.isfinite(v)] for v in clean]
        allv=np.concatenate(clean); lo,hi=np.nanpercentile(allv,[1,99])
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi<=lo: lo,hi=float(np.nanmin(allv)),float(np.nanmax(allv))
        rng=max(hi-lo,1e-9)
        trimmed=[]
        for v in clean:
            vt=v[(v>=lo)&(v<=hi)]; trimmed.append(vt if len(vt)>=2 else v)
        pos=np.arange(1,len(usable)+1)
        fig,ax=plt.subplots(figsize=(1.4*len(usable)+1.5,5.5), constrained_layout=True)
        vp=ax.violinplot(trimmed, positions=pos, showmedians=False, showextrema=False, widths=0.5)
        for body,p in zip(vp["bodies"], usable):
            body.set_facecolor(p["color"]); body.set_edgecolor("black"); body.set_alpha(0.85); body.set_linewidth(0.6)
        if show_box:
            for x,v in zip(pos,clean):
                if len(v)<5: continue
                q1,med,q3=np.percentile(v,[25,50,75]); wlo,whi=np.percentile(v,[5,95])
                ax.plot([x,x],[max(wlo,lo),min(whi,hi)], color="black", lw=1.0, zorder=3)
                ax.plot([x,x],[q1,q3], color="black", lw=5.0, solid_capstyle="butt", zorder=4)
                ax.scatter([x],[med], color="white", s=20, zorder=5, edgecolor="black", linewidth=0.6)
        ax.axhline(0, color="grey", ls=":", lw=1, zorder=1)
        ref=next((p for p in usable if p.get("ref")), None)
        if ref is not None:
            rv=np.asarray(ref["values"],float); rv=rv[np.isfinite(rv)]
            if len(rv): ax.axhline(np.median(rv), color=CAT_COLORS["stable"], ls="--", lw=1, alpha=0.7, zorder=1, label="stable median")
        for x,p in zip(pos,usable):
            if p.get("ref"): continue
            d=p.get("delta",np.nan); pv=p.get("p",np.nan)
            ax.text(x, hi+0.05*rng, p_stars(pv), ha="center", va="bottom", fontsize=12, fontweight="bold")
            ax.text(x, hi+0.18*rng, f"δ={d:+.2f}\n({cliff_tag(d,pv)})", ha="center", va="bottom", fontsize=10)
        ax.set_xticks(pos)
        ax.set_xticklabels([_panel_xlabel(p,len(c)) for p,c in zip(usable,clean)], fontsize=11)
        ax.set_ylim(lo-0.12*rng, hi+0.42*rng); ax.set_ylabel(value_desc)
        ttl = "UNC H3K27me3 effect by reprogrammed set"
        if cell_tag: ttl = f"{ttl}   [{cell_tag}]"
        ax.set_title(ttl + "\n(δ<0 = loses more than its stable reference; δ>0 = spared; * = perm p)")
        fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] violin plot failed: {e}")


def plot_violin_grouped(panels, value_desc, outpath, group_order, title="", group_colors=None,
                        subtitle="(within-group δ vs that group's own stable reference)"):
    """Combined violin: several GROUPS on ONE axis, each group carrying its OWN stable
    reference. Every panel dict must have panel['group'] (a label in group_order); the
    group's reference panel is flagged panel['ref']=True. Draws a dashed median line per
    group spanning only that group's x-range, a divider between groups, and annotates
    each non-ref panel with its δ vs ITS OWN group's stable.

    Groups can be cell lines ("B10"/"Parental") or compartments ("A-compartment"/...).

    INTERPRETATION: when groups are different ARMS, cross-group violin HEIGHTS are not
    comparable (each arm is its own normalization frame). The readable quantity is the
    WITHIN-group target-minus-stable gap, which is what the δ annotations report.
    """
    try:
        default_gc={"B10":"#1f77b4","Parental":"#ff7f0e"}
        gc=dict(default_gc); gc.update(group_colors or {})
        fallback=["#1f77b4","#ff7f0e","#2ca02c","#9467bd"]
        ordered=[]
        for g in group_order:
            ordered+= [p for p in panels if p.get("group")==g]
        usable=[p for p in ordered if len(np.asarray(p["values"],float)[np.isfinite(np.asarray(p["values"],float))])>=2]
        if len(usable)<2: print("[warn] grouped violin: <2 usable panels, skipping"); return
        clean=[np.asarray(p["values"],float)[np.isfinite(np.asarray(p["values"],float))] for p in usable]
        allv=np.concatenate(clean); lo,hi=np.nanpercentile(allv,[1,99])
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi<=lo: lo,hi=float(np.nanmin(allv)),float(np.nanmax(allv))
        rng=max(hi-lo,1e-9)
        trimmed=[]
        for v in clean:
            vt=v[(v>=lo)&(v<=hi)]; trimmed.append(vt if len(vt)>=2 else v)
        pos=np.arange(1,len(usable)+1)
        fig,ax=plt.subplots(figsize=(1.35*len(usable)+2.0,5.8), constrained_layout=True)
        vp=ax.violinplot(trimmed, positions=pos, showmedians=False, showextrema=False, widths=0.5)
        for body,p in zip(vp["bodies"], usable):
            body.set_facecolor(p["color"]); body.set_edgecolor("black"); body.set_alpha(0.85); body.set_linewidth(0.6)
        for x,v in zip(pos,clean):
            if len(v)<5: continue
            q1,med,q3=np.percentile(v,[25,50,75]); wlo,whi=np.percentile(v,[5,95])
            ax.plot([x,x],[max(wlo,lo),min(whi,hi)], color="black", lw=1.0, zorder=3)
            ax.plot([x,x],[q1,q3], color="black", lw=5.0, solid_capstyle="butt", zorder=4)
            ax.scatter([x],[med], color="white", s=20, zorder=5, edgecolor="black", linewidth=0.6)
        ax.axhline(0, color="grey", ls=":", lw=1, zorder=1)
        groups_present=[g for g in group_order if any(p.get("group")==g for p in usable)]
        for i,g in enumerate(groups_present):
            gc.setdefault(g, fallback[i%len(fallback)])
        xpos_by_group={g:[x for x,p in zip(pos,usable) if p.get("group")==g] for g in groups_present}
        for g in groups_present:
            gx=xpos_by_group[g]
            refp=next((p for p in usable if p.get("group")==g and p.get("ref")), None)
            if refp is not None:
                rv=np.asarray(refp["values"],float); rv=rv[np.isfinite(rv)]
                if len(rv):
                    ax.plot([min(gx)-0.45,max(gx)+0.45],[np.median(rv)]*2, color=gc[g],
                            ls="--", lw=1.3, alpha=0.8, zorder=2)
        for ga,gb in zip(groups_present[:-1],groups_present[1:]):
            xr=(max(xpos_by_group[ga])+min(xpos_by_group[gb]))/2.0
            ax.axvline(xr, color="grey", lw=0.8, alpha=0.5, ls="-")
        for x,p in zip(pos,usable):
            if p.get("ref"): continue
            d=p.get("delta",np.nan); pv=p.get("p",np.nan)
            ax.text(x, hi+0.05*rng, p_stars(pv), ha="center", va="bottom", fontsize=11, fontweight="bold")
            ax.text(x, hi+0.18*rng, f"δ={d:+.2f}\n({cliff_tag(d,pv)})", ha="center", va="bottom", fontsize=9)
        ax.set_xticks(pos)
        ax.set_xticklabels([_panel_xlabel(p,len(c)) for p,c in zip(usable,clean)], fontsize=10)
        ax.set_ylim(lo-0.12*rng, hi+0.44*rng); ax.set_ylabel(value_desc)
        from matplotlib.lines import Line2D
        handles=[Line2D([0],[0], color=gc[g], ls="--", lw=1.3, label=f"{g} stable median") for g in groups_present]
        ax.legend(handles=handles, fontsize=9, frameon=False, loc="lower right")
        base_ttl=title or "Compartment-matched UNC H3K27me3 selectivity"
        ax.set_title(base_ttl + "\n" + subtitle)
        fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] grouped violin failed: {e}")


def plot_balance(master, target_mask, baseline_col, name, outpath):
    try:
        fig,ax=plt.subplots(figsize=(6,4), constrained_layout=True)
        ax.hist(master[master["is_stable"]][baseline_col].dropna(), bins=40, density=True, alpha=0.5, label="stable")
        ax.hist(master[target_mask][baseline_col].dropna(), bins=40, density=True, alpha=0.5, label=name)
        ax.set_xlabel("baseline abundance  log2(none signal)"); ax.set_ylabel("density")
        ax.set_title(f"Baseline-abundance balance: {name} vs stable")
        ax.legend(fontsize=12, frameon=False); fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] balance plot failed: {e}")


def plot_forest(results, outpath, title=None):
    try:
        df=results.copy()
        df["_c"]=df["cell_line"].map({"B10":0,"Parental":1}).fillna(9); df["_o"]=range(len(df))
        df=(df.sort_values(["_c","_o"], kind="stable").drop(columns=["_c","_o"]).reset_index(drop=True))
        df=df.iloc[::-1].reset_index(drop=True)
        colors={"B10":"#1f77b4","Parental":"#ff7f0e"}; y=np.arange(len(df))
        fig,ax=plt.subplots(figsize=(11, 0.52*len(df)+1.6), constrained_layout=True)
        for i,r in df.iterrows():
            c=colors.get(r["cell_line"],"k")
            lo,hi=r.get("ci_low",np.nan),r.get("ci_high",np.nan)
            if np.isfinite(lo) and np.isfinite(hi): ax.plot([lo,hi],[i,i], color=c, lw=2.4, alpha=0.85, zorder=1)
            ax.scatter(r["cliffs_delta"], i, color=c, s=50, zorder=2)
            star=" *" if (np.isfinite(r["perm_p_two_sided"]) and r["perm_p_two_sided"]<0.05) else ""
            ax.text(r["cliffs_delta"], i+0.24, f"δ={r['cliffs_delta']:+.3f}{star} (n={r['n_target']})", fontsize=9, ha="center", va="bottom", color=c)
        ax.axvline(0, color="grey", ls=":", lw=1); ax.set_yticks(y)
        ax.set_yticklabels([f"{r['cell_line']}: {r['target_set']}" for _,r in df.iterrows()], fontsize=9)
        ax.set_ylim(-0.7, len(df)-1+0.95)
        xmin=float(np.nanmin([results["ci_low"].min(), results["cliffs_delta"].min()]))
        xmax=float(np.nanmax([results["ci_high"].max(), results["cliffs_delta"].max()]))
        if not (np.isfinite(xmin) and np.isfinite(xmax)) or xmax<=xmin: xmin,xmax=-1.0,1.0
        pad=0.10*(xmax-xmin); ax.set_xlim(xmin-pad, xmax+pad)
        ax.set_xlabel("Cliff's δ   (negative = H3K27me3 depletion enriched at target)", fontsize=11)
        ax.set_title(title or "UNC selective-erasure effect sizes\n(point = δ, bar = bootstrap 95% CI, * = perm p<0.05)")
        from matplotlib.lines import Line2D
        handles=[Line2D([0],[0], marker="o", color="w", label=k, markerfacecolor=v, markersize=9) for k,v in colors.items()]
        ax.legend(handles=handles, fontsize=9, frameon=True, framealpha=0.9, loc="lower left")
        fig.savefig(outpath, dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] forest plot failed: {e}")


# =============================================================================
# COMPARTMENT-MATCHED REFERENCE ANALYSIS
# Each target set is contrasted ONLY against stable bins of the SAME baseline
# compartment IN THE ARM BEING MEASURED, cancelling the dUNC~PC1 gradient that D5
# exposes. A switcher flips compartment between cell lines, so matching flips per arm.
# =============================================================================
def run_compartment_matched(master, masks, args, combined=True, gainers=True, per_arm=True):
    stabA=masks["stabA"]; stabB=masks["stabB"]
    atb=masks["atb"]; bta=masks["bta"]; rev_atb=masks["rev_atb"]; rev_bta=masks["rev_bta"]
    # H3K27me3 gainers = the most DIRECT proxy for "migration-installed mark" (loci that
    # actually gained H3K27me3 parent->B10), unlike switchers which are only compartment
    # proxies. They are NOT compartment-defined, so we split them by their OWN baseline
    # compartment IN EACH ARM (gainA_b10 = gainer & A-in-B10, etc.) and match to same-
    # compartment stable. Migration-specific erasure would read: negative delta in B10
    # (gainers depleted vs stable) but ~0 in Parental (no installed mark to erase yet).
    has_gainers = gainers and all(masks.get(k) is not None for k in ("gainA_b10","gainB_b10","gainA_par","gainB_par"))
    gainA_b10=masks.get("gainA_b10"); gainB_b10=masks.get("gainB_b10")
    gainA_par=masks.get("gainA_par"); gainB_par=masks.get("gainB_par")

    # (cell, ref_label, ref_mask, target_set, target_mask, value_col, baseline_col)
    #   B10 arm : A_to_B is B-in-B10 -> Stable_B ; B_to_A is A-in-B10 -> Stable_A
    #             rev_atb starts A (UNC A->B) -> Stable_A ; rev_bta starts B -> Stable_B
    #   Parental: A_to_B is A-in-Parental -> Stable_A ; B_to_A is B-in-Parental -> Stable_B (no reverters)
    plan=[
        ("B10","Stable_A",stabA,"switchers_B_to_A", bta,    "dUNC_b10","baseline_b10"),
        ("B10","Stable_A",stabA,"reverters_A_to_B", rev_atb,"dUNC_b10","baseline_b10"),
        ("B10","Stable_B",stabB,"switchers_A_to_B", atb,    "dUNC_b10","baseline_b10"),
        ("B10","Stable_B",stabB,"reverters_B_to_A", rev_bta,"dUNC_b10","baseline_b10"),
        ("Parental","Stable_A",stabA,"switchers_A_to_B", atb,"dUNC_parental","baseline_parental"),
        ("Parental","Stable_B",stabB,"switchers_B_to_A", bta,"dUNC_parental","baseline_parental"),
    ]
    if has_gainers:
        plan += [
            ("B10","Stable_A",stabA,"gainers_A", gainA_b10,"dUNC_b10","baseline_b10"),
            ("B10","Stable_B",stabB,"gainers_B", gainB_b10,"dUNC_b10","baseline_b10"),
            ("Parental","Stable_A",stabA,"gainers_A", gainA_par,"dUNC_parental","baseline_parental"),
            ("Parental","Stable_B",stabB,"gainers_B", gainB_par,"dUNC_parental","baseline_parental"),
        ]

    rows=[]
    for i,(cell,ref_label,ref_mask,tname,tmask,vcol,bcol) in enumerate(plan):
        res=run_contrast(master, tmask, ref_mask, vcol, bcol, args.n_strata, args.n_perm,
                         seed=args.seed+100+i, n_boot=args.n_boot,
                         baseline_match=args.baseline_match, n_ref_bins=args.n_ref_bins)
        res.update(cell_line=cell, ref_compartment=ref_label, target_set=tname, value_col=vcol)
        res["interpretation"]=interpret(res["cliffs_delta"], res["perm_p_two_sided"])
        res["reading"]=cliff_tag(res["cliffs_delta"], res["perm_p_two_sided"])
        rows.append(res)
        print(f"  [cmatch {cell:8s}] {ref_label} vs {tname:18s} delta={res['cliffs_delta']:+.3f}  "
              f"p={res['perm_p_two_sided']:.4f}  (nT={res['n_target']}, nRef={res['n_stable']})  -> {res['interpretation']}")

    cm=pd.DataFrame(rows)[["cell_line","ref_compartment","target_set","value_col","n_target","n_stable",
                           "cliffs_delta","ci_low","ci_high","perm_p_two_sided","median_target",
                           "median_stable","n_perm_valid","reading","interpretation"]]
    cm_path=os.path.join(args.outdir,"compartment_matched_results.tsv")
    cm.to_csv(cm_path, sep="\t", index=False)
    print(f"[info] compartment-matched results -> {cm_path}")

    # forest (reuse plot_forest; fold the reference compartment into the row label)
    cm_forest=cm.copy()
    cm_forest["target_set"]=cm_forest["ref_compartment"]+" vs "+cm_forest["target_set"]
    plot_forest(cm_forest, os.path.join(args.outdir,"forest_effect_sizes_cmatched.png"),
                title="Compartment-matched UNC selectivity (each target vs same-compartment stable)\n"
                      "(point = δ, bar = bootstrap 95% CI, * = perm p<0.05)")

    cm_lookup={(r["cell_line"],r["target_set"]):r for _,r in cm.iterrows()}

    def _cpanel(label, mask, cat, set_name, cell, value_col, ref=False, group=None, style=None):
        if mask is None: return None
        mask=np.asarray(mask,bool); n=int(mask.sum())
        if not ref and n<args.min_plot_n: return None
        row=cm_lookup.get((cell,set_name)) if set_name else None
        return {"label":label,"color":cat_color(cat,style),"ref":ref,"n":n,"group":group,
                "values":master.loc[mask,value_col].dropna().values,
                "delta":(row["cliffs_delta"] if row is not None else np.nan),
                "p":(row["perm_p_two_sided"] if row is not None else np.nan)}

    def _add(lst, panel):
        if panel is not None: lst.append(panel)
        return lst

    # ---- 4 single-arm violins (one arm, one compartment; unambiguous -> supplement) ----
    p1=[p for p in [_cpanel("Stable A",stabA,"stable",None,"B10","dUNC_b10",ref=True),
                    _cpanel("B→A switch",bta,"switcher","switchers_B_to_A","B10","dUNC_b10"),
                    _cpanel("A→B revert",rev_atb,"reverter","reverters_A_to_B","B10","dUNC_b10")] if p is not None]
    if has_gainers: _add(p1,_cpanel("gainer",gainA_b10,"gainer","gainers_A","B10","dUNC_b10"))
    plot_violin(p1, "dUNC = log2(UNC/none)  [B10]", os.path.join(args.outdir,"violin_cmatched_b10_stableA.png"),
                cell_tag="B10 · A-compartment vs Stable A")

    p2=[p for p in [_cpanel("Stable B",stabB,"stable",None,"B10","dUNC_b10",ref=True),
                    _cpanel("A→B switch",atb,"switcher","switchers_A_to_B","B10","dUNC_b10"),
                    _cpanel("B→A revert",rev_bta,"reverter","reverters_B_to_A","B10","dUNC_b10")] if p is not None]
    if has_gainers: _add(p2,_cpanel("gainer",gainB_b10,"gainer","gainers_B","B10","dUNC_b10"))
    plot_violin(p2, "dUNC = log2(UNC/none)  [B10]", os.path.join(args.outdir,"violin_cmatched_b10_stableB.png"),
                cell_tag="B10 · B-compartment vs Stable B")

    p3=[p for p in [_cpanel("Stable A",stabA,"stable",None,"Parental","dUNC_parental",ref=True,style="Parental"),
                    _cpanel("A→B switch",atb,"switcher","switchers_A_to_B","Parental","dUNC_parental",style="Parental")] if p is not None]
    if has_gainers: _add(p3,_cpanel("gainer",gainA_par,"gainer","gainers_A","Parental","dUNC_parental",style="Parental"))
    plot_violin(p3, "dUNC = log2(UNC/none)  [Parental]", os.path.join(args.outdir,"violin_cmatched_parental_stableA.png"),
                cell_tag="Parental · A-compartment vs Stable A")

    p4=[p for p in [_cpanel("Stable B",stabB,"stable",None,"Parental","dUNC_parental",ref=True,style="Parental"),
                    _cpanel("B→A switch",bta,"switcher","switchers_B_to_A","Parental","dUNC_parental",style="Parental")] if p is not None]
    if has_gainers: _add(p4,_cpanel("gainer",gainB_par,"gainer","gainers_B","Parental","dUNC_parental",style="Parental"))
    plot_violin(p4, "dUNC = log2(UNC/none)  [Parental]", os.path.join(args.outdir,"violin_cmatched_parental_stableB.png"),
                cell_tag="Parental · B-compartment vs Stable B")

    # ---- 2 COMBINED panels: one compartment, BOTH arms (each arm keeps its own stable) ----
    if combined:
        cA=[p for p in [
            _cpanel("Stable A",stabA,"stable",None,"B10","dUNC_b10",ref=True,group="B10"),
            _cpanel("B→A switch",bta,"switcher","switchers_B_to_A","B10","dUNC_b10",group="B10"),
            _cpanel("A→B revert",rev_atb,"reverter","reverters_A_to_B","B10","dUNC_b10",group="B10"),
        ] if p is not None]
        if has_gainers: _add(cA,_cpanel("gainer",gainA_b10,"gainer","gainers_A","B10","dUNC_b10",group="B10"))
        _add(cA,_cpanel("Stable A",stabA,"stable",None,"Parental","dUNC_parental",ref=True,group="Parental",style="Parental"))
        _add(cA,_cpanel("A→B switch",atb,"switcher","switchers_A_to_B","Parental","dUNC_parental",group="Parental",style="Parental"))
        if has_gainers: _add(cA,_cpanel("gainer",gainA_par,"gainer","gainers_A","Parental","dUNC_parental",group="Parental",style="Parental"))
        plot_violin_grouped(cA, "dUNC = log2(UNC/none)", os.path.join(args.outdir,"violin_cmatched_combined_stableA.png"),
                            group_order=["B10","Parental"], title="Compartment-matched selectivity · A-compartment sets",
                            subtitle="(within-arm δ vs same-compartment stable; heights not comparable across arms)")

        cB=[p for p in [
            _cpanel("Stable B",stabB,"stable",None,"B10","dUNC_b10",ref=True,group="B10"),
            _cpanel("A→B switch",atb,"switcher","switchers_A_to_B","B10","dUNC_b10",group="B10"),
            _cpanel("B→A revert",rev_bta,"reverter","reverters_B_to_A","B10","dUNC_b10",group="B10"),
        ] if p is not None]
        if has_gainers: _add(cB,_cpanel("gainer",gainB_b10,"gainer","gainers_B","B10","dUNC_b10",group="B10"))
        _add(cB,_cpanel("Stable B",stabB,"stable",None,"Parental","dUNC_parental",ref=True,group="Parental",style="Parental"))
        _add(cB,_cpanel("B→A switch",bta,"switcher","switchers_B_to_A","Parental","dUNC_parental",group="Parental",style="Parental"))
        if has_gainers: _add(cB,_cpanel("gainer",gainB_par,"gainer","gainers_B","Parental","dUNC_parental",group="Parental",style="Parental"))
        plot_violin_grouped(cB, "dUNC = log2(UNC/none)", os.path.join(args.outdir,"violin_cmatched_combined_stableB.png"),
                            group_order=["B10","Parental"], title="Compartment-matched selectivity · B-compartment sets",
                            subtitle="(within-arm δ vs same-compartment stable; heights not comparable across arms)")

    # ---- 2 PER-ARM panels: one arm, BOTH compartment groups side by side ----
    # Parental : [Stable A | A→B switch | gainerA]  ||  [Stable B | B→A switch | gainerB]
    # B10      : [Stable A | B→A switch | A→B revert | gainerA] || [Stable B | A→B switch | B→A revert | gainerB]
    # Both groups here are the SAME cell line, so heights ARE comparable within the panel.
    if per_arm:
        GRP_A="A-compartment"; GRP_B="B-compartment"
        comp_colors={GRP_A:"#e6550d", GRP_B:"#3182bd"}

        par=[]
        _add(par,_cpanel("Stable A",stabA,"stable",None,"Parental","dUNC_parental",ref=True,group=GRP_A,style="Parental"))
        _add(par,_cpanel("A→B switch",atb,"switcher","switchers_A_to_B","Parental","dUNC_parental",group=GRP_A,style="Parental"))
        if has_gainers: _add(par,_cpanel("gainer",gainA_par,"gainer","gainers_A","Parental","dUNC_parental",group=GRP_A,style="Parental"))
        _add(par,_cpanel("Stable B",stabB,"stable",None,"Parental","dUNC_parental",ref=True,group=GRP_B,style="Parental"))
        _add(par,_cpanel("B→A switch",bta,"switcher","switchers_B_to_A","Parental","dUNC_parental",group=GRP_B,style="Parental"))
        if has_gainers: _add(par,_cpanel("gainer",gainB_par,"gainer","gainers_B","Parental","dUNC_parental",group=GRP_B,style="Parental"))
        plot_violin_grouped(par, "dUNC = log2(UNC/none)  [Parental]",
                            os.path.join(args.outdir,"violin_cmatched_arm_parental.png"),
                            group_order=[GRP_A,GRP_B], group_colors=comp_colors,
                            title="Compartment-matched selectivity · PARENTAL (both compartments)",
                            subtitle="(δ vs that compartment's own stable; same arm throughout, so heights are comparable)")

        b10=[]
        _add(b10,_cpanel("Stable A",stabA,"stable",None,"B10","dUNC_b10",ref=True,group=GRP_A))
        _add(b10,_cpanel("B→A switch",bta,"switcher","switchers_B_to_A","B10","dUNC_b10",group=GRP_A))
        _add(b10,_cpanel("A→B revert",rev_atb,"reverter","reverters_A_to_B","B10","dUNC_b10",group=GRP_A))
        if has_gainers: _add(b10,_cpanel("gainer",gainA_b10,"gainer","gainers_A","B10","dUNC_b10",group=GRP_A))
        _add(b10,_cpanel("Stable B",stabB,"stable",None,"B10","dUNC_b10",ref=True,group=GRP_B))
        _add(b10,_cpanel("A→B switch",atb,"switcher","switchers_A_to_B","B10","dUNC_b10",group=GRP_B))
        _add(b10,_cpanel("B→A revert",rev_bta,"reverter","reverters_B_to_A","B10","dUNC_b10",group=GRP_B))
        if has_gainers: _add(b10,_cpanel("gainer",gainB_b10,"gainer","gainers_B","B10","dUNC_b10",group=GRP_B))
        plot_violin_grouped(b10, "dUNC = log2(UNC/none)  [B10]",
                            os.path.join(args.outdir,"violin_cmatched_arm_b10.png"),
                            group_order=[GRP_A,GRP_B], group_colors=comp_colors,
                            title="Compartment-matched selectivity · BOTTOM-10 (both compartments)",
                            subtitle="(δ vs that compartment's own stable; same arm throughout, so heights are comparable)")

    return cm


def map_genes(master, gene_table_path, outdir, genes_per_bin="one"):
    gt=pd.read_csv(gene_table_path, sep="\t")
    if not {"gene_name","coordinate"}.issubset(gt.columns):
        sys.exit(f"[ERROR] gene table needs columns gene_name, coordinate; has {list(gt.columns)}")
    gt=gt[gt["coordinate"].isin(set(master["coordinate"]))].copy()
    if not len(gt): sys.exit("[ERROR] no genes fall in the tested bins (coordinate mismatch?)")
    if genes_per_bin=="one" and {"TSS","bin_start","bin_end"}.issubset(gt.columns):
        centre=(gt["bin_start"].astype(float)+gt["bin_end"].astype(float))/2.0
        gt["_dist"]=(gt["TSS"].astype(float)-centre).abs()
        gt=(gt.sort_values("_dist").drop_duplicates("coordinate", keep="first").drop(columns="_dist"))
    elif genes_per_bin=="one":
        gt=gt.drop_duplicates("coordinate", keep="first")
    gldir=os.path.join(outdir,"genelists"); os.makedirs(gldir, exist_ok=True)
    def genes_in(mask):
        coords=set(master.loc[mask,"coordinate"])
        return sorted(set(gt.loc[gt["coordinate"].isin(coords),"gene_name"]))
    sets={"migration_switchers":master["is_mig_switcher"].values,
          "switchers_A_to_B":(master["mig_dir"]=="A_to_B").values,
          "switchers_B_to_A":(master["mig_dir"]=="B_to_A").values,
          "UNC_reverters":master["is_reverter_b10"].values,
          "UNC_reverters_A_to_B":(master["is_reverter_b10"]&(master["unc_dir_b10"]=="A_to_B")).values,
          "UNC_reverters_B_to_A":(master["is_reverter_b10"]&(master["unc_dir_b10"]=="B_to_A")).values,
          "UNC_reverters_strict":master["is_reverter_b10_strict"].values,
          "H3K27me3_gainers":master["is_gainer"].values}
    with open(os.path.join(gldir,"background_universe.txt"),"w") as fh:
        fh.write("\n".join(sorted(set(gt["gene_name"])))+"\n")
    for name,mask in sets.items():
        gl=genes_in(mask)
        with open(os.path.join(gldir,f"{name}.txt"),"w") as fh: fh.write("\n".join(gl)+"\n")
        print(f"    {name:22s} {int(mask.sum()):>5} bins -> {len(gl):>5} genes")
    gt[["gene_name","coordinate"]].to_csv(os.path.join(gldir,"gene_assignment.tsv"), sep="\t", index=False)
    print(f"[info] gene mapping -> {gldir}")
    return gldir


# =============================================================================
# 11. ANOMALY DIAGNOSTICS (D1-D5) -- localize where a directional delta comes from
# =============================================================================
def _binned_median(x, y, nbins=40, min_n=8):
    x=np.asarray(x,float); y=np.asarray(y,float)
    keep=np.isfinite(x)&np.isfinite(y); x,y=x[keep],y[keep]
    if len(x)<min_n*2: return np.array([]),np.array([])
    edges=np.unique(np.quantile(x, np.linspace(0,1,nbins+1)))
    cx,cy=[],[]
    for j in range(len(edges)-1):
        m=(x>=edges[j])&(x<=edges[j+1])
        if m.sum()>=min_n: cx.append(np.median(x[m])); cy.append(np.median(y[m]))
    return np.array(cx),np.array(cy)


def _pc1_cols(migration_axis="dmso"):
    """PC1 column pair for the requested migration axis. Kept in one place so the
    diagnostics use the SAME coordinate system that defined Stable_A/Stable_B."""
    return (PC1_PAR_DMSO,PC1_B10_DMSO) if migration_axis=="dmso" else (PC1_PAR_ORIG,PC1_B10_ORIG)


def _arms(master, migration_axis="dmso"):
    pcp,pcb=_pc1_cols(migration_axis)
    yield ("Parental","dUNC_parental","baseline_parental",pcp)
    yield ("B10","dUNC_b10","baseline_b10",pcb)


def _dir_masks(master):
    return {"A_to_B":(master["mig_dir"]=="A_to_B").values,
            "B_to_A":(master["mig_dir"]=="B_to_A").values,
            "gainer":master["is_gainer"].values}


def diag_dunc_vs_baseline(master, outdir, migration_axis="dmso"):
    """D1: dUNC vs baseline with the stable trend + set medians. Set median ON the
    stable trend -> delta is baseline-trend (artifact); OFF it (below for depletion)
    -> real effect beyond baseline. Dotted drop = median minus stable-trend value."""
    try:
        stab=master["is_stable"].values; dirs=_dir_masks(master)
        colors={"A_to_B":"#2166ac","B_to_A":"#b2182b","gainer":"#1a9850"}
        fig,axes=plt.subplots(1,2, figsize=(13,5.4), constrained_layout=True)
        for ax,(cell,dcol,bcol,_) in zip(axes,_arms(master,migration_axis)):
            base=master[bcol].values; dunc=master[dcol].values
            keep=np.isfinite(base)&np.isfinite(dunc)
            ax.hexbin(base[keep], dunc[keep], gridsize=45, cmap="Greys", bins="log", mincnt=1, zorder=0)
            cx,cy=_binned_median(base[stab], dunc[stab], nbins=40)
            if len(cx): ax.plot(cx,cy, color="black", lw=2.2, zorder=3, label="stable trend")
            if keep.sum()>50:
                beta=np.polyfit(base[keep], dunc[keep], 1)[0]
                ax.text(0.03,0.04, f"global slope β = {beta:+.2f}", transform=ax.transAxes, fontsize=10,
                        bbox=dict(boxstyle="round", fc="white", alpha=0.8))
            for dname,dmask in dirs.items():
                bx=base[dmask]; by=dunc[dmask]; ok=np.isfinite(bx)&np.isfinite(by); bx,by=bx[ok],by[ok]
                if len(bx)<8: continue
                mb,md=np.median(bx),np.median(by)
                exp=np.interp(mb,cx,cy) if len(cx) else np.nan
                ax.scatter([mb],[md], s=120, color=colors[dname], edgecolor="k", zorder=5, label=f"{dname} median")
                if np.isfinite(exp): ax.plot([mb,mb],[exp,md], color=colors[dname], lw=2, ls=":", zorder=4)
            ax.axhline(0, color="grey", ls=":", lw=1)
            ax.set_xlabel(f"baseline = log2(none)  [{cell}]"); ax.set_ylabel("dUNC = log2(UNC/none)")
            ax.set_title(f"{cell}: dUNC vs baseline\ndotted drop = set median minus stable trend")
            ax.legend(fontsize=8, frameon=True, framealpha=0.85, loc="upper right")
        fig.suptitle("D1  Is the directional effect ON the baseline trend (artifact) or OFF it (real)?", fontsize=12)
        fig.savefig(os.path.join(outdir,"D1_dunc_vs_baseline.png"), dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] D1 failed: {e}")


def diag_within_stratum_balance(master, outdir, n_strata=10, migration_axis="dmso"):
    """D2: within each baseline decile, target-minus-stable median baseline. Bars
    away from 0 = decile matching leaves a residual baseline gap (the leak source)."""
    try:
        stab=master["is_stable"].values; dirs=_dir_masks(master)
        contrasts=[(cell,dname,dmask,bcol) for cell,_,bcol,_ in _arms(master,migration_axis) for dname,dmask in dirs.items()]
        ncol=3; nrow=int(np.ceil(len(contrasts)/ncol))
        fig,axes=plt.subplots(nrow,ncol, figsize=(4.4*ncol,3.2*nrow), constrained_layout=True, squeeze=False)
        for ax,(cell,dname,dmask,bcol) in zip(axes.ravel(),contrasts):
            sel=dmask|stab; b=master[bcol].values[sel]; t=dmask[sel]
            kk=np.isfinite(b); b,t=b[kk],t[kk]
            edges=np.quantile(b, np.linspace(0,1,n_strata+1)); edges[0]-=1e-9; edges[-1]+=1e-9
            strata=np.clip(np.digitize(b, edges[1:-1]), 0, n_strata-1)
            gaps=[]
            for s in range(n_strata):
                m=strata==s; tb=b[m&t]; sb=b[m&~t]
                gaps.append(np.median(tb)-np.median(sb) if len(tb) and len(sb) else np.nan)
            col="#2166ac" if dname=="A_to_B" else ("#b2182b" if dname=="B_to_A" else "#1a9850")
            ax.bar(np.arange(n_strata), gaps, color=col, alpha=0.85); ax.axhline(0, color="k", lw=0.8)
            ax.set_title(f"{cell}: {dname}", fontsize=10)
            ax.set_xlabel("baseline decile"); ax.set_ylabel("Δbaseline (target−stable)")
        for ax in axes.ravel()[len(contrasts):]: ax.axis("off")
        fig.suptitle("D2  Within-stratum baseline imbalance (the leakage source)\nbars away from 0 = decile matching leaves a residual baseline gap", fontsize=12)
        fig.savefig(os.path.join(outdir,"D2_within_stratum_balance.png"), dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] D2 failed: {e}")


def diag_decile_vs_residual(master, outdir, args):
    """D3 (decisive): delta under decile matching vs continuous residualization. If
    |delta| collapses toward 0 under 'residual', the decile result was baseline-trend
    leakage; if it survives, the effect is real beyond baseline."""
    rows=[]; stab=master["is_stable"].values; dirs=_dir_masks(master)
    contrasts=[(cell,dname,dmask,dcol,bcol) for cell,dcol,bcol,_ in _arms(master,args.migration_axis) for dname,dmask in dirs.items()]
    nb=min(args.n_boot,300)
    for i,(cell,dname,dmask,dcol,bcol) in enumerate(contrasts):
        out={"cell_line":cell,"target_set":dname,"n_target":int(dmask.sum())}
        for mode in ("decile","residual"):
            r=run_contrast(master, dmask, stab, dcol, bcol, n_strata=args.n_strata, n_perm=args.n_perm,
                           seed=args.seed+i, n_boot=nb, baseline_match=mode, n_ref_bins=args.n_ref_bins)
            out[f"delta_{mode}"]=r["cliffs_delta"]; out[f"p_{mode}"]=r["perm_p_two_sided"]
            out[f"cilo_{mode}"]=r["ci_low"]; out[f"cihi_{mode}"]=r["ci_high"]
        dd,dr=out["delta_decile"],out["delta_residual"]
        out["frac_remaining"]=(dr/dd) if (np.isfinite(dd) and abs(dd)>1e-6) else np.nan
        rows.append(out)
    res=pd.DataFrame(rows)
    res.to_csv(os.path.join(outdir,"D3_decile_vs_residual.tsv"), sep="\t", index=False)
    try:
        df=res.iloc[::-1].reset_index(drop=True); y=np.arange(len(df))
        fig,ax=plt.subplots(figsize=(10, 0.6*len(df)+1.6), constrained_layout=True)
        for i,r in df.iterrows():
            c="#ff7f0e" if r["cell_line"]=="Parental" else "#1f77b4"
            if np.isfinite(r["cilo_decile"]) and np.isfinite(r["cihi_decile"]):
                ax.plot([r["cilo_decile"],r["cihi_decile"]],[i+0.16,i+0.16], color=c, lw=2, alpha=0.85)
            ax.scatter(r["delta_decile"], i+0.16, color=c, s=55, zorder=3)
            if np.isfinite(r["cilo_residual"]) and np.isfinite(r["cihi_residual"]):
                ax.plot([r["cilo_residual"],r["cihi_residual"]],[i-0.16,i-0.16], color=c, lw=2, alpha=0.5, ls="--")
            ax.scatter(r["delta_residual"], i-0.16, facecolors="none", edgecolors=c, s=70, linewidths=2, zorder=3)
            ax.text(r["delta_decile"], i+0.40, f"{r['delta_decile']:+.3f}{p_stars(r['p_decile'])}", ha="center", va="bottom", fontsize=8, color=c)
            ax.text(r["delta_residual"], i-0.52, f"{r['delta_residual']:+.3f}{p_stars(r['p_residual'])}", ha="center", va="top", fontsize=8, color=c)
        ax.axvline(0, color="grey", ls=":", lw=1); ax.set_yticks(y)
        ax.set_yticklabels([f"{r['cell_line']}: {r['target_set']}" for _,r in df.iterrows()], fontsize=9)
        ax.set_ylim(-0.9, len(df)-1+0.9)
        ax.set_xlabel("Cliff's δ   (filled = decile match, open = continuous residual)")
        ax.set_title("D3  Decisive test: does the directional δ survive CONTINUOUS baseline control?\nopen markers collapsing toward 0 = decile δ was baseline-trend leakage")
        from matplotlib.lines import Line2D
        ax.legend(handles=[Line2D([0],[0], marker="o", color="grey", lw=0, label="decile match"),
                           Line2D([0],[0], marker="o", markerfacecolor="none", markeredgecolor="grey", lw=0, label="continuous residual")],
                  fontsize=8, frameon=True, loc="lower right")
        fig.savefig(os.path.join(outdir,"D3_decile_vs_residual.png"), dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] D3 forest failed: {e}")
    return res


def diag_num_denom(master, outdir):
    """D4: is depletion driven by the UNC numerator (real loss) or the none
    denominator? Per-set median log2(none) vs log2(unc) against y=x."""
    try:
        signal={"Parental":("parental_none","parental_unc"),"B10":("b10_none","b10_unc")}
        dirs=_dir_masks(master); stab=master["is_stable"].values
        fig,axes=plt.subplots(1,2, figsize=(12,5.4), constrained_layout=True)
        colors={"A_to_B":"#2166ac","B_to_A":"#b2182b","gainer":"#1a9850","stable":"#7f7f7f"}
        rows=[]
        for ax,(cell,(ncol,ucol)) in zip(axes,signal.items()):
            ln=np.log2(master[ncol].values); lu=np.log2(master[ucol].values)
            sets=dict(dirs); sets["stable"]=stab
            for sname,smask in sets.items():
                xn=ln[smask]; yu=lu[smask]; kk=np.isfinite(xn)&np.isfinite(yu)
                if kk.sum()<8: continue
                mx,my=np.median(xn[kk]),np.median(yu[kk])
                ax.scatter([mx],[my], s=130, color=colors[sname], edgecolor="k", zorder=4, label=sname)
                rows.append({"cell_line":cell,"set":sname,"median_log2_none":mx,"median_log2_unc":my,"median_dUNC":my-mx})
            lo=np.nanpercentile(np.r_[ln,lu],2); hi=np.nanpercentile(np.r_[ln,lu],98)
            ax.plot([lo,hi],[lo,hi], color="k", ls="--", lw=1, label="no change (y=x)")
            ax.set_xlabel(f"median log2(none)  [{cell}]"); ax.set_ylabel("median log2(UNC)")
            ax.set_title(f"{cell}: below y=x = depletion"); ax.legend(fontsize=8, frameon=True, framealpha=0.85)
        pd.DataFrame(rows).to_csv(os.path.join(outdir,"D4_num_denom.tsv"), sep="\t", index=False)
        fig.suptitle("D4  Numerator (UNC) vs denominator (none) decomposition", fontsize=12)
        fig.savefig(os.path.join(outdir,"D4_num_denom.png"), dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] D4 failed: {e}")


def diag_pc1(master, outdir, migration_axis="dmso"):
    """D5: compartment strength (PC1) per set + dUNC vs PC1. If switchers are
    boundary loci and dUNC tracks PC1, the effect is a position confound."""
    pcp,pcb=_pc1_cols(migration_axis)
    if pcp not in master.columns or pcb not in master.columns:
        print(f"[info] D5 skipped: PC1 columns ({pcp},{pcb}) not in master"); return
    try:
        stab=master["is_stable"].values; dirs=_dir_masks(master)
        fig,axes=plt.subplots(2,2, figsize=(12,9), constrained_layout=True)
        arms=[("Parental","dUNC_parental",pcp),("B10","dUNC_b10",pcb)]
        colors={"A_to_B":"#2166ac","B_to_A":"#b2182b","gainer":"#1a9850","stable":"#7f7f7f"}
        for j,(cell,dcol,pcol) in enumerate(arms):
            pc1=pd.to_numeric(master[pcol], errors="coerce").values; dunc=master[dcol].values
            ax=axes[0,j]; keep=np.isfinite(pc1)&np.isfinite(dunc)
            ax.hexbin(pc1[keep], dunc[keep], gridsize=45, cmap="Greys", bins="log", mincnt=1)
            cx,cy=_binned_median(pc1[stab], dunc[stab], nbins=40)
            if len(cx): ax.plot(cx,cy, color="black", lw=2, label="stable trend")
            ax.axvline(0, color="grey", ls=":", lw=1); ax.axhline(0, color="grey", ls=":", lw=1)
            ax.set_xlabel(f"PC1 (>0 = A) [{cell}]"); ax.set_ylabel("dUNC"); ax.set_title(f"{cell}: dUNC vs PC1")
            ax.legend(fontsize=8, frameon=True)
            ax=axes[1,j]; sets=dict(dirs); sets["stable"]=stab
            for sname,smask in sets.items():
                v=pc1[smask]; v=v[np.isfinite(v)]
                if len(v)<8: continue
                ax.hist(v, bins=40, density=True, histtype="step", lw=2, color=colors[sname], label=sname)
            ax.axvline(0, color="grey", ls=":", lw=1)
            ax.set_xlabel(f"PC1 [{cell}]"); ax.set_ylabel("density"); ax.set_title(f"{cell}: compartment strength by set")
            ax.legend(fontsize=8, frameon=True)
        fig.suptitle("D5  Is the effect a PC1/boundary-position confound?", fontsize=12)
        fig.savefig(os.path.join(outdir,"D5_pc1_position.png"), dpi=150); plt.close(fig)
    except Exception as e: print(f"[warn] D5 failed: {e}")


def run_all_diagnostics(master, args, parent_outdir):
    ddir=os.path.join(parent_outdir,"diagnostics"); os.makedirs(ddir, exist_ok=True)
    print(f"[diag] writing anomaly diagnostics -> {ddir}")
    diag_dunc_vs_baseline(master, ddir, args.migration_axis)
    diag_within_stratum_balance(master, ddir, n_strata=args.n_strata, migration_axis=args.migration_axis)
    d3=diag_decile_vs_residual(master, ddir, args)
    diag_num_denom(master, ddir)
    diag_pc1(master, ddir, args.migration_axis)
    print("[diag] D3 decile vs continuous-residual delta (the decisive test):")
    for _,r in d3.iterrows():
        fr=r["frac_remaining"]
        if np.isfinite(fr) and abs(fr)<0.4: tag="COLLAPSES -> decile δ was baseline-trend leakage"
        elif np.isfinite(fr): tag="survives -> effect persists beyond baseline trend"
        else: tag="n/a"
        frs=f"{fr*100:5.0f}% remains" if np.isfinite(fr) else "  n/a"
        print(f"    {r['cell_line']:8s} {r['target_set']:8s}  decile δ={r['delta_decile']:+.3f}{p_stars(r['p_decile'])}  ->  "
              f"residual δ={r['delta_residual']:+.3f}{p_stars(r['p_residual'])}  ({frs})  {tag}")
    return d3


def git_hash():
    try:
        return subprocess.check_output(["git","rev-parse","HEAD"], stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "NA"


def write_manifest(outdir, args, samples, n_bins, results_path, norm):
    lines=["# MANIFEST  --  unc_selective_erasure.py",
           f"generated   : {_dt.datetime.now().isoformat(timespec='seconds')}",
           f"git_commit  : {git_hash()}", f"normalization : {norm}",
           f"baseline_match : {args.baseline_match}",
           f"python      : {sys.version.split()[0]}", "", "## parameters"]
    for k,v in sorted(vars(args).items()): lines.append(f"  {k}: {v}")
    lines+=["", "## samples (role, size factor [greenlist only], BAM)"]
    for _,r in samples.iterrows():
        sf=r["deseq2_size_factor"]; sf_str=f"{sf:.4f}" if pd.notna(sf) else "NA(cpm)"
        lines.append(f"  {r['SampleID']:28s} {r['role']:14s} sf={sf_str}  {r['bam']}")
    lines+=["", f"## bins analysed : {n_bins}", f"## results       : {results_path}", ""]
    with open(os.path.join(outdir,"MANIFEST.md"),"w") as fh: fh.write("\n".join(lines))


def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    keep_chroms=[f"chr{i}" for i in range(1,23)]+(["chrX"] if args.include_chrx else [])
    norm=args.norm
    if norm=="greenlist" and not args.size_factors:
        print("[warn] --norm greenlist requested but no --size-factors; falling back to cpm"); norm="cpm"
    drop_samples=[s for s in (args.drop_samples or "").split(",")]
    samples=load_samples(args.size_factors, args.sample_sheet, args.project_root, factor=args.factor, drop_samples=drop_samples)
    print(f"[info] normalization: {norm} | baseline-match: {args.baseline_match}")
    comp=load_bins(args.compartments, keep_chroms)
    if args.counts_table:
        raw=load_raw_counts(args.counts_table, samples)
    else:
        bins_bed=write_bins_bed(comp, os.path.join(args.outdir,"bins_250kb.bed"))
        raw_path=os.path.join(args.outdir,"raw_bin_counts.tsv")
        count_reads(bins_bed, samples, out_npz=os.path.join(args.outdir,"bin_counts.npz"), out_raw=raw_path,
                    blacklist=args.blacklist, threads=args.threads, min_mapq=args.min_mapq,
                    extend_reads=args.extend_reads, use_deeptools_blacklist=args.deeptools_blacklist)
        raw=load_raw_counts(raw_path, samples)
    norm_signal=normalize_and_aggregate(raw, samples, args.pseudocount, norm)
    master=build_master(comp, norm_signal, args.migration_axis, args.top_gainer_pct)
    n_lost=len(comp)-len(master)
    if n_lost: print(f"[info] {n_lost} compartment bins had no counts and are absent.")
    if args.blacklist and not args.deeptools_blacklist:
        bl_frac=compute_blacklist_fraction(master, args.blacklist); master["blacklist_frac"]=bl_frac
        print(f"[info] blacklist (python): {int((bl_frac>=args.blacklist_max_frac).sum())} bins >= {args.blacklist_max_frac:.0%} covered -> dropped")
        master=master[bl_frac<args.blacklist_max_frac].copy()
    if args.min_baseline_quantile>0:
        fb=np.nanquantile(master["baseline_b10"], args.min_baseline_quantile)
        fp=np.nanquantile(master["baseline_parental"], args.min_baseline_quantile)
        master=master[(master["baseline_b10"]>=fb)|(master["baseline_parental"]>=fp)].copy()
    master=master.reset_index(drop=True)
    master.to_csv(os.path.join(args.outdir,"master_bins_H3K27me3_250kb.tsv"), sep="\t", index=False)

    mig=master["mig_dir"]; unc=master["unc_dir_b10"]
    moved=mig.isin(["A_to_B","B_to_A"]) & unc.isin(["A_to_B","B_to_A"])
    print(f"[info] master: {len(master)} bins ({master['is_mig_switcher'].sum()} switchers, "
          f"{master['is_reverter_b10'].sum()} reverters, {int((moved&(mig==unc)).sum())} progressors, "
          f"{master['is_gainer'].sum()} gainers, {master['is_stable'].sum()} stable)")

    reverter_col="is_reverter_b10_strict" if args.reverter_strict else "is_reverter_b10"
    sw=master["is_mig_switcher"].values; bta=(master["mig_dir"]=="B_to_A").values; atb=(master["mig_dir"]=="A_to_B").values
    rev=master[reverter_col].values
    rev_atb=rev & (master["unc_dir_b10"]=="A_to_B").values
    rev_bta=rev & (master["unc_dir_b10"]=="B_to_A").values
    gain=master["is_gainer"].values; stab=master["is_stable"].values
    # compartment-split stable references (same definition as is_stable, split by side)
    stabA=(master["mig_dir"]=="Stable_A").values & (~gain)
    stabB=(master["mig_dir"]=="Stable_B").values & (~gain)
    # gainer compartment split (per arm) for the compartment-matched gainer contrasts.
    # Use the migration-axis PC1 (same coordinate system that defines Stable_A/Stable_B).
    _pc1_par_col,_pc1_b10_col=_pc1_cols(args.migration_axis)
    if _pc1_par_col in master.columns and _pc1_b10_col in master.columns:
        _p1p=pd.to_numeric(master[_pc1_par_col],errors="coerce").values
        _p1b=pd.to_numeric(master[_pc1_b10_col],errors="coerce").values
        gainA_par=gain&(_p1p>0); gainB_par=gain&(_p1p<0)
        gainA_b10=gain&(_p1b>0); gainB_b10=gain&(_p1b<0)
    else:
        gainA_par=gainB_par=gainA_b10=gainB_b10=None

    plan=[("B10","migration_switchers",sw,"dUNC_b10","baseline_b10"),
          ("B10","switchers_A_to_B",atb,"dUNC_b10","baseline_b10"),
          ("B10","switchers_B_to_A",bta,"dUNC_b10","baseline_b10"),
          ("B10","UNC_reverters",rev,"dUNC_b10","baseline_b10"),
          ("B10","reverters_A_to_B",rev_atb,"dUNC_b10","baseline_b10"),
          ("B10","reverters_B_to_A",rev_bta,"dUNC_b10","baseline_b10"),
          ("B10","H3K27me3_gainers",gain,"dUNC_b10","baseline_b10"),
          ("Parental","migration_switchers",sw,"dUNC_parental","baseline_parental"),
          ("Parental","switchers_A_to_B",atb,"dUNC_parental","baseline_parental"),
          ("Parental","switchers_B_to_A",bta,"dUNC_parental","baseline_parental"),
          ("Parental","H3K27me3_gainers",gain,"dUNC_parental","baseline_parental")]

    rows=[]
    for i,(cell,name,tmask,vcol,bcol) in enumerate(plan):
        res=run_contrast(master, tmask, stab, vcol, bcol, args.n_strata, args.n_perm, seed=args.seed+i,
                         n_boot=args.n_boot, baseline_match=args.baseline_match, n_ref_bins=args.n_ref_bins)
        res.update(cell_line=cell, target_set=name, value_col=vcol)
        res["interpretation"]=interpret(res["cliffs_delta"], res["perm_p_two_sided"])
        rows.append(res)
        print(f"  [{cell:8s}] {name:20s} delta={res['cliffs_delta']:+.3f}  p={res['perm_p_two_sided']:.4f}  "
              f"(nT={res['n_target']}, nS={res['n_stable']})  -> {res['interpretation']}")

    results=pd.DataFrame(rows)[["cell_line","target_set","value_col","n_target","n_stable","cliffs_delta",
                                "ci_low","ci_high","perm_p_two_sided","median_target","median_stable",
                                "n_perm_valid","interpretation"]]
    results_path=os.path.join(args.outdir,"enrichment_results.tsv")
    results.to_csv(results_path, sep="\t", index=False)

    lookup={(r["cell_line"],r["target_set"]):r for _,r in results.iterrows()}

    # Panels for ONE arm. Each arm uses its own dUNC column and looks up its own
    # stats. Reverter sets exist only for B10 (reversion is defined on the B10
    # DMSO->UNC axis), so the Parental arm omits them by construction.
    def _panel(label, mask, cat, set_name, cell, value_col, ref=False, style=None):
        mask=np.asarray(mask,bool); n=int(mask.sum())
        if not ref and n<args.min_plot_n: return None
        row=lookup.get((cell,set_name)) if set_name else None
        return {"label":label,"color":cat_color(cat,style),"ref":ref,"n":n,
                "values":master.loc[mask,value_col].dropna().values,
                "delta":(row["cliffs_delta"] if row is not None else np.nan),
                "p":(row["perm_p_two_sided"] if row is not None else np.nan)}

    # --- B10 arm: switchers + reverters + optional gainers ---
    panels_b10=[p for p in [_panel("stable",stab,"stable",None,"B10","dUNC_b10",ref=True),
                        _panel("A→B switch",atb,"switcher","switchers_A_to_B","B10","dUNC_b10"),
                        _panel("B→A switch",bta,"switcher","switchers_B_to_A","B10","dUNC_b10"),
                        _panel("A→B revert",rev_atb,"reverter","reverters_A_to_B","B10","dUNC_b10"),
                        _panel("B→A revert",rev_bta,"reverter","reverters_B_to_A","B10","dUNC_b10")] if p is not None]
    if args.plot_gainers:
        g=_panel("gainer",gain,"gainer","H3K27me3_gainers","B10","dUNC_b10")
        if g is not None: panels_b10.append(g)
    vdesc_b10="dUNC  =  log2(UNC / none)   [B10]"
    plot_ecdf(panels_b10, vdesc_b10, os.path.join(args.outdir,"ecdf_dUNC_b10.png"))
    plot_violin(panels_b10, vdesc_b10, os.path.join(args.outdir,"violin_dUNC_b10.png"), cell_tag="B10")

    # --- Parental arm (migration-naive control): switchers + optional gainers, NO reverters ---
    panels_par=[p for p in [_panel("stable",stab,"stable",None,"Parental","dUNC_parental",ref=True,style="Parental"),
                        _panel("A→B switch",atb,"switcher","switchers_A_to_B","Parental","dUNC_parental",style="Parental"),
                        _panel("B→A switch",bta,"switcher","switchers_B_to_A","Parental","dUNC_parental",style="Parental")] if p is not None]
    if args.plot_gainers:
        g=_panel("gainer",gain,"gainer","H3K27me3_gainers","Parental","dUNC_parental",style="Parental")
        if g is not None: panels_par.append(g)
    vdesc_par="dUNC  =  log2(UNC / none)   [Parental]"
    plot_violin(panels_par, vdesc_par, os.path.join(args.outdir,"violin_dUNC_parental.png"), cell_tag="Parental")

    plot_balance(master, sw, "baseline_b10", "switchers", os.path.join(args.outdir,"balance_switchers.png"))
    plot_forest(results, os.path.join(args.outdir,"forest_effect_sizes.png"))

    # --- compartment-matched selectivity analysis ---
    if args.compartment_matched:
        print("[cmatch] compartment-matched reference analysis (each target vs same-compartment stable):")
        print(f"[cmatch] gainer contrasts: {'ON (--cmatched-gainers given)' if args.cmatched_gainers else 'OFF (opt-in; pass --cmatched-gainers to enable)'}")
        _cmasks=dict(stabA=stabA, stabB=stabB, atb=atb, bta=bta, rev_atb=rev_atb, rev_bta=rev_bta)
        if args.cmatched_gainers and gainA_b10 is not None:
            _cmasks.update(gainA_b10=gainA_b10, gainB_b10=gainB_b10, gainA_par=gainA_par, gainB_par=gainB_par)
        elif args.cmatched_gainers and gainA_b10 is None:
            print(f"[cmatch] gainer split skipped: PC1 columns ({_pc1_par_col},{_pc1_b10_col}) not in master")
        run_compartment_matched(master, _cmasks, args, combined=args.cmatched_combined,
                                gainers=args.cmatched_gainers, per_arm=args.cmatched_per_arm)

    if args.diagnose:
        run_all_diagnostics(master, args, args.outdir)
    if args.gene_table:
        map_genes(master, args.gene_table, args.outdir, genes_per_bin=args.genes_per_bin)
    write_manifest(args.outdir, args, samples, len(master), results_path, norm)
    print(f"[done] results -> {results_path}")


def self_test():
    print("=== SELF-TEST: Q-A enrichment core ===")
    rng=np.random.default_rng(42); N=6000; baseline=rng.normal(0,1,N)
    def make(c, inject, confound=False):
        is_t=np.zeros(N,bool); is_t[rng.choice(N,700,replace=False)]=True
        base=baseline+np.where(is_t,1.2,0.0) if confound else baseline.copy()
        d=rng.normal(0,0.5,N)+c; d[is_t]+=inject
        return d,is_t,base
    def test(values,is_t,base): return stratified_matched_test(values,is_t,base,n_strata=10,n_perm=500,seed=1)
    d,t,b=make(0.0,-0.8); r=test(d,t,b)
    print(f"(a) inject -0.8, c=0   -> delta={r['cliffs_delta']:+.3f} p={r['perm_p_two_sided']:.4f}")
    assert r["cliffs_delta"]<-0.1 and r["perm_p_two_sided"]<0.05
    d2,t2,b2=make(-5.0,-0.8); r2=test(d2,t2,b2)
    print(f"(b) inject -0.8, c=-5  -> delta={r2['cliffs_delta']:+.3f} (offset cancellation)")
    assert r2["cliffs_delta"]<-0.1
    d3,t3,b3=make(-0.3,0.0); r3=test(d3,t3,b3)
    print(f"(c) inject  0.0        -> delta={r3['cliffs_delta']:+.3f} (null)")
    assert abs(r3["cliffs_delta"])<0.1 or r3["perm_p_two_sided"]>=0.05
    d4,t4,b4=make(0.0,0.0,confound=True); r4=test(d4,t4,b4)
    print(f"(d) confound, no effect-> delta={r4['cliffs_delta']:+.3f} (matching null)")
    assert abs(r4["cliffs_delta"])<0.12 or r4["perm_p_two_sided"]>=0.05
    assert cliffs_delta([1,2,3],[4,5,6])<0

    print("--- sign convention + normalization invariance ---")
    # delta<0 must mean "target loses MORE than reference"
    rr=np.random.default_rng(3)
    ref=rr.normal(0.0,0.4,3000); tgt=rr.normal(-0.5,0.4,800)   # target more depleted
    dd=cliffs_delta(tgt,ref)
    print(f"    target more depleted -> delta={dd:+.3f} (must be < 0), tag='{cliff_tag(dd,0.001)}'")
    assert dd<0 and "extra loss" in cliff_tag(dd,0.001)
    dd2=cliffs_delta(ref,tgt)
    assert dd2>0 and "spared" in cliff_tag(dd2,0.001)
    # a GLOBAL shift (what any normalization applies) must not move delta
    shift=1.7
    dd_shift=cliffs_delta(tgt+shift, ref+shift)
    print(f"    +{shift} global shift  -> delta={dd_shift:+.3f} (unchanged: normalization-invariant)")
    assert abs(dd_shift-dd)<1e-9
    # NEGATING the arm flips EVERY delta in it (why 'flip the sign to match biology' is invalid)
    dd_neg=cliffs_delta(-tgt,-ref)
    print(f"    negating the arm     -> delta={dd_neg:+.3f} (flips sign of EVERY set, incl. gainers)")
    assert abs(dd_neg+dd)<1e-9

    print("--- diagnostic: baseline-trend leakage recovery ---")
    rng=np.random.default_rng(7); M=8000; base=rng.normal(0,1,M)
    dunc=-0.7*base+rng.normal(0,0.25,M)            # dUNC obeys a baseline slope, NO real target effect
    is_t=np.zeros(M,bool); prob=(base-base.min())/(base.max()-base.min())
    is_t[rng.random(M)<0.10*(0.3+prob)]=True       # target sits at HIGH baseline edge within strata
    dfm=pd.DataFrame({"v":dunc,"b":base,"is_stable":~is_t})
    dd_dec=run_contrast(dfm,is_t,~is_t,"v","b",n_strata=10,n_perm=300,seed=1,baseline_match="decile")["cliffs_delta"]
    dd_res=run_contrast(dfm,is_t,~is_t,"v","b",n_strata=10,n_perm=300,seed=1,baseline_match="residual",n_ref_bins=50)["cliffs_delta"]
    print(f"    leakage artifact: decile δ={dd_dec:+.3f}  residual δ={dd_res:+.3f}")
    assert abs(dd_res)<abs(dd_dec)*0.6+0.02, "residualization should shrink leakage"
    print("--- compartment-matched: pairing + plumbing ---")
    _test_compartment_matched()
    print("ALL SELF-TESTS PASSED")


def _test_compartment_matched():
    """Synthetic check that run_compartment_matched maps each target to the correct
    same-arm compartment (A_to_B -> Stable_A in Parental but Stable_B in B10, etc.)
    and emits the table + forest + all violins (single-arm, combined, per-arm)."""
    import tempfile
    rng=np.random.default_rng(0); N=6000
    mig=np.random.default_rng(1).choice(["Stable_A","Stable_B","A_to_B","B_to_A","None"],size=N,p=[0.35,0.30,0.12,0.12,0.11])
    pc1_par=np.where(np.isin(mig,["Stable_A","A_to_B"]),0.04,-0.04)+rng.normal(0,0.02,N)
    pc1_b10=np.where(np.isin(mig,["Stable_A","B_to_A"]),0.04,-0.04)+rng.normal(0,0.02,N)
    unc_dir=np.array([None]*N,dtype=object)
    pool=np.where(np.isin(mig,["A_to_B","B_to_A"]))[0]
    for idx in rng.choice(pool,size=int(0.4*len(pool)),replace=False):
        unc_dir[idx]="B_to_A" if mig[idx]=="A_to_B" else "A_to_B"
    master=pd.DataFrame({"mig_dir":mig,"unc_dir_b10":unc_dir,
        "dUNC_parental":-6.0*pc1_par+rng.normal(0,0.3,N),
        "dUNC_b10":-1.2*pc1_b10+rng.normal(0,0.3,N),
        "baseline_parental":6+2*pc1_par+rng.normal(0,0.5,N),
        "baseline_b10":6+2*pc1_b10+rng.normal(0,0.5,N),
        "is_gainer":rng.random(N)<0.08})
    master["is_stable"]=np.isin(mig,["Stable_A","Stable_B"])&~master["is_gainer"]
    gain=master["is_gainer"].values
    masks=dict(stabA=(master["mig_dir"]=="Stable_A").values&~gain,
               stabB=(master["mig_dir"]=="Stable_B").values&~gain,
               atb=(master["mig_dir"]=="A_to_B").values,
               bta=(master["mig_dir"]=="B_to_A").values,
               rev_atb=(master["unc_dir_b10"]=="A_to_B").values,
               rev_bta=(master["unc_dir_b10"]=="B_to_A").values,
               gainA_par=gain&(pc1_par>0), gainB_par=gain&(pc1_par<0),
               gainA_b10=gain&(pc1_b10>0), gainB_b10=gain&(pc1_b10<0))
    # purity assertions (the whole point of the corrected pairing)
    assert (pc1_par[masks["atb"]]>0).mean()>0.9, "A_to_B should be A in Parental"
    assert (pc1_b10[masks["atb"]]>0).mean()<0.1, "A_to_B should be B in B10"
    assert (pc1_par[masks["bta"]]>0).mean()<0.1, "B_to_A should be B in Parental"
    assert (pc1_b10[masks["bta"]]>0).mean()>0.9, "B_to_A should be A in B10"
    assert (pc1_b10[masks["rev_atb"]]>0).mean()>0.9, "rev_atb (UNC A->B, starts A) should be A in B10"
    assert (pc1_b10[masks["rev_bta"]]>0).mean()<0.1, "rev_bta (UNC B->A, starts B) should be B in B10"
    class _A: pass
    a=_A(); a.n_strata=8; a.n_perm=150; a.n_boot=80; a.seed=0
    a.baseline_match="decile"; a.n_ref_bins=40; a.min_plot_n=15
    with tempfile.TemporaryDirectory() as td:
        a.outdir=td
        cm=run_compartment_matched(master, masks, a, combined=True, gainers=True, per_arm=True)
        need=["compartment_matched_results.tsv","forest_effect_sizes_cmatched.png",
              "violin_cmatched_b10_stableA.png","violin_cmatched_b10_stableB.png",
              "violin_cmatched_parental_stableA.png","violin_cmatched_parental_stableB.png",
              "violin_cmatched_combined_stableA.png","violin_cmatched_combined_stableB.png",
              "violin_cmatched_arm_parental.png","violin_cmatched_arm_b10.png"]
        missing=[f for f in need if not os.path.exists(os.path.join(td,f))]
        assert not missing, f"missing cmatched outputs: {missing}"
        assert len(cm)==10 and set(cm["cell_line"])=={"B10","Parental"}, f"expected 10 contrasts, got {len(cm)}"
        assert set(cm["target_set"])>={"gainers_A","gainers_B"}, "gainer contrasts missing"
        # gainers must be looked up per (arm, set): two distinct rows share target_set
        assert len(cm[cm["target_set"]=="gainers_A"])==2, "gainers_A must exist for both arms"
    print("    compartment-matched: 10 contrasts (incl. gainers) + 10 figures OK, pairing verified")


def build_parser():
    p=argparse.ArgumentParser(description="Does UNC1999 selectively erase migration-acquired H3K27me3? "
                              "(scale-free between-set enrichment + compartment-matched selectivity + anomaly diagnostics).",
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--sample-sheet"); p.add_argument("--norm", choices=["cpm","greenlist"], default="cpm")
    p.add_argument("--size-factors", default=None); p.add_argument("--factor", default="H3K27me3")
    p.add_argument("--drop-samples", default="Parental-H3K27me3-R3"); p.add_argument("--compartments")
    p.add_argument("--blacklist", default=None); p.add_argument("--deeptools-blacklist", action="store_true")
    p.add_argument("--blacklist-max-frac", type=float, default=0.5); p.add_argument("--counts-table", default=None)
    p.add_argument("--project-root", default="."); p.add_argument("--outdir", default="results/exploratory/unc_selective_erasure")
    p.add_argument("--threads", type=int, default=8); p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--extend-reads", action="store_true"); p.add_argument("--pseudocount", type=float, default=1.0)
    p.add_argument("--migration-axis", choices=["dmso","orig"], default="dmso")
    p.add_argument("--top-gainer-pct", type=float, default=10.0); p.add_argument("--reverter-strict", action="store_true")
    p.add_argument("--min-baseline-quantile", type=float, default=0.0); p.add_argument("--include-chrx", action="store_true")
    p.add_argument("--gene-table", default=None); p.add_argument("--genes-per-bin", choices=["one","all"], default="one")
    p.add_argument("--plot-gainers", action="store_true"); p.add_argument("--min-plot-n", type=int, default=15)
    p.add_argument("--baseline-match", choices=["decile","fine","residual"], default="decile",
                   help="how the target-vs-stable contrast controls baseline abundance: 'decile' (10 strata, "
                        "original), 'fine' (>=50 strata), or 'residual' (continuous: subtract the stable baseline "
                        "trend, then single-stratum Cliff's delta). If the diagnostics show decile-vs-residual "
                        "collapse, re-run the headline with --baseline-match residual.")
    p.add_argument("--n-ref-bins", type=int, default=50,
                   help="quantile knots for the stable baseline trend (residual mode and D3 diagnostic)")
    p.add_argument("--compartment-matched", action=argparse.BooleanOptionalAction, default=True,
                   help="also run the COMPARTMENT-MATCHED selectivity analysis (each target vs same-compartment "
                        "stable in the arm being measured) into its own forest + violins. --no-compartment-matched "
                        "to skip. This is the selectivity-clean companion to the compartment-MIXED headline "
                        "contrast, and removes the dUNC~PC1 gradient D5 flags.")
    p.add_argument("--cmatched-combined", action=argparse.BooleanOptionalAction, default=True,
                   help="within the compartment-matched block, emit the 2 COMBINED panels (ONE compartment, BOTH "
                        "arms on one axis, each arm keeping its own stable reference; within-arm δ only). "
                        "--no-cmatched-combined to skip.")
    p.add_argument("--cmatched-per-arm", action=argparse.BooleanOptionalAction, default=True,
                   help="within the compartment-matched block, emit the 2 PER-ARM panels (ONE arm, BOTH compartment "
                        "groups side by side: Stable A vs its A-sets, then Stable B vs its B-sets). Because a per-arm "
                        "panel never mixes normalization frames, heights ARE comparable within it. "
                        "--no-cmatched-per-arm to skip.")
    p.add_argument("--cmatched-gainers", action=argparse.BooleanOptionalAction, default=False,
                   help="OPT-IN. Within the compartment-matched block, also test H3K27me3 GAINERS (the most direct "
                        "proxy for migration-installed mark) split by their own-arm compartment vs same-compartment "
                        "stable. Pass --cmatched-gainers to add 4 contrasts (gainers_A/gainers_B x B10/Parental) and a "
                        "gainer violin to each panel. Needs PC1 columns in the compartment file. Off by default, so the "
                        "compartment-matched block emits the 6 switcher/reverter contrasts only.")
    p.add_argument("--diagnose", action=argparse.BooleanOptionalAction, default=True,
                   help="run the anomaly diagnostics suite (D1-D5) into <outdir>/diagnostics/. --no-diagnose to skip.")
    p.add_argument("--n-strata", type=int, default=10); p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--n-boot", type=int, default=1000); p.add_argument("--seed", type=int, default=0)
    return p


if __name__=="__main__":
    args=build_parser().parse_args()
    if args.self_test: self_test(); sys.exit(0)
    miss=[r for r in ["sample_sheet","compartments"] if getattr(args,r) is None]
    if miss: sys.exit(f"[ERROR] missing required args: {miss} (or use --self-test)")
    main(args)