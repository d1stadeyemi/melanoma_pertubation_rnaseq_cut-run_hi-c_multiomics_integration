import pandas as pd
g   = pd.read_csv("results/exploratory/cut_run_pc1_integration/genelists/track1_gene_assignments.txt", sep="\t")
rev = pd.read_csv("results/exploratory/reversal/Reversal_classified.txt", sep="\t")[
        ["coordinate","mig_switch"]]          # <-- drop reversal_status here
gg = g.merge(rev, on="coordinate", how="left")

sw  = gg["mig_switch"].astype(str).str.contains("_to_")   # A_to_B or B_to_A
rvd = gg["reversal_status"] == "reverted"                 # from g, no collision now

print("Q1 switchers: genes", int(sw.sum()),
      "| of those in csaw-down:", int((sw & g.erode_csaw).sum()),
      "| in delta-loss:", int((sw & g.erode_delta).sum()))
print("Q2 reverted : genes", int(rvd.sum()),
      "| in csaw-down:", int((rvd & g.erode_csaw).sum()),
      "| in delta-loss:", int((rvd & g.erode_delta).sum()))

md  = pd.read_csv("data/processed/mark_signal_250kb/marks_delta_250kb.txt", sep="\t")
rev = pd.read_csv("results/exploratory/reversal/Reversal_classified.txt", sep="\t")
d = rev.merge(md, on="coordinate")
for grp,name in [(d.reversal_status=="reverted","reverted"),
                (d.mig_switch=="A_to_B","A_to_B"),
                (d.mig_switch=="B_to_A","B_to_A"),
                (d.mig_switch=="Stable","stable")]:
    print(name, "median ΔH3K27me3_migration =",
            round(d.loc[grp,"H3K27me3_migration"].median(),3))
    
ab = d[d.mig_switch == "A_to_B"]
revd = ab.reversal_status == "reverted"
print("A→B reverted   median Δ =", round(ab.loc[revd,  "H3K27me3_migration"].median(), 3), "n=", int(revd.sum()))
print("A→B NOT reverted median Δ =", round(ab.loc[~revd, "H3K27me3_migration"].median(), 3), "n=", int((~revd).sum()))

for mk in ["H3K27me3","H3K9me3","H3K4me1","H3K27ac"]:
    col = f"{mk}_migration"
    if col not in d.columns:
        print(mk, "— no column"); continue
    ab = d[d.mig_switch == "A_to_B"]; r = ab.reversal_status == "reverted"
    print(f"{mk:9s} A→B reverted {ab.loc[r,col].median():+.3f}  vs non-reverted {ab.loc[~r,col].median():+.3f}")