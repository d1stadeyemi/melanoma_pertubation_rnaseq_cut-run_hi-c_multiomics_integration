#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=logs/integrate_%j.out
#SBATCH --error=logs/integrate_%j.err
#SBATCH --mail-user=majagbe@vols.utk.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration

module load anaconda3/2024.06
source $ANACONDA3_SH

conda activate multiomics-py

python scripts/analysis/integrate_hic_cutrun_rnaseq.py \
  --hic data/processed/compartment_change_18AUG2026.txt \
  --hic_coord_col Coord --hic_chr_col Chr \
  --hic_delta_col "Bottom10_UNC-DMSO" \
  --hic_category_col "Bottom10_UNC-DMSO_category" \
  --restrict_to "Switch AtoB" "Switch BtoA" \
  --cutrun data/processed/cut_run_signal_table_250kb.eps_fixed.txt \
  --cutrun_coord_col extra_1 \
  --mark_cols log2ratio_p-b10_h3k27me3_mean \
  --mark_names H3K27me3 \
  --deseq2 /lustre/isaac24/proj/UTK0009/Taiwo/2026JUL-Taiwo-A375UNC-RNA/New_Deseq2_Output/BU_vs_BD/tables/20_results_Bottom-10-UNC_vs_Bottom-10-DMSO_full.csv \
  --deseq2_lfc_col log2FoldChange_shrunk \
  --gene_map data/processed/GeneListTSShg19_2kb_up_5kb_down.txt \
  --deseq2_gene_col gene_name \
  --outdir results/exploratory/integration_hic_cutrun_rnaseq/

# The code below is used to check how many of the DE genes 
# are mapped to the gene list used for integration. 
# It is not part of the main analysis, but it provides 
#useful information about the overlap between the DE genes and the gene list.

python - << 'PY'
import pandas as pd

DE = "/lustre/isaac24/proj/UTK0009/Taiwo/2026JUL-Taiwo-A375UNC-RNA/New_Deseq2_Output/BU_vs_BD/tables/20_results_Bottom-10-UNC_vs_Bottom-10-DMSO_full.csv"
GM = "data/processed/GeneListTSShg19_2kb_up_5kb_down.txt"

de = pd.read_csv(DE)
gm = pd.read_csv(GM, sep="\t")

mapped = set(gm["gene_name"].astype(str).str.strip())
hit  = de[de.gene_name.astype(str).str.strip().isin(mapped)]
miss = de[~de.gene_name.astype(str).str.strip().isin(mapped)]

print(f"mapped {len(hit):,} / {len(de):,}   unmapped {len(miss):,} ({len(miss)/len(de):.1%})\n")

if "biotype" in de.columns:
    print("unmapped by biotype:")
    print(miss.biotype.value_counts().head(10).to_string(), "\n")
    print("mapped by biotype:")
    print(hit.biotype.value_counts().head(10).to_string(), "\n")

col = "log2FoldChange_shrunk"
print(f"median |{col}|   unmapped {miss[col].abs().median():.3f}"
      f"   mapped {hit[col].abs().median():.3f}")
if "padj" in de.columns:
    print(f"fraction padj<0.05   unmapped {(miss.padj<0.05).mean():.1%}"
          f"   mapped {(hit.padj<0.05).mean():.1%}")

if "chrom" in de.columns:
    print("\nunmapped by chromosome (top 8):")
    print(miss.chrom.value_counts().head(8).to_string())
PY

# The code below is used to calculate the Spearman correlation between the bin categories.

python - << 'PY'
import pandas as pd
from scipy.stats import spearmanr
d = pd.read_csv("results/exploratory/integration_hic_cutrun_rnaseq/integration_bin_table.txt", sep="\t")
m = "log2ratio_p-b10_h3k27me3_mean"
print("all:", spearmanr(d[m], d.expr).statistic.round(3), len(d))
k = d[d[m] > -1.5]
print("excluding the far-left tail:", spearmanr(k[m], k.expr).statistic.round(3), len(k))
for g, s in d.groupby("_cat_"):
    print(f"  within {g}: rho={spearmanr(s[m], s.expr).statistic:+.3f}  n={len(s)}")
PY