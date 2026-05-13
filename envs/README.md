# Conda environments

Three environments cover this project:

| Env              | Purpose                                          | Used by                                       |
|------------------|--------------------------------------------------|-----------------------------------------------|
| `multiomics-py`  | Python data wrangling, bigwig, plotting          | most `scripts/preprocessing/`, `scripts/analysis/`, `scripts/figures/*.py` |
| `multiomics-ml`  | ML, SHAP, deep learning                          | `scripts/analysis/RF_SHAP.py`, future ML work |
| `multiomics-r`   | R + Bioconductor (mixOmics, DESeq2, circlize)    | `scripts/figures/circos_plot.R`, `scripts/figures/PLSDA.R` |

## Setup (first time)

```bash
module load anaconda3/2024.06
source $ANACONDA3_SH
conda env create -f envs/multiomics-py.yml
conda env create -f envs/multiomics-ml.yml
conda env create -f envs/multiomics-r.yml
```

## Exact reproduction (from lock files)

```bash
conda env create -f envs/multiomics-py.lock.yml
```

## Activating

```bash
conda activate multiomics-py    # or -ml, or -r
```

## Installing an extra tool into an existing env

```bash
conda install -c bioconda <tool-name>
# or for a Python package not on conda:
pip install <package-name>
```

### Then update your env files so the addition is reproducible:
### Edit envs/multiomics-py.yml and add the tool under dependencies: (or under pip: if you used pip)
### Keep it unpinned for now, just the name.

## Updating after editing a .yml

```bash
conda env update -f envs/multiomics-py.yml --prune
conda env export --no-builds | grep -v "^prefix:" > envs/multiomics-py.lock.yml
```
## Quick way to find if a tool is available and where:

```bash
conda search -c bioconda <tool-name>
conda search -c conda-forge <tool-name>
# If neither channel has it, try
pip search <tool-name>
```

## Steps to follow when creating a new dedicated env (e.g., for MultiQC):

```bash
cat > envs/multiqc.yml << 'EOF'
name: multiqc
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - multiqc
EOF

conda env create -f envs/multiqc.yml
conda activate multiqc
multiqc --version

conda env export --no-builds | grep -v "^prefix:" > envs/multiqc.lock.yml
conda deactivate

git add envs/multiqc.yml envs/multiqc.lock.yml
git commit -m "Add multiqc env"
```

## Creating an env for Snakemake

```bash
cat > envs/workflow.yml << 'EOF'
name: workflow
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - snakemake-minimal
  - snakemake-wrapper-utils
  - pulp
EOF

conda env create -f envs/workflow.yml
```