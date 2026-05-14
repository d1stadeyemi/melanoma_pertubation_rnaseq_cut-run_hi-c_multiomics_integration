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

## Patches

### `patch_upsetplot.sh` — required for `multiomics-py`

`upsetplot 0.9.0` (the latest available on PyPI) has two bugs that crash any
script using `UpSet.plot()` under the pandas and matplotlib versions in this
env. The patch script monkey-patches the installed library file in place.

**You must run this once after creating or updating the env:**

```bash
conda activate multiomics-py
bash envs/patch_upsetplot.sh
```

The script is **idempotent** — re-running on an already-patched install is a
no-op. It also creates a timestamped backup of the original file before
modifying it.

#### What the patch fixes

**Bug 1 — `ValueError: Invalid RGBA argument: nan`**

In pandas ≥3.0, Copy-on-Write is always on. Calling
`styles["edgecolor"].fillna(..., inplace=True)` operates on a throwaway copy of
the column, so `edgecolor` stays as `NaN`. When matplotlib later tries to draw
a dot with color `NaN`, it raises.

**Bug 2 — `TypeError: only 0-dimensional arrays can be converted to Python scalars`**

`np.diff(ax.get_xlim())` returns a 1-element numpy array, not a scalar. Newer
matplotlib rejects this when used as a text coordinate.

#### Manual diff (for the record)

Bug 1 — in the styles-fill block of `upsetplot/plotting.py` (~line 795):

```diff
- styles["linewidth"].fillna(1, inplace=True)
- styles["facecolor"].fillna(self._facecolor, inplace=True)
- styles["edgecolor"].fillna(styles["facecolor"], inplace=True)
- styles["linestyle"].fillna("solid", inplace=True)
+ styles["linewidth"] = styles["linewidth"].fillna(1)
+ styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)
+ styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])
+ styles["linestyle"] = styles["linestyle"].fillna("solid")
```

Bug 2 — in `_label_sizes` (~lines 904–926, three occurrences total):

```diff
- margin = 0.01 * abs(np.diff(ax.get_xlim()))
+ margin = 0.01 * abs(np.diff(ax.get_xlim())).item()
```

Both `get_xlim()` instances and the one `get_ylim()` instance.

#### Affected analysis scripts

- `scripts/analysis/trimodal_coupling_flex.py` (Fig 2a UpSet plot)

#### When to re-run the patch

- After `conda env create -f envs/multiomics-py.yml`
- After `conda env update -f envs/multiomics-py.yml`
- After any `pip install --force-reinstall upsetplot`
- After cloning the repo on a new machine

#### How to verify the patch is active

```bash
conda activate multiomics-py
PLOTTING=$(pip show upsetplot | awk '/^Location:/ {print $2}')/upsetplot/plotting.py
grep -c 'styles\[.*\]\.fillna.*inplace=True' "$PLOTTING"   # expect: 0
grep -c '\.item()' "$PLOTTING"                              # expect: ≥3
```