#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# patch_upsetplot.sh
#
# Patch upsetplot 0.9.0 to work with pandas ≥3.0 (Copy-on-Write) and
# matplotlib ≥3.9. Without this, the UpSet plots in
# scripts/analysis/trimodal_coupling_flex.py crash.
#
# See envs/README.md for the full rationale.
#
# USAGE:
#   conda activate multiomics-py
#   bash envs/patch_upsetplot.sh
#
# This script is idempotent — safe to re-run. It will detect an already-patched
# install and exit cleanly.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Sanity check: are we in the right conda env? ─────────────────────────────
if [[ "${CONDA_DEFAULT_ENV:-}" != "multiomics-py" ]]; then
    echo "ERROR: activate the multiomics-py env first:" >&2
    echo "         conda activate multiomics-py" >&2
    exit 1
fi

# ── Locate upsetplot's plotting.py without importing it ──────────────────────
# (Import via 'python -c' is slow on Lustre; pip show is instant.)
UPSET_SITE=$(pip show upsetplot 2>/dev/null | awk '/^Location:/ {print $2}')
if [[ -z "${UPSET_SITE}" ]]; then
    echo "ERROR: upsetplot not found in this env. Did pip install fail?" >&2
    exit 1
fi
PLOTTING="${UPSET_SITE}/upsetplot/plotting.py"

if [[ ! -f "${PLOTTING}" ]]; then
    echo "ERROR: expected file not found: ${PLOTTING}" >&2
    exit 1
fi

# ── Version check ────────────────────────────────────────────────────────────
UPSET_VERSION=$(pip show upsetplot | awk '/^Version:/ {print $2}')
if [[ "${UPSET_VERSION}" != "0.9.0" ]]; then
    echo "WARNING: patch was written for upsetplot 0.9.0, found ${UPSET_VERSION}." >&2
    echo "         Inspect envs/README.md and verify the patch targets still apply." >&2
    echo "         Continuing anyway in 5 seconds — Ctrl-C to abort..."
    sleep 5
fi

# ── Idempotency check ────────────────────────────────────────────────────────
# Bug 1's signature: presence of inplace=True on the styles fill block.
# If absent, the patch has already been applied.
if ! grep -q 'styles\["linewidth"\]\.fillna(1, inplace=True)' "${PLOTTING}"; then
    echo "[OK] upsetplot already patched. Nothing to do."
    exit 0
fi

# ── Back up before modifying ─────────────────────────────────────────────────
BACKUP="${PLOTTING}.bak.$(date +%Y%m%d_%H%M%S)"
cp "${PLOTTING}" "${BACKUP}"
echo "[INFO] Backup written to: ${BACKUP}"

# ── BUG 1: inplace=True silently fails under pandas Copy-on-Write ────────────
sed -i \
  's/styles\["linewidth"\]\.fillna(1, inplace=True)/styles["linewidth"] = styles["linewidth"].fillna(1)/' \
  "${PLOTTING}"

sed -i \
  's/styles\["facecolor"\]\.fillna(self\._facecolor, inplace=True)/styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)/' \
  "${PLOTTING}"

sed -i \
  's/styles\["edgecolor"\]\.fillna(styles\["facecolor"\], inplace=True)/styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])/' \
  "${PLOTTING}"

sed -i \
  's/styles\["linestyle"\]\.fillna("solid", inplace=True)/styles["linestyle"] = styles["linestyle"].fillna("solid")/' \
  "${PLOTTING}"

# ── BUG 2: np.diff() returns 1-element array, not scalar ─────────────────────
sed -i \
  's/np\.diff(ax\.get_xlim()))/np.diff(ax.get_xlim())).item()/g' \
  "${PLOTTING}"

sed -i \
  's/np\.diff(ax\.get_ylim()))/np.diff(ax.get_ylim())).item()/g' \
  "${PLOTTING}"

# ── Clear cached bytecode so the new source is actually loaded ──────────────
find "${UPSET_SITE}/upsetplot" -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

# ── Verify ──────────────────────────────────────────────────────────────────
REMAINING_INPLACE=$(grep -c 'styles\[.*\]\.fillna.*inplace=True' "${PLOTTING}" || true)
ITEM_HITS=$(grep -c '\.item()' "${PLOTTING}" || true)

echo
echo "[VERIFY] inplace=True style hits remaining: ${REMAINING_INPLACE} (expected 0)"
echo "[VERIFY] .item() occurrences:                ${ITEM_HITS} (expected ≥ 3)"

if [[ "${REMAINING_INPLACE}" != "0" || "${ITEM_HITS}" -lt 3 ]]; then
    echo "ERROR: patch verification failed. Restoring backup." >&2
    cp "${BACKUP}" "${PLOTTING}"
    exit 1
fi

echo
echo "[SUCCESS] upsetplot ${UPSET_VERSION} patched in: ${PLOTTING}"
echo "          Backup retained at: ${BACKUP}"