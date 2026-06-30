#!/bin/bash
#SBATCH -J net_check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --output=net_check_%j.out
#SBATCH --error=net_check_%j.err

###############################################################################
# Compute-node internet connectivity check.
#
# Submit this as a JOB so it runs on a compute node (the whole point is to test
# the compute node, NOT the login node):
#     sbatch check_internet.sh
#
# Then read the report:
#     cat net_check_<jobid>.out
#
# Each test uses a short timeout so blocked endpoints fail fast instead of
# hanging. A reachable endpoint prints its HTTP status line; a blocked one
# prints FAILED.
###############################################################################

TIMEOUT=8   # seconds per test; blocked endpoints give up after this

echo "============================================================"
echo " ISAAC compute-node connectivity report"
echo " Node : $(hostname)"
echo " Date : $(date)"
echo " Job  : ${SLURM_JOB_ID:-NA}  Partition: ${SLURM_JOB_PARTITION:-NA}"
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# 0. Proxy configuration (if a proxy is set, outbound traffic may route via it)
# ---------------------------------------------------------------------------
echo "--- Proxy environment variables ---"
proxy_vars=$(env | grep -i proxy || true)
if [[ -n "${proxy_vars}" ]]; then
    echo "${proxy_vars}"
    echo "NOTE: a proxy is configured — outbound internet likely goes through it."
else
    echo "(none set)"
fi
echo ""

# ---------------------------------------------------------------------------
# Helper: test an HTTPS endpoint, print its status line or FAILED
# ---------------------------------------------------------------------------
test_url () {
    local label="$1"
    local url="$2"
    local status
    # -s silent, -I headers only, -L follow redirects, --max-time hard timeout
    status=$(curl -sIL --max-time "${TIMEOUT}" "${url}" 2>/dev/null | head -1 | tr -d '\r')
    if [[ -n "${status}" ]]; then
        printf "  %-28s REACHABLE   (%s)\n" "${label}" "${status}"
        return 0
    else
        printf "  %-28s FAILED\n" "${label}"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# 1. DNS resolution (separate from HTTP — DNS can work while HTTP is blocked)
# ---------------------------------------------------------------------------
echo "--- DNS resolution ---"
if getent hosts github.com >/dev/null 2>&1; then
    echo "  DNS resolves github.com      OK"
else
    echo "  DNS resolves github.com      FAILED (no DNS on this node)"
fi
echo ""

# ---------------------------------------------------------------------------
# 2. General internet (HTTPS)
# ---------------------------------------------------------------------------
echo "--- General HTTPS endpoints ---"
general_ok=0
test_url "Google"            "https://www.google.com"      && general_ok=$((general_ok+1))
test_url "GitHub"            "https://github.com"          && general_ok=$((general_ok+1))
test_url "GitHub raw"        "https://raw.githubusercontent.com" && general_ok=$((general_ok+1))
test_url "PyPI"              "https://pypi.org"            && general_ok=$((general_ok+1))
echo ""

# ---------------------------------------------------------------------------
# 3. Bioinformatics-relevant endpoints (the ones your pipelines actually hit)
# ---------------------------------------------------------------------------
echo "--- Bioinformatics endpoints ---"
bio_ok=0
test_url "Ensembl FTP"       "https://ftp.ensembl.org"            && bio_ok=$((bio_ok+1))
test_url "Ensembl REST"      "https://rest.ensembl.org/info/ping" && bio_ok=$((bio_ok+1))
test_url "NCBI FTP"          "https://ftp.ncbi.nlm.nih.gov"       && bio_ok=$((bio_ok+1))
test_url "UCSC downloads"    "https://hgdownload.soe.ucsc.edu"    && bio_ok=$((bio_ok+1))
test_url "Bioconductor"      "https://bioconductor.org"           && bio_ok=$((bio_ok+1))
test_url "KEGG REST"         "https://rest.kegg.jp"               && bio_ok=$((bio_ok+1))
echo ""

# ---------------------------------------------------------------------------
# 4. Conda/package channels (relevant if you ever install from a job)
# ---------------------------------------------------------------------------
echo "--- Package channels ---"
test_url "conda-forge (anaconda)" "https://conda.anaconda.org/conda-forge/noarch/repodata.json"
test_url "bioconda (anaconda)"    "https://conda.anaconda.org/bioconda/noarch/repodata.json"
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "============================================================"
echo " SUMMARY"
echo "   General HTTPS reachable : ${general_ok} / 4"
echo "   Bioinformatics reachable: ${bio_ok} / 6"
if [[ ${general_ok} -eq 0 && ${bio_ok} -eq 0 ]]; then
    echo "   VERDICT: NO outbound internet from this compute node."
    echo "            -> Pre-download all references/data on the LOGIN node,"
    echo "               then have jobs read from local project storage."
elif [[ ${general_ok} -gt 0 || ${bio_ok} -gt 0 ]]; then
    echo "   VERDICT: SOME / FULL outbound internet available."
    echo "            -> Note which specific endpoints FAILED above; those"
    echo "               services must be pre-fetched on the login node."
fi
echo "============================================================"