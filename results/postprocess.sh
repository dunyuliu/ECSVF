#!/usr/bin/env bash

set -euo pipefail

# Post-process an EQquasi case: run all plotting scripts and save to a.results/
#
# Usage:
#   bash postprocess.sh <case_folder_name>
#
# Example:
#   bash postprocess.sh nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_NAME="${1:?Usage: postprocess.sh <case_folder_name>}"
CASE_DIR="${SCRIPT_DIR}/${CASE_NAME}"

if [ ! -d "${CASE_DIR}" ]; then
    echo "Error: case directory not found: ${CASE_DIR}"
    exit 1
fi

# Count cycles
NCYC=$(ls -d "${CASE_DIR}"/Q* 2>/dev/null | wc -l | tr -d ' ')
if [ "${NCYC}" -eq 0 ]; then
    echo "Error: no Q* cycle folders found in ${CASE_DIR}"
    exit 1
fi
echo "Case: ${CASE_NAME}"
echo "Cycles: ${NCYC}"

# Create output folder
OUTDIR="${CASE_DIR}/a.results"
rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}"

export PYTHONPATH="${CASE_DIR}:${PYTHONPATH:-}"

# --- 1. post.process.eq.stats.py (fast) ---
echo ""
echo "=== post.process.eq.stats.py ==="
(cd "${CASE_DIR}" && python "${SCRIPT_DIR}/post.process.eq.stats.py")
[ -f "${CASE_DIR}/earthquake_stats.csv" ] && mv "${CASE_DIR}/earthquake_stats.csv" "${OUTDIR}/"

# --- 2. plotSlipAtPaleoSite: station at 50% along strike (moderate) ---
echo ""
echo "=== plotSlipAtPaleoSite ==="
(cd "${CASE_DIR}" && python "${SCRIPT_DIR}/plotSlipAtPaleoSite" 50 "${NCYC}" 1)
mv "${CASE_DIR}"/c.paleosite_slip_*.png "${OUTDIR}/" 2>/dev/null || true

# --- 3. plotAccumulated: all variables, both modes, surface (100%) (slow) ---
echo ""
echo "=== plotAccumulated ==="
(cd "${CASE_DIR}" && python "${SCRIPT_DIR}/plotAccumulated" all 0 "${NCYC}" 1 88)
mv "${CASE_DIR}"/cVerticalProfile.*.png "${OUTDIR}/" 2>/dev/null || true
mv "${CASE_DIR}"/cHorizontalProfile.*.png "${OUTDIR}/" 2>/dev/null || true

# --- 4. plotOnFaultVars: run in each Q* folder (slowest) ---
echo ""
echo "=== plotOnFaultVars ==="
for q in $(ls -d "${CASE_DIR}"/Q* | sort -V); do
    qname=$(basename "${q}")
    echo "  Processing ${qname}..."
    (cd "${q}" && python "${SCRIPT_DIR}/plotOnFaultVars")
    [ -f "${q}/on_fault_vars.gif" ] && cp "${q}/on_fault_vars.gif" "${OUTDIR}/${qname}_on_fault_vars.gif"
done

echo ""
echo "=== Done ==="
echo "Results in: ${OUTDIR}"
ls "${OUTDIR}"
