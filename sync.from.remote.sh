#!/usr/bin/env bash
# Pull the latest git commits and sync EQquasi simulation cases from Knox.
# Usage:
#   bash sync.from.remote.sh                   # uses UTIG_HOST env var or default
#   bash sync.from.remote.sh dliu@knox.ig.utexas.edu

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_HOST="${1:-${UTIG_HOST:-dliu@knox.ig.utexas.edu}}"
REMOTE_BASE="/home/staff/dliu/0.Dunyu/EQquasi/nz.work"

# Both benchmark cases live on Knox
CASES=(
    "nz.bp5.qdc.dip50.2000.norm_6mm_yr"
    "nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad"
)

# ── 1. git pull ────────────────────────────────────────────────────────────────
echo "=== git pull ==="
git -C "${ROOT_DIR}" pull

# ── 2. rsync each case from Knox, then zip locally ────────────────────────────
for CASE in "${CASES[@]}"; do
    LOCAL_PATH="${ROOT_DIR}/results/${CASE}"
    LOCAL_ZIP="${ROOT_DIR}/results/${CASE}.zip"

    mkdir -p "${LOCAL_PATH}"
    echo ""
    echo "=== Syncing ${CASE} ==="
    rsync -av --progress "${REMOTE_HOST}:${REMOTE_BASE}/${CASE}/" "${LOCAL_PATH}/"

    echo "=== Zipping ${CASE} ==="
    zip -r "${LOCAL_ZIP}" "${LOCAL_PATH}"
    echo "  -> ${LOCAL_ZIP}"
done

echo ""
echo "=== Done. Run bash process.eqquasi.sh to rebuild benchmark bundles. ==="
echo ""
echo "Upload zips to Google Drive:"
echo "  https://drive.google.com/drive/folders/11oD23J5W3nmJigIHYhmviJK0_77hvalL"
