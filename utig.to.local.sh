#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_NAME="nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad"
REMOTE_HOST="${1:-${UTIG_HOST:-dliu@knox.ig.utexas.edu}}"
REMOTE_PATH="/home/staff/dliu/0.Dunyu/EQquasi/nz.work/${MODEL_NAME}/"
LOCAL_PATH="${ROOT_DIR}/results/${MODEL_NAME}/"

mkdir -p "${LOCAL_PATH}"

echo "Syncing ${REMOTE_HOST}:${REMOTE_PATH} -> ${LOCAL_PATH}"
rsync -av --progress "${REMOTE_HOST}:${REMOTE_PATH}" "${LOCAL_PATH}"
