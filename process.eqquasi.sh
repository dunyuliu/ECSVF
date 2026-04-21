#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MPLCONFIGDIR="${ROOT_DIR}/benchmark/.mplconfig"

PLANAR_CASE="${ROOT_DIR}/results/nz.bp5.qdc.dip50.2000.norm_6mm_yr"
PLANAR_GEOM="${ROOT_DIR}/geometry/processed_geometry/202603_alpine_planar_dip50_case.res2km/alpine_planar_dip50_case.res2km_grid.npz"
NO_DIP_CASE="${ROOT_DIR}/results/nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad"
VARY_DIP_CASE="${ROOT_DIR}/results/nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad"
EQQ_OUTPUT="${ROOT_DIR}/benchmark/Simulation_results/EQquasi_results"

rm -rf "${EQQ_OUTPUT}"

python "${ROOT_DIR}/benchmark/utils/process.eqquasi.py" \
  --input-dir "${PLANAR_CASE}" \
  --geometry-npz "${PLANAR_GEOM}" \
  --output-dir "${EQQ_OUTPUT}" \
  --prefix alpine_planar \
  --slip-rate-threshold 0.1

python "${ROOT_DIR}/benchmark/utils/verify.eqquasi.geometry.conversion.py" \
  --case-dir "${PLANAR_CASE}" \
  --geometry-npz "${PLANAR_GEOM}" \
  --export-dir "${EQQ_OUTPUT}" \
  --prefix alpine_planar \
  --slip-rate-threshold 0.1

python "${ROOT_DIR}/benchmark/utils/process.eqquasi.py" \
  --input-dir "${NO_DIP_CASE}" \
  --output-dir "${EQQ_OUTPUT}" \
  --slip-rate-threshold 0.1

python "${ROOT_DIR}/benchmark/utils/verify.eqquasi.geometry.conversion.py" \
  --case-dir "${NO_DIP_CASE}" \
  --export-dir "${EQQ_OUTPUT}" \
  --prefix alpine_no_dip_change \
  --slip-rate-threshold 0.1

python "${ROOT_DIR}/benchmark/utils/process.eqquasi.py" \
  --input-dir "${VARY_DIP_CASE}" \
  --output-dir "${EQQ_OUTPUT}" \
  --slip-rate-threshold 0.1

python "${ROOT_DIR}/benchmark/utils/verify.eqquasi.geometry.conversion.py" \
  --case-dir "${VARY_DIP_CASE}" \
  --export-dir "${EQQ_OUTPUT}" \
  --prefix alpine_varying_dip \
  --slip-rate-threshold 0.1
