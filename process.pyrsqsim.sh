#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MPLCONFIGDIR="${ROOT_DIR}/benchmark/.mplconfig"

PYRSQSIM_WORK_DIR="${ROOT_DIR}/../PyRSQSim.dev/work/alpine_vary_dip"
PYRSQSIM_RESULTS_DIR="${PYRSQSIM_WORK_DIR}/results"
PYRSQSIM_STL="${PYRSQSIM_WORK_DIR}/variable_dip_5km.stl"
PYRSQSIM_OUTPUT="${ROOT_DIR}/benchmark/Simulation_results/PyRSQSim_results"

rm -rf "${PYRSQSIM_OUTPUT}"

python "${ROOT_DIR}/benchmark/utils/process.pyrsqsim.py" \
  --results-dir "${PYRSQSIM_RESULTS_DIR}" \
  --stl "${PYRSQSIM_STL}" \
  --output-dir "${PYRSQSIM_OUTPUT}" \
  --prefix alpine_varying_dip_5km \
  --hypocenter-policy onset_centroid
