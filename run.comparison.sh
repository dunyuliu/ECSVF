#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MPLCONFIGDIR="${ROOT_DIR}/benchmark/.mplconfig"

ANALYTICS_OUTPUT="${ROOT_DIR}/benchmark/analytics"

rm -rf "${ANALYTICS_OUTPUT}"

bash "${ROOT_DIR}/process.eqquasi.sh"
bash "${ROOT_DIR}/process.pyrsqsim.sh"

python "${ROOT_DIR}/benchmark/utils/benchmark.comparison.analytics.py" \
  --simulation-root "${ROOT_DIR}/benchmark/Simulation_results" \
  --outdir "${ANALYTICS_OUTPUT}"
