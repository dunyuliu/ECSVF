# PROJECT_RULES.md

## Rule 1 — All verification checks must pass

Every `[PASS]` / `[FAIL]` check emitted by `verify.eqquasi.geometry.conversion.py` must be `[PASS]`.
The five checks are:
- `geometry_round_trip`
- `station_patch_mapping`
- `site_export_consistency`
- `event_summary_consistency`
- `geometry_visualization`

Run for every EQquasi case via `process.eqquasi.sh`. A single `[FAIL]` is a blocking issue.

## Rule 2 — Processing pipeline runs without error

The full pipeline must complete without a non-zero exit code:
```
bash process.eqquasi.sh
python benchmark/utils/process.hbi.py
python benchmark/utils/benchmark.comparison.analytics.py \
  --simulation-root benchmark/Simulation_results \
  --outdir benchmark/analytics
```

PyRSQSim processing (`process.pyrsqsim.sh`) requires the sibling `PyRSQSim.dev` repo and is exempt when that repo is unavailable.

## Rule 3 — Bundle completeness

Every directory under `benchmark/Simulation_results/*_results/` that contains any `*_event_info.csv` file must also contain the matching `*_site1_event_info.csv`, `*_site2_event_info.csv`, and `*_site3_event_info.csv`. Incomplete bundles are blocking.

## Rule 4 — No stale duplicate bundles

No two bundles in the same `*_results/` directory may cover the same model run. If a manually converted legacy bundle (e.g. `alpine_varying_dip_*`) and an auto-generated replacement (e.g. `alpine_hbi_*`) both exist, the legacy bundle must be removed.

## Rule 5 — No hardcoded user-specific absolute paths in scripts

Python scripts and shell wrappers must not contain `/Users/<username>/` paths. Documentation (README.md) uses `${PROJECT_ROOT}/` as the placeholder.
