# Benchmark

This folder is the cross-code benchmark layer for the Alpine Fault intercomparison. The root project is your `EQquasi` workspace; this folder only handles format conversion, verification, and comparison against other simulators.

All benchmark scripts follow one rule: no fallback behavior. If a required file, format, case, or bundle is missing or inconsistent, the script stops with an error instead of silently switching to an alternate path.

The current source of truth for the benchmark format is the latest local draft:

- `benchmark/docs/Earthquake Simulators Intercomparison file format_v0.1`

That PDF currently states:

- File 1 uses `Mw > 7.0` in the detailed section
- rupture area uses `slip > 0 m`
- southern and northern rupture bounds use `slip > 0.1 m`
- site files include all earthquakes that rupture the patch closest to the site

> **Implementation note:** the EQquasi exporter uses `0.5 m` for the rupture bound threshold (overriding the draft `0.1 m`). This suppresses sub-threshold noise fragments at the domain edges that would otherwise inflate rupture extents to nearly full-fault length. Rupture bounds are also computed from the **surface row** of the fault mesh only, so that the down-dip displacement of deep patches on the dipping fault does not bias the along-strike northing bounds.

## Layout

- `Simulation_results/`: Avi-style benchmark bundles grouped by code. **Local-only, not tracked in git.**
- `analytics/`: generated summary tables and plots. **Local-only, not tracked in git.**
- `utils/process.eqquasi.py`: export one `EQquasi` case into an Avi-style bundle.
- `utils/process.pyrsqsim.py`: export the sibling `PyRSQSim.dev/examples/alpine` run directly from its STL, catalog, and history files.
- `utils/verify.eqquasi.geometry.conversion.py`: verify geometry back-transform, surface-site mapping, site slips, and event summary values by recomputing them from the raw `Q*` outputs.
- `utils/benchmark.comparison.analytics.py`: compare every complete bundle found under `Simulation_results/*_results`.
- `utils/plot.accumulated.slip.eqquasi.py`: standalone script that produces three publication-quality figures for a single `EQquasi` case — individual earthquake slip profiles (colored by time), running accumulated coseismic slip, and long-term slip rate vs. along-strike distance. Reads `fault.r.nc` from each `Q*` folder at 88% depth; applies the 0.5 m slip threshold. Accepts `--case-dir` and optional `--outdir`.
- `/Users/dliu/scratch/alpine/process.eqquasi.sh`: shared shell wrapper for rebuilding the `EQquasi` benchmark bundles.
- `/Users/dliu/scratch/alpine/process.pyrsqsim.sh`: shared shell wrapper for rebuilding the `PyRSQSim` Alpine example bundle.
- `/Users/dliu/scratch/alpine/run.comparison.sh`: shared shell runner for rebuilding benchmark bundles and analytics together.

Each complete bundle must contain exactly four CSVs:

- `<prefix>_event_info.csv`
- `<prefix>_site1_event_info.csv`
- `<prefix>_site2_event_info.csv`
- `<prefix>_site3_event_info.csv`

If any required file is missing, the scripts stop instead of skipping it.

## Export PyRSQSim Alpine Example

The `PyRSQSim` converter works directly from the sibling example run:

- `/Users/dliu/scratch/PyRSQSim.dev/examples/alpine`

It reads:

- the example STL mesh
- `earthquakes.csv`
- `rsfsim_history.npz`
- `summary_data.npz`

and reconstructs each event from the saved history arrays without any mesh-coordinate conversion step.

Example:

```bash
bash /Users/dliu/scratch/alpine/process.pyrsqsim.sh
```

The shared shell wrapper writes to:

- `Simulation_results/PyRSQSim_results`

and explicitly uses:

```bash
--hypocenter-policy onset_centroid
```

because the current Alpine example does not save a unique first-rupturing element for many large events.

The Python converter itself stays strict by default:

- site mapping uses the nearest shallowest STL element centroid to each benchmark station
- event slip is reconstructed from the history arrays, not from the summary PNG/NPZ plots
- `--hypocenter-policy strict` is the default
- if a large event has an ambiguous onset set rather than a unique first-rupturing element, the script stops instead of inventing a hypocentre

If you explicitly want a proxy hypocentre for such cases, pass:

```bash
--hypocenter-policy onset_centroid
```

## Export EQquasi

Case naming for this benchmark layer:

- `nz.bp5.qdc.dip50.2000.norm_6mm_yr` exports as `alpine_planar` — uses `geometry/raw_geometry/post_workshop_geoms/alpine_planar.stl`
- `nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad` exports as `alpine_no_dip_change`
- `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad` exports as `alpine_varying_dip` — uses `geometry/raw_geometry/post_workshop_geoms/variable_dip_1km.stl`
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad` is the older varying-dip case and is no longer the shared benchmark source

Planar example:

```bash
python /Users/dliu/scratch/alpine/benchmark/utils/process.eqquasi.py \
  --input-dir /Users/dliu/scratch/alpine/results/nz.bp5.qdc.dip50.2000.norm_6mm_yr \
  --geometry-npz /Users/dliu/scratch/alpine/geometry/processed_geometry/202603_alpine_planar_dip50_case.res2km/alpine_planar_dip50_case.res2km_grid.npz \
  --output-dir /Users/dliu/scratch/alpine/benchmark/Simulation_results/EQquasi_results \
  --prefix alpine_planar \
  --slip-rate-threshold 0.1
```

This writes:

- `Simulation_results/EQquasi_results/alpine_planar_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_planar_site1_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_planar_site2_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_planar_site3_event_info.csv`

No-dip-change example:

```bash
python /Users/dliu/scratch/alpine/benchmark/utils/process.eqquasi.py \
  --input-dir /Users/dliu/scratch/alpine/results/nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad \
  --output-dir /Users/dliu/scratch/alpine/benchmark/Simulation_results/EQquasi_results \
  --slip-rate-threshold 0.1
```

This writes:

- `Simulation_results/EQquasi_results/alpine_no_dip_change_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_no_dip_change_site1_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_no_dip_change_site2_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_no_dip_change_site3_event_info.csv`

Variable-dip example:

```bash
python /Users/dliu/scratch/alpine/benchmark/utils/process.eqquasi.py \
  --input-dir /Users/dliu/scratch/alpine/results/nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad \
  --output-dir /Users/dliu/scratch/alpine/benchmark/Simulation_results/EQquasi_results \
  --slip-rate-threshold 0.1
```

This writes:

- `Simulation_results/EQquasi_results/alpine_varying_dip_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_varying_dip_site1_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_varying_dip_site2_event_info.csv`
- `Simulation_results/EQquasi_results/alpine_varying_dip_site3_event_info.csv`

You can override the bundle name with `--prefix`, and the paths with `--input-dir` and `--output-dir`.

Defaults in the exporter:

- File 1 threshold: `Mw > 7.0`
- rupture area threshold: `slip > 0.0 m`
- rupture bound threshold: `slip > 0.5 m` (overrides draft `0.1 m` to suppress noise fragments)
- site inclusion threshold: `slip > 0.0 m` on the mapped closest patch
- rupture bound algorithm: largest contiguous along-strike segment above threshold, bounds taken from the surface row of the fault mesh

The exporter and verifier only support the direct draft-compliant workflow. They do not substitute proxy hypocentres or hidden alternate paths. If a future case lacks an in-window dynamic snapshot for a File 1 event, the exporter will stop with an error.

## Verify Export

```bash
python /Users/dliu/scratch/alpine/benchmark/utils/verify.eqquasi.geometry.conversion.py \
  --case-dir /Users/dliu/scratch/alpine/results/nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad \
  --export-dir /Users/dliu/scratch/alpine/benchmark/Simulation_results/EQquasi_results \
  --prefix alpine_no_dip_change \
  --slip-rate-threshold 0.1
  # --rupture-bound-slip-threshold defaults to 0.5 (matches process.eqquasi.py)
```

This checks:

- geometry round-trip
- station-to-surface-patch mapping
- site slip export against direct recomputation from the raw `Q*` snapshots
- event summary export against direct recomputation from the raw `Q*` snapshots
- geometry conversion figure with:
  - top panel: original STL in UTM with the three benchmark stations
  - bottom panel: processed EQquasi geometry with transformed stations and selected nearest surface patches

The verifier infers the processed geometry and output prefix from `--case-dir` unless you override them explicitly.

## Compare Bundles

```bash
python /Users/dliu/scratch/alpine/benchmark/utils/benchmark.comparison.analytics.py \
  --simulation-root /Users/dliu/scratch/alpine/benchmark/Simulation_results \
  --outdir /Users/dliu/scratch/alpine/benchmark/analytics
```

The analytics script scans all complete bundles under `Simulation_results/*_results`, writes:

- `analytics/bundle_summary.csv`
- `analytics/pairwise_event_comparison.csv`
- `analytics/pairwise_site_comparison.csv`
- `analytics/event_mw_histogram.png`
- `analytics/recurrence_cdf.png`
- `analytics/site1_slip_cdf.png`
- `analytics/site2_slip_cdf.png`
- `analytics/site3_slip_cdf.png`
- `analytics/rupture_extents.png` — time vs. along-strike northing for all models, 3000-yr window, one varying-dip bundle per source
- `analytics/long_term_slip_rate.png` — slip rate vs. along-strike northing for all bundles; EQquasi uses actual `fault.r.nc` slip profiles (via `plot.accumulated.slip.eqquasi` functions), all other codes use `0.5 × max_slip_m` distributed uniformly over the rupture extent

All figures are publication-style (8 × 6 in, 300 dpi, bold fonts/axes/lines). Colors are consistent across all figures via a shared per-bundle palette.

Pairwise CSV rows are only produced when two bundles share the same benchmark scenario name, for example `alpine_varying_dip`.

## Root Runner

The shared shell workflow is:

- [process.eqquasi.sh](/Users/dliu/scratch/alpine/process.eqquasi.sh): rebuild `EQquasi_results`
- [process.pyrsqsim.sh](/Users/dliu/scratch/alpine/process.pyrsqsim.sh): rebuild `PyRSQSim_results`
- [run.comparison.sh](/Users/dliu/scratch/alpine/run.comparison.sh): run both processors, then rebuild `benchmark/analytics`

[run.comparison.sh](/Users/dliu/scratch/alpine/run.comparison.sh) executes the full benchmark workflow:

- export `alpine_planar`
- verify `alpine_planar`
- export `alpine_no_dip_change`
- verify `alpine_no_dip_change`
- export `alpine_varying_dip`
- verify `alpine_varying_dip`
- export `PyRSQSim alpine_varying_dip_5km`
- rebuild `benchmark/analytics`

Run it from the project root with:

```bash
bash /Users/dliu/scratch/alpine/run.comparison.sh
```

At the moment this means:

- `EQquasi alpine_planar` does pair directly with `RSQSim alpine_planar`
- `EQquasi alpine_no_dip_change` does not have a direct RSQSim pair because Avi supplied `alpine_planar`, not `alpine_no_dip_change`
- `EQquasi alpine_varying_dip` does pair directly with `RSQSim alpine_varying_dip`

If a case cannot satisfy the draft requirements from the saved outputs, the strict Python converter will stop with an error instead of silently switching to a proxy rule.

## Inventory

Current benchmark bundle inventory:

- `Simulation_results/EQquasi_results`
  - contains the exported `EQquasi` benchmark bundles
  - `alpine_planar_*`
  - `alpine_no_dip_change_*`
  - `alpine_varying_dip_*`
  - `alpine_planar_geometry_conversion_check.png`
  - `alpine_no_dip_change_geometry_conversion_check.png`
  - `alpine_varying_dip_geometry_conversion_check.png`
- `Simulation_results/RSQSim_results`
  - `alpine_planar_*`
  - `alpine_varying_dip_*`
- `Simulation_results/HBI_results`
  - STL geometry files and a PPTX
  - `alpine_varying_dip_*` bundle (variable backslip run)
  - note: `m0_Nm`, `area_m2`, and `dt_s` columns are `nan` in the event file — HBI does not provide them
- `Simulation_results/PyRSQSim_results`
  - rebuilt from the sibling `PyRSQSim.dev/examples/alpine` workflow
  - current shared shell wrapper exports `alpine_varying_dip_5km_*` using the explicit `onset_centroid` hypocentre policy
- `Simulation_results/MCQsim_results`
  - `alpine_planar_*`
  - `alpine_varying_dip_*`
- `Simulation_results/Tandem_results`
  - empty

## Third-Party Bundle Format Notes

MCQsim and HBI deliver results in formats that differ from the benchmark standard and required normalization on import. These fixes are applied once when unzipping and are not re-applied by any automated script.

### MCQsim

Source zip: `MCQsim_results-*.zip`

- **Whitespace-padded headers**: MCQsim uses fixed-width CSV formatting. All column names and values had leading/trailing whitespace stripped.
- **Site file column casing**: MCQsim site files use `mw`, `slip_m`, `rake_deg`. Renamed to `Mw`, `Slip_m`, `Rake_deg` to match the benchmark standard.
- Event file columns are otherwise identical to the benchmark format.

### HBI

Source zip: `HBI_results-*.zip`

- **File naming**: HBI delivers results as `alpine_catalogue_summary_full_variable_backslip.csv` and `alpine_catalogue_summary_site[1-3]_variable_backslip.csv`. Renamed to `alpine_varying_dip_event_info.csv` and `alpine_varying_dip_site[1-3]_event_info.csv`.
- **Missing event columns**: HBI does not provide `m0_Nm`, `area_m2`, or `dt_s`. These columns are present in the renamed file with value `nan`.
- **Site file column casing**: HBI site files use `mw`, `slip_m`, `Raked_deg`. Renamed to `Mw`, `Slip_m`, `Rake_deg`.

If HBI delivers updated results in the future, apply the same renaming and column-padding steps before running analytics.
