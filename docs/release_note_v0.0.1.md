# Release Note v0.0.1

Date: 2026-04-21

## Release Scope

Initial public release of the New Zealand Alpine Fault EQquasi cross-code benchmark repository (ECSVF). This release establishes the tracked repository structure, geometry assets, benchmark utilities, and workflow scripts. All large simulation archives, benchmark result bundles, and analytics outputs remain local-only and are excluded from the tracked repository via `.gitignore`.

## Repository Structure (Tracked)

```
.
├── .gitignore
├── README.md
├── release_note_v0.0.1.md
├── process.eqquasi.sh          # rebuild EQquasi benchmark bundles
├── process.pyrsqsim.sh         # rebuild PyRSQSim benchmark bundle
├── run.comparison.sh           # rebuild all bundles + analytics
├── utig.to.local.sh            # sync varying-dip case from UTIG Knox
├── benchmark/
│   ├── README.md
│   └── utils/
│       ├── process.eqquasi.py                    # export EQquasi Q-folders → benchmark CSVs
│       ├── process.pyrsqsim.py                   # export PyRSQSim → benchmark CSVs
│       ├── verify.eqquasi.geometry.conversion.py # verify EQquasi export geometry
│       ├── benchmark.comparison.analytics.py     # cross-code comparison figures
│       └── plot.accumulated.slip.eqquasi.py      # EQquasi slip profiles and slip rate
├── geometry/
│   ├── README.md
│   ├── processing_scenarios_summary.txt
│   ├── transfer.geometry.to.ls6.sh
│   ├── earthquake_model_case{1,2,3}.png
│   ├── utils/
│   │   ├── stl_to_eqquasi_fault_grid.py
│   │   ├── fault_trace_analysis.py
│   │   ├── earthquake_model_3d.py
│   │   └── stl_transform.py
│   ├── raw_geometry/
│   │   ├── geoms_202508/          # pre-workshop STL meshes
│   │   └── post_workshop_geoms/   # final workshop STL meshes (active)
│   └── processed_geometry/        # processed grids (.npz, .csv), geometry reports
└── results/
    ├── README.md
    └── <shared post-processing scripts>  # plotAccumulated, plotOnFaultVars, etc.
```

## Local-Only (Not Tracked)

| Path | Reason |
|------|--------|
| `results/<case>/` | Full simulation archives (NetCDF, Q* folders) |
| `benchmark/Simulation_results/` | Result CSVs and STL files from all codes |
| `benchmark/analytics/` | Generated comparison figures and summary tables |
| `archive/` | Legacy scripts, old release notes, proxy bundles |
| `CLAUDE.md`, `AGENTS.md` | Local AI instruction files |

## Audit Findings and Fixes Applied

- `geometry/src/` renamed to `geometry/utils/` for consistency with `benchmark/utils/`
- Three throwaway geometry scripts deleted: `test_stl.py`, `view.fault.geometry.py`, `visualize.stl.py`
- `archive/`, `benchmark/Simulation_results/`, `benchmark/analytics/` added to `.gitignore` and removed from index
- `AGENTS.md` and `CLAUDE.md` added to `.gitignore`
- `benchmark/docs/` and `benchmark/archive/` moved to local `archive/` (gitignored)
- All READMEs updated to reflect `geometry/utils/` rename, removal of archive references, and local-only status of `Simulation_results/` and `analytics/`
- All prior git commits stripped of `Co-Authored-By: Claude` trailer
- Git history squashed to single init commit

## Active Simulation Cases (Local)

| Case directory | Benchmark label |
|---------------|-----------------|
| `nz.bp5.qdc.dip50.2000.norm_6mm_yr` | `alpine_planar` |
| `nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad` | `alpine_no_dip_change` |
| `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad` | `alpine_varying_dip` (active, running on UTIG Knox) |

## Benchmark State (Local)

Active codes in the comparison: `EQquasi`, `RSQSim`, `PyRSQSim`, `MCQsim`, `HBI`.

Key EQquasi exporter settings:
- Rupture bound threshold: 0.5 m (overrides draft 0.1 m to suppress noise fragments)
- Rupture bounds computed from surface row only (avoids dip-induced northing bias)
- Largest contiguous along-strike segment algorithm for rupture extent

## Known Limitations / Open Items

- `alpine_varying_dip` EQquasi case is still running on UTIG Knox. Sync with `bash utig.to.local.sh`, then re-run `bash process.eqquasi.sh` and analytics to pick up new Q folders.
- `PyRSQSim_results` bundle is not current; the sibling PyRSQSim Alpine example has not been updated to match the current mesh.
- `benchmark/Simulation_results/Tandem_results/` is empty — Tandem results not yet received.
- Raw HBI and MCQsim source ZIPs are local-only; format normalization steps are manual (documented in `benchmark/README.md`).
