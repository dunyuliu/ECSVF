# Release Note v0.0.2

Date: 2026-05-11

## Summary of Scope

Scientific audit pass over the benchmark toolchain. Key additions: HBI multi-model
automated conversion (`process.hbi.py`), publication-quality figure suite overhaul,
project rules definition, and a full consistency audit across scripts, constants, and
documentation.

## Files Added

| File | Purpose |
|------|---------|
| `benchmark/utils/process.hbi.py` | Automated HBI raw-CSV → benchmark bundle converter (6 models) |
| `AUDIT.md` | Scientific audit findings (units, constants, method, reproducibility) |
| `PROJECT_RULES.md` | Formal project rules used as acceptance checks on every release |
| `docs/release_note_v0.0.1.md` | v0.0.1 archived from repo root |

## Files Modified

| File | Change |
|------|--------|
| `benchmark/utils/benchmark.comparison.analytics.py` | Full figure suite overhaul (see below) |
| `benchmark/utils/plot.accumulated.slip.eqquasi.py` | Removed "coseismic" from titles/labels; correct term is geologic slip rate |
| `benchmark/utils/process.pyrsqsim.py` | Fixed `SECONDS_PER_YEAR` constant; relaxed time-validation tolerance for legacy catalog |
| `process.eqquasi.sh` | Added `--site-slip-threshold 0.5` to all three EQquasi export calls |
| `process.pyrsqsim.sh` | Added `--site-slip-threshold 0.5` and `--hypocenter-policy onset_centroid` |
| `run.comparison.sh` | Added `process.hbi.py` step before analytics |
| `benchmark/README.md` | Updated figure names, documented HBI automation, fixed absolute paths → `${PROJECT_ROOT}/` |
| `.gitignore` | Added `package-lock.json` |

## Figure Suite Overhaul (`benchmark.comparison.analytics.py`)

- Replaced empirical CDF plots → relative-frequency histograms (bars sum to 1)
- Added `plot_mfd` — magnitude-frequency distribution (raw event counts, log y-axis)
- Added `plot_gr_law` — cumulative N(Mw ≥ M) Gutenberg-Richter plot (log y-axis)
- Added `plot_accumulated_slip_along_strike` — accumulated slip vs. along-strike northing, first 3000 yr
- Added `plot_long_term_slip_rate` — long-term slip rate vs. along-strike northing
- Added `plot_rupture_extents` — time vs. along-strike northing, 3000-yr window
- Renamed all output PNGs: `mw_histogram.png`, `mw_frequency.png`, `gutenberg_richter.png`,
  `recurrence_histogram.png`, `slip_accumulated_along_strike.png`, `site{1,2,3}_slip_histogram.png`,
  `rupture_extents.png`, `long_term_slip_rate.png`
- Legend overflow fix: legends placed outside axes (`bbox_to_anchor=(1.02,1)`) with `bbox_inches="tight"`
- Short bundle labels via `shorten_label()`: EQquasi→EQ, MCQsim→MCQ, RSQSim→RSQ, PyRSQSim→pRSQ, HBI→HBI

## HBI Multi-Model Support

`process.hbi.py` replaces the manual HBI normalization workflow:
- Reads `Alpine_jobid_readme.txt` to map 6 job IDs to scenario names
- Converts `Alpine_catalogue_summary_full_*` → `alpine_hbi_{jobid}_{scenario}_event_info.csv`
- Converts `Alpine_catalogue_summary_station{n}_*` → `alpine_hbi_{jobid}_{scenario}_site{n}_event_info.csv`
- Handles empty HBI columns (no `m0_Nm`, `area_m2`, `dt_s`) → `nan`
- `run.comparison.sh` now calls `process.hbi.py` automatically

## Audit Findings and Fixes

| # | Finding | Fix Applied |
|---|---------|-------------|
| 1 | `SECONDS_PER_YEAR = 3.15e7` in `process.pyrsqsim.py` (0.18% off from SI year) | Changed to `365.25 * 24.0 * 3600.0`; validation relaxed to `rel_tol=5e-3` for existing catalog |
| 2 | Site slip threshold not applied to PyRSQSim export | Added `--site-slip-threshold 0.5` to `process.pyrsqsim.sh` |
| 3 | EQquasi site files included near-zero slip events (distant site, fault barely grazed) | Added `--site-slip-threshold 0.5` to all three EQquasi cases in `process.eqquasi.sh` |
| 4 | "Coseismic slip rate" label incorrect — should be geologic/long-term slip rate | Removed "coseismic" from all titles, y-labels, and docstrings |
| 5 | Absolute `/Users/dliu/` paths in `benchmark/README.md` | Replaced with `${PROJECT_ROOT}/` placeholder throughout |
| 6 | HBI conversion was manual (undocumented normalization steps) | `process.hbi.py` automates the full conversion |
| 7 | `mw_histogram` used `density=True` (bin-width normalized, confusing for Mw bins) | Fixed to relative frequency = `counts / counts.sum()` |
| 8 | Legend boxes overflowed figure area | Moved outside axes; `bbox_inches="tight"` on all savefig calls |
| 9 | No formal project rules for CI-style acceptance | `PROJECT_RULES.md` created with 5 rules |

## PROJECT_RULES.md Audit Results

| Rule | Check | Status |
|------|-------|--------|
| Rule 3 — Bundle completeness | All `*_event_info.csv` have matching site1/2/3 files | PASS |
| Rule 4 — No stale duplicates | No two bundles cover the same model run in the same `*_results/` directory | PASS |
| Rule 5 — No hardcoded user paths | No `/Users/<username>/` in any `.py` or `.sh` script | PASS |
| Rule 1 — Verify checks all PASS | Requires live EQquasi pipeline run | Pending (see open items) |
| Rule 2 — Pipeline runs without error | Requires live pipeline run | Pending (see open items) |

## Remaining Open Items

- **Rules 1 & 2 live verification**: Run `bash run.comparison.sh` and
  `python benchmark/utils/verify.eqquasi.geometry.conversion.py` on each EQquasi case to
  confirm all five `[PASS]` checks pass. Cannot be automated at release time without the full
  simulation runtime.
- **PyRSQSim bundle absent**: Local Alpine example catalog has only 1 event (Mw 5.06). Full
  simulation is on the remote Knox server. Re-sync with `bash utig.to.local.sh` and re-run
  `bash process.pyrsqsim.sh` when available.
- **Tandem_results empty**: No Tandem results received yet.
- **`alpine_varying_dip` EQquasi case**: Still running on UTIG Knox. Sync periodically.
