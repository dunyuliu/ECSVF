# AUDIT.md — Alpine Fault EQquasi / Benchmark Project

Date: 2026-05-14
Scope: pre-release audit covering the cross-code benchmark layer (`benchmark/`) and root-level orchestration scripts. The `results/<case>/` archive is in scope only for inventory/stale checks per AGENTS.md (do not reorganize cases).

---

## Executive Summary

Top-3 wins (observed):
1. Benchmark pipeline is internally strict and explicit: `process.eqquasi.py`, `process.pyrsqsim.py`, `verify.eqquasi.geometry.conversion.py`, and `process.hbi.py` all raise on missing inputs rather than silently skipping — matches AGENTS.md "no fallback behavior" rule.
2. `SECONDS_PER_YEAR` is consistently `365.25 * 24.0 * 3600.0` in all three exporters (`process.eqquasi.py:16`, `process.pyrsqsim.py:14`, `plot.accumulated.slip.eqquasi.py:56`). Prior 0.18% drift fixed in v0.0.2.
3. Verifier (`verify.eqquasi.geometry.conversion.py`) recomputes event/site rows from raw `Q*` outputs and cross-checks against the exported CSV with a 1e-6 tolerance — a real round-trip check, not a checksum.

Top-3 risks:
1. **Magic constants duplicated across files (Critical).** `DEPTH_PCT = 88.0`, `GRID_SPACING_M = 2000.0`, and `SLIP_THRESHOLD_M = 0.5` are defined in `plot.accumulated.slip.eqquasi.py:55-57,105`; the analytics script imports `find_q_dirs`/`load_slip_profile`/`cumulative_end_times_yr` from it and reuses these implicitly. There is no asserted link to the actual fault.r.nc grid resolution or the surface row used by `process.eqquasi.py`. If a case is run at a different mesh resolution or different depth convention, the analytics figures will silently misalign.
2. **HBI bundle Mw inconsistency vs. EQquasi/PyRSQSim (Critical).** EQquasi computes Mw from `m0_Nm` via `(2/3)*log10(M0*1e7)-10.7` (`process.eqquasi.py:251-254`). HBI ingests its own `magnitude_Mw` column and leaves `m0_Nm = nan` (`process.hbi.py:111-115`). Pairwise comparisons (`benchmark.comparison.analytics.py`) consume `mw` directly, so the cross-code Mw distribution mixes EQquasi's NM-based Mw with HBI's source-defined Mw without explicit reconciliation of the Mw=⅔log10(M0)−c convention. Suspected discrepancy if HBI uses a different c.
3. **Rupture extent / slip-rate plots for non-EQquasi codes use a uniform-slip proxy.** `benchmark.comparison.analytics.py:643-667` distributes `0.5 * max_slip_m` uniformly between `sbound` and `nbound` for RSQSim, MCQsim, HBI, PyRSQSim. This is documented in `benchmark/README.md:202-204` but produces apples-to-oranges comparisons on `slip_accumulated_along_strike.png` and `long_term_slip_rate.png`. The figures will under-estimate non-EQquasi slip rates in proportion to how peaked their actual rupture is.

---

## 1. Goal & Implementation

**Project in one sentence:** Cross-simulator benchmark layer that exports per-code earthquake catalogs (EQquasi, RSQSim, PyRSQSim, MCQsim, HBI) into a shared Avi-style CSV bundle (`<prefix>_event_info.csv` + three `_siteN_event_info.csv`) and produces comparison statistics and figures for Alpine Fault scenarios (planar, no-dip-change, varying-dip).

**End-to-end trace (orchestrated by `run.comparison.sh:12-18`):**
1. `process.eqquasi.sh:14` clears `EQquasi_results/`.
2. For each of three EQquasi cases (`process.eqquasi.sh:16,32,45`) → `process.eqquasi.py:469` → `build_catalog_rows` (line 388) → per-`Q*` `extract_cycle_details` (line 305) → event/site CSV (line 489-520).
3. After each export, `verify.eqquasi.geometry.conversion.py:327` runs 5 checks (round_trip, station_patch_mapping, site_export_consistency, event_summary_consistency, geometry_visualization).
4. `process.pyrsqsim.sh:15` → `process.pyrsqsim.py:208` reads STL + `rsfsim_history.npz` + `earthquakes.csv` → emits `alpine_varying_dip_5km_*` bundle.
5. `process.hbi.py:129` reads `Alpine_jobid_readme.txt` → emits 6 bundles (jobids 811033, 81102019, 81103019, 821033, 82102019, 82103019).
6. `benchmark.comparison.analytics.py:846` discovers bundles (line 198), writes 3 CSVs + 10 PNGs to `benchmark/analytics/`.

**README vs implementation:**
- `benchmark/README.md:218-227` lists run.comparison.sh sub-steps but omits `process.hbi.py` — the actual script (`run.comparison.sh:14`) runs it. **Doc drift.**
- `benchmark/README.md:262-264` says PyRSQSim is "rebuilt from sibling `PyRSQSim.dev/examples/alpine`". The actual `process.pyrsqsim.sh:8-10` points at `../PyRSQSim.dev/work/alpine_vary_dip`, not `examples/alpine`. **Doc drift.** The Python default (`process.pyrsqsim.py:16`) also disagrees: `PyRSQSim.dev/examples/alpine/results/variable_dip_1km_5km`.
- `README.md:124-126` lists `results/sync_models.py` as a TACC/LS6 sync utility; that exists, but `utig.to.local.sh` is the only sync covered for UTIG/Knox — the README mentions it under "varying-dip" but `sync_models.py` is not the relevant tool for the active varying-dip case.
- `benchmark/README.md:282-289` documents MCQsim/HBI normalization but does **not** say `MCQsim_results/` is still imported manually (no script processes the source zip). MCQsim normalization step is undocumented in code.

---

## 2. Inventory & Stale Items

Top-level (observed):
| Path | One-liner | Recommendation |
|---|---|---|
| `benchmark/` | active cross-code benchmark layer | keep |
| `geometry/` | STL → EQquasi fault-grid prep | keep |
| `results/` | per-case archive, gitignored except listed helpers | keep |
| `archive/post_utility_dev_varyDip_legacy/` | legacy scratch (AGENTS.md acknowledges) | keep (archive) |
| `archive/Earthquake Simulators Intercomparison file format.docx` | older format draft; PDF v0.1 is current | archive (already in archive/) |
| `archive/release_notes/release_note_v0.0.1.md` | duplicate of `docs/release_note_v0.0.1.md` | **delete one** (see below) |
| `docs/release_note_v0.0.1.md` | archived release note per project rule | keep |
| `release_note_v0.0.2.md` | current release note at root | keep |
| `paper/slides` | unspecified | keep (untouched) |
| `__pycache__/` (root and benchmark/utils and every results/*) | bytecode cruft | delete (covered by `.gitignore:6` but committed in tree) |
| `.DS_Store` (10 instances incl. tracked benchmark/, archive/, results/) | macOS metadata | delete (already in `.gitignore:4-5`) |

Stale / duplicated bundles in `benchmark/Simulation_results/` (Rule 4 of PROJECT_RULES.md):
- **No duplicate bundles observed.** HBI emits `alpine_hbi_<jobid>_<scenario>_*` (process.hbi.py:149), avoiding collision with `alpine_planar`/`alpine_varying_dip` from RSQSim/MCQsim.

Stale case directories under `results/`:
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad` and `…slowInitialLoad.casedir.only` — `benchmark/README.md:95` says "no longer the shared benchmark source." Recommendation: **archive** (do not delete; modified within retention window suspected, and AGENTS.md forbids casual reorg).
- `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad.zip` — zip of the active case. Recommendation: keep if it pre-dates the live sync; otherwise archive.

Other:
- `benchmark/__pycache__` / `benchmark/utils/__pycache__` / root `__pycache__` — observed but `.gitignore` covers them; remove from working tree if committed.
- `benchmark/Simulation_results/Tandem_results/` is empty (`benchmark/README.md:268-269` notes this). Keep as placeholder.
- `benchmark/Simulation_results/PyRSQSim_results/` is **empty in the working tree** (sibling repo not present locally). Per PROJECT_RULES.md Rule 2 this is exempt.

---

## 3. Reproducibility

| Check | Status | Detail |
|---|---|---|
| Setup command pinned? | **Missing** | `README.md:147-156` lists `numpy`/`matplotlib`/`xarray`/`netCDF4` with no version pins, no `requirements.txt`, no `pyproject.toml` |
| Sample inputs committed? | Partial | `.gitignore:23-34` excludes `results/*` and `benchmark/Simulation_results/` — sample bundles are local-only, pipeline cannot run on a clean clone without external data |
| Hardcoded paths in pipeline scripts? | **Clean** for benchmark/ utils and root .sh (`PROJECT_RULES.md:36-38` satisfied) |
| Hardcoded paths elsewhere | **Violations in case dirs** — `/Users/dliu/scratch/eqquasi.hbi/bp5.qdc.2000.dipNN` literals in 9 copies of `results/*/plotPeakSliprateTime.py` (lines 40-42 in each) |
| Random seeds | n/a — exporters are deterministic; no seed needed |
| Pinned deps | None |
| Multi-process file locking | n/a (single-process Python exporters); shell wrappers serialize |
| Atomic writes | No — `write_csv` (`process.eqquasi.py:380`, `process.hbi.py:120`) writes in place; if interrupted, partial CSV remains. `run.comparison.sh:10-13` mitigates by `rm -rf` of output dirs at start. |
| Cleanup on Ctrl-C | None explicit |

Suspected issue: `MPLCONFIGDIR` exported in all three .sh wrappers (`process.eqquasi.sh:6`, `process.pyrsqsim.sh:6`, `run.comparison.sh:6`) but the directory is gitignored (`.gitignore:8`) and not pre-created by any script. On a fresh checkout matplotlib will create it lazily — fine — but if `benchmark/.mplconfig` exists from a prior user with wrong perms, ambiguous failure mode. Suggest `mkdir -p` in the .sh wrappers.

---

## 4. Physics & Numerics

**Units walk (observed):**
| Variable | Units | Source |
|---|---|---|
| `local_x_km`, `local_y_km`, `local_z_km` | km | `process.eqquasi.py:129-131` (npz `x`,`y`,`z`) |
| `global_x_m`, `global_y_m`, `global_z_m` | m | line 160-162 (`*1000` from km) |
| `area_m2` | m² | line 136 (`(res_km*1000)**2`) |
| `slips`, `slipd`, `end_slip_mag` | m | xarray dataset, hypot in line 330 |
| `shear_modulus_pa` | Pa | line 241 (`rho*Vs²`, kg/m³ × (m/s)²) |
| `scalar_moment_nm` | N·m | line 339 (`mu * area * slip`) |
| `compute_mw` | dimensionless | line 254 — converts N·m → dyne·cm via `*1e7` then `(2/3)log10 - 10.7` (Kanamori) |
| `t0_s` | s | line 322 |
| `t0_year` | yr | line 436 (s / 365.25·86400) |
| `sbound_NZTM_m` / `nbound_NZTM_m` | m (NZTM northing) | line 374-375 |
| PyRSQSim `time_years` | yr | catalog input — see catalog cross-check at line 251 (5e-3 rel tol for legacy 3.15e7 constant) |

No implicit Pa↔MPa or °C↔K issues observed.

**Sign conventions:**
- `compute_mw` (line 252): explicitly returns `nan` for `moment_nm <= 0`. Good.
- `end_slip_mag = hypot(slips, slipd)` (line 330) is non-negative by construction.
- PyRSQSim event_slip clipped to ≥0 (`process.pyrsqsim.py:244`) after a sign sanity check; explicit.
- Rake computed in **0–360°** (process.eqquasi.py:359) but `circular_mean_deg` (analytics line 270-277) treats it as a circular variable, so the convention is internally consistent. Suspected drift with MCQsim/HBI if they report rake in −180..180° (no observed normalization step in `process.hbi.py`).

**Magic numbers / hardcoded literals (observed):**
| Value | Where | Concern |
|---|---|---|
| `EVENT_WINDOW_SIZE = 20` (steps) | `process.eqquasi.py:17` | undocumented; used as 20-step persistence test for event start/end (line 207-216). No comment on why 20. |
| `slip-rate-threshold default 0.1` (m/s) | `process.eqquasi.py:65` | also re-set to 0.1 in `process.eqquasi.sh`; consistent |
| `rupture-bound-slip-threshold default 0.5` (m) | `process.eqquasi.py:66`, mirrored in verifier line 92, `plot.accumulated.slip.eqquasi.py:105` | duplicated 3x — see Section 5 |
| `site-slip-threshold default 0.0` in .py but 0.5 in .sh | `process.eqquasi.py:68` vs `process.eqquasi.sh:22,30,36,43,49,56` | Default disagrees with shipped pipeline. If user runs the .py directly per README example, they get a different result than the canonical shell wrapper. **Medium risk.** |
| `min-mw default 7.0` | `process.eqquasi.py:64`, `process.pyrsqsim.py:38` | consistent; matches draft spec (`benchmark/README.md:12`) |
| `DEFAULT_RAKE_DEG = 160.0` | `process.pyrsqsim.py:18` | PyRSQSim does not store rake, all events get 160° as a fixed value. **Observed.** This will artificially flatten any cross-code rake comparison involving PyRSQSim. |
| `DEFAULT_SHEAR_MODULUS_PA = 30.0e9` | `process.pyrsqsim.py:19` | Hardcoded 30 GPa for PyRSQSim. EQquasi reads `rho*Vs²` from `model.txt` (line 241). Cross-code Mw comparisons assume both use the same μ. |
| Site coordinates (`x_nztm_m`, `y_nztm_m`) | `process.eqquasi.py:31-35` and `process.pyrsqsim.py:22-26` | **Duplicated literal triples.** If site definitions ever change, two files must be edited. |
| `DEPTH_PCT = 88.0` | `plot.accumulated.slip.eqquasi.py:57` | row index = round(0.88 * (nrows-1)). For a 25-row mesh this is row 21; for a different mesh it picks a different depth. Embedded assumption about mesh layout. |
| `GRID_SPACING_M = 2000.0` | `plot.accumulated.slip.eqquasi.py:55` | Hardcodes 2 km resolution. Geometry npz names contain `.res2km`; if a future case is run at 1 km the x-axis is wrong. |
| `BIN_WIDTH_M = 10e3` (slip rate) | `benchmark.comparison.analytics.py:596,741` | Defined in two functions. |
| `WINDOW_YR = 3000.0` | `benchmark.comparison.analytics.py:530,742` | Defined twice. |
| `1.0e7` (N·m→dyne·cm) | `process.eqquasi.py:254` | Kanamori scaling; correct value, but not symbolic |
| `-10.7` | `process.eqquasi.py:254` | Kanamori constant; correct |
| `CATALOG_SLIP_EPS = 1.0e-6` | `process.pyrsqsim.py:20` | reasonable for "slipped" detection |

**Physical bounds / NaN-Inf risk:**
- `compute_mw` guarded for `M0 <= 0` (line 252-253).
- `np.hypot(slips, slipd)` cannot produce NaN unless inputs are NaN — no observed guard against NetCDF NaNs.
- `bound_cols` empty case handled (line 374-375 returns nan; line 430-431 raises if mw>min_mw triggers). Good.
- `circular_mean_deg`: returns nan if mean vector is zero-length (`benchmark.comparison.analytics.py:275-276`); good.
- `empirical_cdf_distance`: explicit zero-size guard (line 263); good.
- `surface_row_index` (`process.eqquasi.py:286-291`): raises if surface row is not unique. Good.
- `compute_mw` for negative slip would still be fine because `end_slip_mag` is non-negative; suspected: if a Q-folder has an all-zero slip frame, `np.max(end_slip_mag[rupture_mask])` would still be fine because mask is empty → returns None at line 332-333. Good.
- `plot_long_term_slip_rate` divides by `time_span_yr`, guarded `<=0` (line 654-655). Good.
- `time_yr is None or t_yr > WINDOW_YR` filter in `plot_accumulated_slip_along_strike` line 760 — if no profile is loaded `accumulated` stays None and the bundle is silently skipped. Suspected silent failure if a case has all empty `global.dat`.

**Conservation:**
- Scalar moment `M0 = μ · Σ A·D` (line 339) is by definition; no further conservation invariants tested.
- No test that Σ(site slip) is bounded by max_slip_m × area or by long-term plate rate.

**Constants provenance:**
- Kanamori (1977) Mw = (2/3) log10 M0 − 10.7 (M0 in dyne·cm) — standard, no inline citation.
- SECONDS_PER_YEAR = 365.25·86400 — Julian year convention; consistent across files now.
- μ = ρ·Vs² for EQquasi (from model.txt); μ=30 GPa hardcoded for PyRSQSim — provenance unclear.

---

## 5. Implementation Consistency

| # | Issue | Files | Severity |
|---|---|---|---|
| 5.1 | Site coordinate triples duplicated | `process.eqquasi.py:31-35` ↔ `process.pyrsqsim.py:22-26` | Medium |
| 5.2 | `rupture-bound-slip-threshold = 0.5` defaulted in three places | `process.eqquasi.py:66`, `verify.eqquasi.geometry.conversion.py:92`, `plot.accumulated.slip.eqquasi.py:105` | Medium |
| 5.3 | `--site-slip-threshold` default mismatch: 0.0 in .py vs 0.5 in .sh | `process.eqquasi.py:68` vs `process.eqquasi.sh:22,30,36,43,49,56`; `process.pyrsqsim.py:39` vs `process.pyrsqsim.sh:21` | Medium |
| 5.4 | EQquasi uses real on-fault slip; non-EQquasi uses `0.5×max_slip` uniform proxy | `benchmark.comparison.analytics.py:643-667, 781-803` | High (semantic) |
| 5.5 | μ provenance differs: EQquasi from model.txt; PyRSQSim hardcoded 30 GPa | `process.eqquasi.py:241` vs `process.pyrsqsim.py:19` | Medium |
| 5.6 | Mw column source differs by code; HBI bundle has nan m0_Nm | `process.hbi.py:111-115` and resulting CSVs (column `m0_Nm` all `nan`) | High |
| 5.7 | `PyRSQSim_results` dir documented as default for `examples/alpine` (`benchmark/README.md:42-46`), but shell wrapper uses `work/alpine_vary_dip` | `process.pyrsqsim.sh:8` | Low (doc drift) |
| 5.8 | `BIN_WIDTH_M`, `WINDOW_YR` duplicated across two analytics functions | `benchmark.comparison.analytics.py:530,596,741,742` | Low |
| 5.9 | Verifier `--site-slip-threshold` default 0.0 (line 94) vs `process.eqquasi.sh` 0.5 → verifier will recompute against 0.0 but compare to a 0.5-thresholded export → site export consistency check WILL FAIL unless shell passes `--site-slip-threshold 0.5` to verifier too. **Observed: `process.eqquasi.sh:30,43,56` does pass `--site-slip-threshold 0.5` to the verifier.** OK in practice but contract is fragile. | `verify.eqquasi.geometry.conversion.py:94` | Low |

**Configured vs actual:**
- `benchmark/README.md:202-204` says EQquasi accumulated-slip uses `fault.r.nc`; code confirms (`benchmark.comparison.analytics.py:614-624`). Match.
- `benchmark/README.md:218-227` describes `run.comparison.sh` flow but does not list `process.hbi.py`; code (line 14) calls it. **Doc drift.**

**Port equivalence:** N/A (single language).

**Data pipeline integrity:**
- Cycle-time offsets: `build_cycle_offsets` (process.eqquasi.py:225-232) accumulates `times_s[-1]` of each Q*. Observed in alpine_planar bundle: event 0 t0_year ≈ 2.1e-7 (first cycle, start at t≈6.6 s) while event 1 t0_year ≈ 440 (after a 13.9 Gs offset). Assumes each Q* ends at the start of the next; if a Q* is truncated or skipped, offsets are wrong. No observed guard.
- Train/val/test leakage: N/A (not an ML pipeline).
- Time-alignment of multi-source inputs: PyRSQSim verifies catalog `time_years` against history `time` with rel_tol 5e-3 (line 251). Good.

---

## 6. Logging & Error Handling

| Aspect | Observation | File:Line |
|---|---|---|
| Exporters fail loudly on missing inputs | Yes — `FileNotFoundError`/`RuntimeError`/`ValueError` raised explicitly | `process.eqquasi.py:121,173,177,238,329,429,431,465`; `process.pyrsqsim.py:63,67,80,82,89,96,126,141,157,160,167,222,224,228,240,243,249,252,261,288,308` |
| Verifier prints `[PASS]/[FAIL]` and exits non-zero on any fail | Yes | `verify.eqquasi.geometry.conversion.py:371-375` |
| Logging library used | No — bare `print()` | `process.eqquasi.py:512,520`; `process.hbi.py:126`; etc. |
| Git SHA / case ID in logs | No | — |
| Timestep / cycle ID in logs | Indirect (prints output path) | — |
| Per-case shell wrapper continue-on-error? | No: `set -euo pipefail` in all .sh scripts | `process.eqquasi.sh:3`, `process.pyrsqsim.sh:3`, `run.comparison.sh:3` |
| Silent skip risks | **Yes** — `extract_cycle_details` returns `None` if no event window or no rupture (process.eqquasi.py:319, 333). The Q* is silently dropped from the catalog. No log line. If a real event was missed because slip-rate threshold was too high, the user has no indication. | `process.eqquasi.py:319,333` |
| Silent skip risks (analytics) | `plot_long_term_slip_rate` and `plot_accumulated_slip_along_strike` `continue` on empty EQquasi profile or missing case dir (lines 626, 652, 654, 658, 769, 787) with no warning | `benchmark.comparison.analytics.py:626,652,769` |
| Cascade behavior on partial bundle | `discover_bundles` raises `FileNotFoundError` if any of site1/2/3 is missing (line 222). Good. | `benchmark.comparison.analytics.py:220-222` |

---

## 7. Performance & Scaling

Observed: dataset sizes are small (event catalogs of 10s to low 100s of rows per bundle; ~25 surface columns). Performance is not a bottleneck.

Specific patterns:
- `extract_cycle_details` opens every `fault.NNNN.nc` snapshot inside `[start_step, end_step]` looking for the first slipping frame (line 342-352). For Q* with many snapshots and a late onset, this is O(N) NetCDF opens per event. Could short-circuit by binary search on snapshot file index, but in practice events are rare and N is small.
- `n_cumulative = [np.sum(mw >= m) for m in mw_unique]` in `plot_gr_law` (line 722) is O(n²). Trivial for current sizes.
- `cycle_dirs` is called multiple times in `build_catalog_rows` (`process.eqquasi.py:401-402`). Filesystem cost, but small.

Memory: all arrays are surface-strike or per-event scalars. No O(N²) array growth observed.

---

## 8. Top-N Priorities (ranked by impact × ease, <1 day each)

1. **[Critical, easy] Unify the `--site-slip-threshold` default between Python and shell.** Pick one (0.5) and set it as the Python default in `process.eqquasi.py:68`, `process.pyrsqsim.py:39`, and `verify.eqquasi.geometry.conversion.py:94`. Reduces silent contract drift (Finding 5.3).
2. **[Critical, easy] Document Mw provenance per code in `benchmark/README.md` and confirm HBI's Mw convention.** Add an explicit table: code → Mw formula → constants used. If HBI's `magnitude_Mw` does not use Kanamori c=10.7, recompute from rupture area + assumed slip or flag in the analytics legend (Finding 5.6, Risk 2).
3. **[High, easy] Extract shared constants to a single module.** Create `benchmark/utils/constants.py` with `SECONDS_PER_YEAR`, `SLIP_THRESHOLD_M`, `SITES`, `KANAMORI_C`, `WINDOW_YR`, `BIN_WIDTH_M`, and import from all utilities. Removes Findings 5.1, 5.2, 5.8.
4. **[High, easy] Fix README ↔ pipeline drift.** Update `benchmark/README.md:42-46,218-227` to list `process.hbi.py` and the actual `PyRSQSim.dev/work/alpine_vary_dip` path. Update `process.pyrsqsim.py:16` default to match the .sh wrapper, or vice versa (Findings 1, 5.7).
5. **[High, medium] Add a `requirements.txt` (or `pyproject.toml`) pinning `numpy`, `xarray`, `netCDF4`, `matplotlib` to current working versions.** Fixes Section 3 reproducibility gap.
6. **[Medium, easy] Replace silent `continue`/`return None` with a warning log line in `extract_cycle_details` and the EQquasi branches of the two slip-rate plotters.** A missed event should not be invisible (Section 6 silent-skip risks).
7. **[Medium, easy] Document `DEPTH_PCT = 88.0` and `GRID_SPACING_M = 2000.0` in `plot.accumulated.slip.eqquasi.py` with an assertion against `fault.r.nc` shape and the npz `resolution` metadata.** Fail loudly if a future case has a different grid (Finding 4 magic numbers; Risk 1).
8. **[Medium, easy] Reconcile rake units across codes.** Add an explicit normalization in `process.hbi.py` and a sanity-check `(0 <= rake < 360)` assertion in the analytics loader, so that MCQsim/HBI rakes match EQquasi's 0–360° convention.
9. **[Medium, medium] Replace the uniform-slip proxy with a more honest plotting choice for non-EQquasi codes.** Either (a) plot only EQquasi on accumulated-slip / slip-rate panels and provide a separate panel for the catalog-derived rate proxy, or (b) annotate the figure caption inside the script (already in README but not on the figure). Risk 3.
10. **[Low, easy] Remove tracked `__pycache__/` and `.DS_Store` from the working tree.** Already in `.gitignore`; just `git rm -r --cached`.
11. **[Low, easy] Decide between `archive/release_notes/release_note_v0.0.1.md` and `docs/release_note_v0.0.1.md` — keep one.** PROJECT_RULES.md (Claude CLAUDE.md release workflow) says `docs/`; remove the archive duplicate.
12. **[Advisory] Add a `mkdir -p "${MPLCONFIGDIR}"` in the three .sh wrappers** to make the matplotlib config path explicit on fresh checkouts.

---

## Distinguishing Observed vs Suspected

- **Observed** (direct inspection of code/CSV/file tree): Findings 1, 5.1–5.5, 5.7–5.9, 6 silent-skip locations, Section 2 inventory, Section 4 magic-number table, Section 3 pinning gap, t0_year jump example.
- **Suspected** (would require running the pipeline to confirm): Finding 5.6 HBI Mw equivalence (depends on HBI's internal constant); HBI rake unit (Finding 8 in priorities); behavior of `cycle_offsets` when a Q* is truncated; whether `MCQsim_results/` has whitespace residues after the manual import documented in benchmark/README.md:278-281.

Sign-off rests with the human reviewer. No code edits applied by this audit.
