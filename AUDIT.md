# AUDIT.md — Alpine Fault Earthquake Simulator Benchmark

**Audit date:** 2026-05-11  
**Auditor:** Claude (automated static analysis)  
**Scope:** benchmark pipeline only (`benchmark/utils/`, shell wrappers, associated docs)

---

## Executive Summary

**Top-3 wins**
1. The end-to-end pipeline is clean and self-contained: each processing script is strict-by-default, raises on missing data rather than substituting silently, and column names are consistent from raw → bundle → analytics.
2. The `mw` formula in `process.eqquasi.py` (Hanks & Kanamori 1979: `(2/3)*log10(M0_Nm×1e7)−10.7`) is correct; `empirical_cdf_distance` correctly implements the KS statistic.
3. `process.hbi.py` column mapping matches every actual raw-CSV header column-for-column (verified against live files).

**Top-3 risks**
1. **SECONDS_PER_YEAR mismatch** — `process.pyrsqsim.py` uses `3.15e7` (0.18% low), all other scripts use `365.25*24*3600`. Over a 10 000-yr simulation this drifts `t0_year` by ~18 yr and silently corrupts the cross-code time comparison.
2. **Legacy HBI bundle `alpine_varying_dip_*` coexists with six new `alpine_hbi_*` bundles** — analytics loads all seven as independent "HBI:varying_dip" entries, inflating HBI representation in every pairwise comparison and all figures.
3. **`benchmark/README.md` filenames are stale** — documents three PNG names (`event_mw_histogram.png`, `recurrence_cdf.png`, `site{1,2,3}_slip_cdf.png`) that do not exist on disk; actual files have different names. Any developer following the README to locate outputs will fail.

---

## Section 1 — Goal & Implementation

**Project summary.** This repository compares earthquake-cycle simulation codes (EQquasi, RSQSim, PyRSQSim, MCQsim, HBI) on the New Zealand Alpine Fault by converting heterogeneous outputs into a common Avi-style CSV bundle format and then generating cross-code comparison figures.

### Main pipeline (traced end-to-end)

```
Raw outputs                 Processor                     Bundle CSVs
─────────────────────────────────────────────────────────────────────
Q*/fault.r.nc + global.dat  process.eqquasi.py            EQquasi_results/*
PyRSQSim NPZ + CSV          process.pyrsqsim.py           PyRSQSim_results/*
HBI raw catalogue CSVs      process.hbi.py                HBI_results/alpine_hbi_*
(MCQsim / RSQSim)           (manual rename + strip)       MCQsim_results/* / RSQSim_results/*
                                    │
                          benchmark.comparison.analytics.py
                                    │
                          analytics/{*.csv, *.png}
```

Shell entry points:
- `process.eqquasi.sh` → runs `process.eqquasi.py` × 3 cases + `verify.eqquasi.geometry.conversion.py` × 3
- `process.pyrsqsim.sh` → runs `process.pyrsqsim.py` once with `--hypocenter-policy onset_centroid`
- `run.comparison.sh` → sources both shells then calls `benchmark.comparison.analytics.py`
- `process.hbi.py` is **not wired into any shell script** — must be run manually.

### README claims vs actual implementation

| README claim | Actual |
|---|---|
| `analytics/event_mw_histogram.png` | `analytics/mw_histogram.png` (MISMATCH) |
| `analytics/recurrence_cdf.png` | `analytics/recurrence_histogram.png` (MISMATCH) |
| `analytics/site{1,2,3}_slip_cdf.png` | `analytics/site{1,2,3}_slip_histogram.png` (MISMATCH) |
| "10 PNGs" (implied by root README) | 10 PNGs generated (CORRECT) |
| `benchmark/README.md` L212: "`run.comparison.sh` executes verify steps" | Verify steps are in `process.eqquasi.sh`, not `run.comparison.sh` (minor) |
| `benchmark/README.md` L29–31: absolute `/Users/dliu/scratch/alpine/` paths in all shell references | Scripts themselves use `${ROOT_DIR}` (portable); only the markdown docs are hardcoded |

### Dead links / missing files referenced in docs

| Location | Reference | Status |
|---|---|---|
| `benchmark/README.md` line 9 | `benchmark/docs/Earthquake Simulators Intercomparison file format_v0.1` | MISSING — no `benchmark/docs/` directory exists |
| `benchmark/README.md` lines 29–31 | Absolute paths `/Users/dliu/scratch/alpine/process.*.sh` | Functional but non-portable (only true on author's machine) |
| `benchmark/README.md` line 46 | `/Users/dliu/scratch/PyRSQSim.dev/examples/alpine` | External repo, not present in this repo — expected but undocumented dependency |

---

## Section 2 — Inventory & Stale Items

### Top-level structure

| Path | Role | Status |
|---|---|---|
| `README.md` | Project overview and workflow | Active |
| `CLAUDE.md` | Project AI-assistant rules | Active |
| `AGENTS.md` | Directory scope notes for AI | Active |
| `process.eqquasi.sh` | Shell wrapper: rebuild EQquasi bundles | Active |
| `process.pyrsqsim.sh` | Shell wrapper: rebuild PyRSQSim bundle | Active |
| `run.comparison.sh` | Full workflow runner | Active |
| `utig.to.local.sh` | HPC sync helper | Active |
| `release_note_v0.0.1.md` | Release notes (should be in `docs/`) | Active (not yet moved) |
| `package-lock.json` | Node.js lockfile | SUSPECT — no `package.json` or JS tooling visible; likely committed by accident |
| `__pycache__/` | Python bytecode cache | Stale artifact; excluded from git via `.gitignore` (good) |
| `benchmark/utils/` | Core Python utilities | Active |
| `benchmark/Simulation_results/` | Bundle data (local-only, not tracked) | Active |
| `benchmark/analytics/` | Generated outputs (local-only, not tracked) | Active |
| `benchmark/Simulation_results/Tandem_results/` | Empty directory | Stale placeholder |
| `benchmark/Simulation_results/HBI_results/alpine_varying_dip_*` | Old manual HBI bundle (Apr 21) | STALE — superseded by `alpine_hbi_81102019_*` generated by `process.hbi.py` on May 11; both are active in analytics simultaneously |

### Flagged items

| File | Finding |
|---|---|
| `package-lock.json` (root) | OBSERVED: no `package.json` exists; this file has no purpose here. |
| `HBI_results/alpine_varying_dip_event_info.csv` + 3 site files | OBSERVED: created 2026-04-21, predates `process.hbi.py`; has `nan` in `m0_Nm`, `area_m2`, `dt_s`; now shadowed by `alpine_hbi_81102019_varying_dip_*` but still loaded by analytics. |
| `benchmark/Simulation_results/Tandem_results/` | OBSERVED: empty directory; `benchmark/README.md` lists it as inventory but no data is present. |

---

## Section 3 — Reproducibility

### Hardcoded absolute paths

| File | Line | Path | Risk |
|---|---|---|---|
| `benchmark/README.md` | 29–31, 46, 100–145 | `/Users/dliu/scratch/alpine/…` and `/Users/dliu/scratch/PyRSQSim.dev/…` | Documentation-only; shell scripts use `${ROOT_DIR}` (portable) |
| `process.eqquasi.py` | 19–20 | `REPO_ROOT`-relative defaults | Resolved at runtime via `Path(__file__).parents[2]` — portable |
| `process.pyrsqsim.py` | 16–17 | `REPO_ROOT.parent / "PyRSQSim.dev/…"` | Portable relative reference, but only works when sibling repo exists at exactly one level up |

OBSERVED: No `~/` expansions anywhere in Python scripts.

### Pinned dependencies

No `requirements.txt`, `pyproject.toml`, `environment.yml`, or `setup.py` exists anywhere in the repository. The root `README.md` lists `numpy`, `matplotlib`, `xarray`, `netCDF4` as minimum requirements but without version pins.

### Sample inputs committed

No sample `Q*` folders, `fault.r.nc`, or simulator CSVs are committed. The `benchmark/Simulation_results/` tree is `.gitignore`-listed as local-only. The `.gitignore` correctly excludes `benchmark/Simulation_results/` and `benchmark/analytics/`.

### Atomic writes

OBSERVED: All `write_csv()` functions open the file with `path.open("w")` and write directly — no temp-file-then-rename pattern. A crash mid-write would leave a truncated CSV with a valid header but partial data, silently accepted on the next analytics run.

---

## Section 4 — Physics & Numerics

### Unit constants traced across files

| Constant | `process.eqquasi.py` | `process.pyrsqsim.py` | `plot.accumulated.slip.eqquasi.py` | `benchmark.comparison.analytics.py` |
|---|---|---|---|---|
| Seconds per year | `365.25 * 24.0 * 3600.0` = 31 557 600 (line 16) | `3.15e7` = 31 500 000 (line 14) | `365.25 * 24 * 3600` = 31 557 600 (line 56) | not defined (reads `t0_year` from CSV) |
| Grid spacing (m) | not used | not used | `GRID_SPACING_M = 2000.0` (line 55) | `BIN_WIDTH_M = 10e3` (lines 596, 739) |
| Depth pct (%) | not used | not used | `DEPTH_PCT = 88.0` (line 57) | not used (calls `plot.accumulated.slip` module) |
| Slip threshold (m) | `--rupture-bound-slip-threshold` default 0.5 (line 66) | `--rupture-bound-slip-threshold` default 0.1 (line 41) | `SLIP_THRESHOLD_M = 0.5` (line 105) | not used directly |
| Bin width (m) | — | — | — | `BIN_WIDTH_M = 10e3` in two functions (lines 596, 739) — consistent |

### SECONDS_PER_YEAR discrepancy — OBSERVED

`process.pyrsqsim.py:14` uses `3.15e7` (31 500 000 s/yr).  
All other scripts use `365.25 * 24 * 3600` = 31 557 600 s/yr.  
**Error: 0.1825%** (−57 600 s/yr). Over a 10 000-yr simulation, `t0_year` for events near the end is ~18 yr early in the PyRSQSim bundle relative to all other codes. This biases recurrence interval comparisons.

Additionally, `process.pyrsqsim.py:250` uses the same constant for a validation check:
```python
if not math.isclose(float(time_s[end_idx] / SECONDS_PER_YEAR), float(catalog_row["time_years"]), rel_tol=0.0, abs_tol=1.0e-9):
```
If `earthquakes.csv` was written with the correct 31 557 600 constant, this check will fail for every event due to the constant mismatch. SUSPECTED: the PyRSQSim catalog was also written with `3.15e7`, so the check happens to pass — but this is fragile.

### Rupture-bound slip threshold inconsistency — OBSERVED

- `process.eqquasi.py` default `--rupture-bound-slip-threshold 0.5` m (line 66)  
- `process.pyrsqsim.py` default `--rupture-bound-slip-threshold 0.1` m (line 41)

The `process.eqquasi.sh` passes `0.5` explicitly (consistent with the Python default). The `process.pyrsqsim.sh` does not pass `--rupture-bound-slip-threshold`, so PyRSQSim uses 0.1 m. The 0.5 m choice is documented in `benchmark/README.md` as suppressing noise, but is applied only to EQquasi, not to PyRSQSim. This inconsistency is silently present.

### Site-slip-threshold asymmetry — OBSERVED

- `process.eqquasi.py` CLI default: `--site-slip-threshold 0.0` (line 68)
- `process.eqquasi.sh` passes `--site-slip-threshold 0.5` m for all three EQquasi cases (lines 22, 36, 49, 56)
- `process.pyrsqsim.py` default: `--site-slip-threshold 0.0` (line 39); `process.pyrsqsim.sh` does not override it
- `process.hbi.py` does no slip filtering at all — it passes through whatever HBI reports
- MCQsim / RSQSim: post-receive manual normalization with no slip filter applied

The 0.5 m threshold is applied to EQquasi site events but **not** to any other code. This reduces EQquasi site-event counts and systematically biases slip-CDF comparisons across codes. The threshold is not documented as intentionally asymmetric anywhere.

### Magic numbers — all hardcoded physics literals

| Value | Location | Meaning | Notes |
|---|---|---|---|
| `365.25` | process.eqquasi.py:16, plot.accumulated.slip.eqquasi.py:56 | days/yr | Consistent |
| `3.15e7` | process.pyrsqsim.py:14 | s/yr (wrong) | Should be 31 557 600 |
| `1.0e7` | process.eqquasi.py:254 | N·m → dyne·cm conversion | Correct (Hanks-Kanamori) |
| `10.7` | process.eqquasi.py:254 | Hanks-Kanamori offset | Standard |
| `2.0 / 3.0` | process.eqquasi.py:254 | Mw scaling | Correct |
| `30.0e9` | process.pyrsqsim.py:19 | default shear modulus (Pa) | Consistent with typical crustal value |
| `160.0` | process.pyrsqsim.py:20 | default rake (deg) | Alpine Fault dextral-reverse; undocumented |
| `1.0e-6` | process.pyrsqsim.py:21 | catalog slip epsilon (m) | No comment |
| `0.5` | process.eqquasi.py:66, plot.accumulated.slip.eqquasi.py:105 | rupture-bound threshold (m) | Consistent between these two files |
| `0.1` | process.pyrsqsim.py:41, process.eqquasi.py:65 (slip-rate) | bound threshold / slip-rate threshold | Different meanings; documented |
| `88.0` | plot.accumulated.slip.eqquasi.py:57 | DEPTH_PCT | Not present in analytics.py; analytics calls the module |
| `2000.0` | plot.accumulated.slip.eqquasi.py:55 | GRID_SPACING_M (m) | Used only for distance axis; correct for 2km resolution |
| `10e3` | benchmark.comparison.analytics.py:596, 739 | BIN_WIDTH_M (m) | Both occurrences consistent |
| `3000.0` | benchmark.comparison.analytics.py:530, 740 | analysis window (yr) | Consistent |
| `0.5` (slip approximation) | benchmark.comparison.analytics.py:661, 799 | 0.5 × max_slip for non-EQquasi | SUSPECTED: undocumented approximation; see below |
| `25` | plot.accumulated.slip.eqquasi.py:172 | hardcoded nrows for subtitle only | OBSERVED: used only in string `depth_row` display, not data; but assumes 25-row grid |

### 0.5 × max_slip approximation — OBSERVED, undocumented

In `benchmark.comparison.analytics.py:661` and `799`:
```python
slip_accum[in_rupture] += 0.5 * slip
```
This approximates mean slip per event as 50% of the reported `max_slip_m`, applied uniformly over the full rupture extent. `benchmark/README.md` line 201 mentions it but calls it a method choice, not an approximation. No source or justification is cited. For a triangular slip distribution the correct factor is ~0.33–0.5; for EQquasi the actual spatial profile is used instead. The approximation is applied to RSQSim, PyRSQSim, HBI, and MCQsim — every non-EQquasi code — and is asymmetric (EQquasi uses real profiles). This is documented in the README but the analytical basis is not justified.

### `compute_moment` in process.pyrsqsim.py — OBSERVED correct

`process.pyrsqsim.py:176–177`:
```python
def compute_moment(area_m2, slip_m, shear_modulus_pa):
    return float(shear_modulus_pa * np.sum(area_m2 * slip_m))
```
Standard scalar seismic moment. Used for output `m0_Nm` but **not used for `mw`** — `mw` is taken directly from `catalog_row["magnitude"]` (line 256). So the computed moment goes into `m0_Nm` but is inconsistent with the reported `mw` unless PyRSQSim's internal moment matches `DEFAULT_SHEAR_MODULUS_PA = 30 GPa`. SUSPECTED: if PyRSQSim uses a different shear modulus internally, `m0_Nm` and `mw` will be physically inconsistent for the same event.

### `empirical_cdf_distance` — OBSERVED correct

`benchmark.comparison.analytics.py:261–267`: uses `np.searchsorted` on pre-sorted arrays, evaluated at the union of both samples. This is the standard KS statistic. Correct.

### `circular_mean_deg` — OBSERVED correct

`benchmark.comparison.analytics.py:270–277`: uses complex exponential representation `np.exp(1j * angles)`, takes mean, returns angle. Handles the wraparound discontinuity correctly. Returns NaN for zero-length resultant (undefined mean direction). Correct.

### `mw` formula — OBSERVED correct

`process.eqquasi.py:254`:
```python
(2.0 / 3.0) * math.log10(moment_nm * 1.0e7) - 10.7
```
Converts N·m to dyne·cm (×10⁷), applies Hanks & Kanamori (1979). Correct.  
Guards against `moment_nm <= 0.0` (line 252), returns NaN. Correct.

### NaN / Inf risks

| Location | Risk | Verdict |
|---|---|---|
| `process.eqquasi.py:254` | `log10(moment_nm)` — guarded by `if moment_nm <= 0.0` | SAFE |
| `process.eqquasi.py:339` | `np.sum(geometry.area_m2[rupture_mask] * end_slip_mag[rupture_mask])` — only reached if `np.any(rupture_mask)` | SAFE |
| `process.pyrsqsim.py:79` | `area = 0.5 * np.linalg.norm(cross)` — raises on degenerate triangle | SAFE (raises) |
| `benchmark.comparison.analytics.py:651` | `time_span_yr = years.max() - years.min()` guarded by `if time_span_yr <= 0` | SAFE |
| `process.pyrsqsim.py:171` | `np.sqrt(...)` — arguments are squared, always ≥0 | SAFE |
| `plot.accumulated.slip.eqquasi.py:221` | `running[-1] / t_arr.max()` — only reached after at least one profile loaded, so `t_arr.max() > 0` | SAFE |
| `benchmark.comparison.analytics.py:633` | `northing_m = np.linspace(min(sbounds_eq), max(nbounds_eq), ...)` — `sbounds_eq` could be empty if all events have NaN bounds | SUSPECTED risk: if `sbounds_eq` is empty, `min([])` raises `ValueError` |

---

## Section 5 — Implementation Consistency

### SECONDS_PER_YEAR across files

| File | Value | Correct? |
|---|---|---|
| `process.eqquasi.py:16` | `365.25 * 24.0 * 3600.0` = 31 557 600 | Yes |
| `process.pyrsqsim.py:14` | `3.15e7` = 31 500 000 | **No — 0.1825% low** |
| `plot.accumulated.slip.eqquasi.py:56` | `365.25 * 24 * 3600` = 31 557 600 | Yes |
| `benchmark.comparison.analytics.py` | Not defined; reads `t0_year` from CSV | N/A |

### Column name mapping in `process.hbi.py` vs actual HBI raw CSV headers

Headers verified against live files on disk:

**Event file** (`Alpine_catalogue_summary_full_*.csv`):
```
event_ID, start_time_sec, start_time_yr, magnitude_M0, magnitude_Mw,
rupture_area_m2, duration_sec, max_slip_m, hypo_x_NZTM, hypo_y_NZTM,
hypo_z_NZTM, s_bound_NZTM, n_bound_NZTM
```
OBSERVED: All 13 source columns in `EVENT_COL_MAP` match the actual headers exactly. No mismatches.

**Site file** (`Alpine_catalogue_summary_station{n}_*.csv`):
```
event_ID, event_time_sec, event_time_yr, magnitude_Mw, event_slip_m, rake
```
OBSERVED: All 5 source columns in `SITE_COL_MAP` match exactly. Note: `benchmark/README.md` "Third-Party Bundle Format Notes" section states HBI uses `"Raked_deg"` (with a 'd') — this is **incorrect documentation**. The actual column is `"rake"`, and `process.hbi.py` correctly maps `"rake"` → `"Rake_deg"`.

### `infer_scenario()` vs actual bundle names

| Bundle ID | Result | Correct? |
|---|---|---|
| `alpine_no_dip_change` | `no_dip_change` | Yes |
| `alpine_planar` | `planar` | Yes |
| `alpine_varying_dip` | `varying_dip` | Yes |
| `alpine_hbi_81102019_varying_dip` | `varying_dip` | Yes |
| `alpine_hbi_81103019_planar` | `planar` | Yes |
| `alpine_hbi_811033_planar` | `planar` | Yes |
| `alpine_hbi_82102019_varying_dip` | `varying_dip` | Yes |
| `alpine_hbi_82103019_planar` | `planar` | Yes |
| `alpine_hbi_821033_planar` | `planar` | Yes |
| `alpine_varying_dip_5km` (PyRSQSim) | `varying_dip` | Yes |

OBSERVED: `infer_scenario` note — the `"planar"` check runs before `"varying_dip"`. A bundle named `"planar_varying_dip"` would misclassify as `"planar"`. No current bundle has this form, so this is not an active bug.

### `shorten_label()` vs all 14 current bundle IDs

OBSERVED: All 14 bundles shorten without error. However, two bundles produce the **same shortened label**:

| Label | Shortened |
|---|---|
| `HBI:alpine_varying_dip` | `HBI:vDip` |
| `PyRSQSim:alpine_varying_dip_5km` | `pRSQ:vDip` |

These are different sources so their `src` prefix differs — no collision. But:

| Label | Shortened |
|---|---|
| `EQquasi:alpine_varying_dip` | `EQ:vDip` |
| `MCQsim:alpine_varying_dip` | `MCQ:vDip` |
| `RSQSim:alpine_varying_dip` | `RSQ:vDip` |
| `HBI:alpine_varying_dip` | `HBI:vDip` |
| `HBI:alpine_hbi_81102019_varying_dip` | `HBI:81102019-vDip` |
| `HBI:alpine_hbi_82102019_varying_dip` | `HBI:82102019-vDip` |

OBSERVED: the old `HBI:alpine_varying_dip` bundle and the new `HBI:alpine_hbi_81102019_varying_dip` bundle both appear in the same figures with different short labels (`HBI:vDip` vs `HBI:81102019-vDip`). They contain overlapping but different data (different row counts, different `m0_Nm` values). This is confusing and likely unintended duplication.

### DEPTH_PCT consistency

`plot.accumulated.slip.eqquasi.py:57`: `DEPTH_PCT = 88.0` — used at line 116 for actual data extraction (`row = int(round(0.88 * (nrows-1)))`).

OBSERVED: line 172 uses `int(round((DEPTH_PCT / 100.0) * 25))` to compute `depth_row` for the subtitle string only. The literal `25` hardcodes the expected row count. If the grid has fewer than 25 rows, the subtitle would display the wrong depth index. This is display-only and does not affect data correctness.

`benchmark.comparison.analytics.py` does not define DEPTH_PCT; it delegates to the imported `plot.accumulated.slip.eqquasi` module, which uses 88%. The value is therefore consistent between the two scripts.

---

## Section 6 — Logging & Error Handling

### Silent vs loud failures

| Script | Behavior |
|---|---|
| `process.eqquasi.py` | Raises `RuntimeError` if no events meet criteria, if hypocenter is missing for Mw>7 event, if geometry mismatch. Loud. |
| `process.pyrsqsim.py` | Raises `RuntimeError` / `ValueError` on all mismatch conditions. Loud. |
| `process.hbi.py` | Prints skip message for unrecognized geometry type; continues to next job (semi-silent). If a column is absent in the source CSV, `remap_rows` inserts `"nan"` without warning. |
| `benchmark.comparison.analytics.py` | `_find_eqquasi_case_dir` returns `None` silently; caller checks `if case_dir is None: continue` (line 750) — silently skips EQquasi long-term slip rate if case dir is missing. No warning to user. |
| `plot.accumulated.slip.eqquasi.py` | Prints per-directory skip messages to stdout; raises if no profiles loaded at all. |

### Missing file handling

| Script | Behavior |
|---|---|
| `process.eqquasi.py` | `validate_geometry_shape` raises on missing Q dirs. Individual `fault.*.nc` missing: silently no hypocenter (NaN), then raises if event is Mw>7 and hypocenter is NaN. |
| `process.pyrsqsim.py` | `infer_paths` raises `FileNotFoundError` for each missing input before processing begins. |
| `process.hbi.py` | `find_file` raises on missing or ambiguous matches. |
| `benchmark.comparison.analytics.py` | `discover_bundles` raises on incomplete bundles (missing site files). |

### Per-item vs batch failure propagation

- `process.eqquasi.py` and `process.pyrsqsim.py`: failure on any event propagates up and stops the entire export. No partial output on error.
- `process.hbi.py`: job-level loop continues past unrecognized geometry types; stops on missing files.
- `benchmark.comparison.analytics.py`: each plot function is independent; a failure in one figure does not prevent others from completing (no try/except, but each function is called sequentially from `main()`).

---

## Section 7 — Performance & Scaling

### O(N²) patterns

**`plot_gr_law` — OBSERVED O(N²)**  
`benchmark.comparison.analytics.py:719–720`:
```python
mw_unique = np.unique(mw)
n_cumulative = np.array([np.sum(mw >= m) for m in mw_unique])
```
For N events with K unique magnitudes, this is O(N × K) ≈ O(N²) in the worst case. At current scale (tens of events per bundle) this is negligible. For catalogs with thousands of events the vectorized replacement is one line:
```python
n_cumulative = mw_unique.size - np.searchsorted(mw, mw_unique)
```

**`pairwise_event_rows` and `pairwise_site_rows` — OBSERVED redundant I/O**  
`benchmark.comparison.analytics.py:314–345` and `348–383`: for each pair `(lhs, rhs)` of bundles sharing a scenario, both event files are re-loaded from disk with `load_rows()`. With 14 bundles currently across two scenarios (21 + 15 = 36 pairs), each event file is loaded once per pair it participates in (up to 7 times for the most-paired bundle). At current scale this is fast. No in-memory caching.

**Site histograms — OBSERVED triple redundant I/O**  
`plot_site_slip_histograms` loops over all bundles for each of the three sites, loading site files independently per call. With 14 bundles × 3 sites = 42 file reads per analytics run. Acceptable now.

### IO that could be batched

- `benchmark.comparison.analytics.py` calls `load_rows(bundle.event_file, ...)` in at least 5 separate functions (`bundle_summary_rows`, `pairwise_event_rows`, `plot_event_mw_histogram`, `plot_mfd`, `plot_gr_law`, `plot_rupture_extents`, `plot_long_term_slip_rate`, `plot_accumulated_slip_along_strike`). Each function re-opens and re-parses the same CSV. A pre-load pass before calling plot functions would eliminate redundant I/O.

---

## Section 8 — Top-N Priorities

Ranked by (impact on correctness × ease of fix), each actionable in <1 day.

| Priority | Finding | Impact | Fix | File:Line |
|---|---|---|---|---|
| **1** | `process.pyrsqsim.py` uses `SECONDS_PER_YEAR = 3.15e7` instead of `365.25 * 24 * 3600` — 0.18% error, ~18 yr drift per 10 000 yr | High: silently biases `t0_year` and recurrence intervals for PyRSQSim cross-code comparisons | Change `3.15e7` → `365.25 * 24.0 * 3600.0` | `process.pyrsqsim.py:14` |
| **2** | Old HBI bundle `alpine_varying_dip_*` (Apr 21) coexists with new `alpine_hbi_81102019_varying_dip_*` (May 11); analytics loads both as separate HBI entries, inflating HBI representation | High: pairwise and histogram figures show HBI twice for the same raw run | Remove the four `alpine_varying_dip_*` files from `HBI_results/` or add an exclusion rule | `HBI_results/alpine_varying_dip_*` |
| **3** | `benchmark/README.md` documents three PNG names (`event_mw_histogram.png`, `recurrence_cdf.png`, `site{1,2,3}_slip_cdf.png`) that do not match the actual generated filenames | Medium: breaks every developer following the docs to find outputs | Update six filename references in the README | `benchmark/README.md:192–204` |
| **4** | `process.hbi.py` is not wired into any shell script; HBI bundles must be rebuilt manually, and the step is absent from `run.comparison.sh` | Medium: HBI results silently go stale after any HBI raw data update | Add `python benchmark/utils/process.hbi.py` invocation to `run.comparison.sh` before analytics | `run.comparison.sh` |
| **5** | Site-slip-threshold 0.5 m applied to EQquasi via shell script but 0.0 m for all other codes — asymmetric filtering biases cross-code site-slip CDFs | Medium: EQquasi site-event counts are artificially reduced | Either document this as intentional with physical justification, or apply the same threshold to all codes | `process.eqquasi.sh:22,36,49,56`; `process.pyrsqsim.sh` |
| **6** | `benchmark.comparison.analytics.py:633`: `min(sbounds_eq)` will raise `ValueError` if all EQquasi events have NaN `sbound_NZTM_m` | Low-medium: rare edge case, but silent failure mode if bounds are missing | Add guard: `if not sbounds_eq: continue` | `benchmark.comparison.analytics.py:631–633` |
| **7** | `benchmark/README.md` hardcodes absolute paths `/Users/dliu/scratch/alpine/` throughout the Export and Verify sections | Low: documentation only, but breaks copy-paste for any other user | Replace all absolute paths with `${PROJECT_ROOT}/` | `benchmark/README.md:29,46,60,100–145` |
| **8** | `package-lock.json` at project root has no corresponding `package.json` and no JS tooling | Low: clutter, potential confusion | Remove `package-lock.json` and add `*.lock` or `package-lock.json` to `.gitignore` | `/package-lock.json` |
