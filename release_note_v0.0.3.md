# Release Note v0.0.3

Date: 2026-05-18

## Summary of Scope

Plot clarity improvements: added `--max-per-source-scenario` and `--exclude` filters to
the analytics script so that sources with multiple runs (HBI has 6) can be capped to one
representative bundle per scenario, and specific cases can be suppressed entirely.
`run.comparison.sh` now defaults to one HBI bundle per scenario and excludes the
EQquasi no-dip-change case from the comparison figures. AUDIT.md refreshed by a full
re-audit pass.

## Files Added

None.

## Files Modified

| File | Change |
|------|--------|
| `benchmark/utils/benchmark.comparison.analytics.py` | Added `--max-per-source-scenario N` and `--exclude PATTERN…` CLI arguments; added `cap_bundles()` helper; filtering applied before color assignment in `main()` |
| `run.comparison.sh` | Added `--max-per-source-scenario 1 --exclude no_dip_change` to analytics call |
| `AUDIT.md` | Replaced with fresh full-project audit (victor-reyes pass, 2026-05-18) |

## Files Archived

| File | Destination |
|------|-------------|
| `release_note_v0.0.2.md` | `docs/release_note_v0.0.2.md` |

## Content Updates to Master Documents

- `AUDIT.md` regenerated — previous version was from the v0.0.2 audit pass (2026-05-11).

## Audit Findings and Fixes (PROJECT_RULES.md)

| Rule | Check | Result |
|------|-------|--------|
| Rule 3 — Bundle completeness | All `*_event_info.csv` bundles have matching site1/2/3 files | PASS |
| Rule 4 — No stale duplicates | No two bundles cover the same model run in the same `*_results/` dir | PASS |
| Rule 5 — No hardcoded user paths | No `/Users/<username>/` in any tracked `.py` or `.sh` | PASS |
| Rule 1 — Verify checks all PASS | Requires live EQquasi pipeline run | Pending |
| Rule 2 — Pipeline runs without error | Requires live pipeline run | Pending |

## Remaining Open Items

- **Rules 1 & 2**: Run `bash run.comparison.sh` end-to-end and
  `python benchmark/utils/verify.eqquasi.geometry.conversion.py` per EQquasi case to
  confirm all five `[PASS]` checks pass. Not automatable at release time.
- **PyRSQSim bundle absent**: Local Alpine example is a stub (1 event). Re-sync from Knox
  with `bash utig.to.local.sh` when the full simulation is available.
- **Tandem_results empty**: No Tandem results received yet.
- **EQquasi alpine_varying_dip**: Still accumulating events on UTIG Knox. Current export
  has 11 events in 3000 yr vs 41 (RSQSim) and 69 (MCQsim) — re-sync periodically.
- **MCQsim varying-dip slip rate outlier**: Long-term slip rate ~80–95 mm/yr vs geological
  target ~27 mm/yr. Flagged for follow-up with the MCQsim team.
- **HBI bundle selection**: Default is alphabetically first per scenario (`81102019_varying_dip`,
  `81103019_planar`). Override via `--exclude` / `--max-per-source-scenario` if a different
  representative run is preferred.
