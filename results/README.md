# Results

This folder is the main `EQquasi` simulation archive for the Alpine Fault project. Each case directory is mostly self-contained and usually includes setup files, restart files, plotting scripts, and cycle outputs in `Q*` folders.

The benchmark layer should read from here, but cross-code conversion and comparison scripts should live in `benchmark/`, not inside these case directories.

## Case Status

- `nz.bp5.qdc.dip50.2000.norm_6mm_yr`
  - planar benchmark case with `24` `Q*` cycle folders
  - this is the active planar source for `benchmark/Simulation_results/EQquasi_results/alpine_planar_*`
- `nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad`
  - no-dip-change benchmark case with `21` `Q*` cycle folders
  - this is the active no-dip-change source for `benchmark/Simulation_results/EQquasi_results/alpine_no_dip_change_*`
- `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad`
  - varying-dip benchmark case with `18` `Q*` cycle folders
  - this is the active varying-dip source for `benchmark/Simulation_results/EQquasi_results/alpine_varying_dip_*`
  - this case uses the newer varying-dip mesh
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad`
  - older varying-dip case with `23` `Q*` cycle folders
  - no longer the shared benchmark source
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad.casedir.only`
  - case-directory shell only
  - no `Q*` results

Other directories here are older or alternate Alpine Fault runs and should be treated as archived cases unless you explicitly want to benchmark or post-process them.

## Post-Processing

`postprocess.sh` runs the full post-processing pipeline for a single case and collects all figure outputs into `<case>/a.results/`.

```bash
bash results/postprocess.sh <case_folder_name>
```

Example:

```bash
bash results/postprocess.sh nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad
```

It runs the following scripts in order:

1. `plotAccumulated` — all variables (slips, slipd, stresses, slip rate), both vertical and horizontal profiles, at 100% depth
2. `plotOnFaultVars` — per-snapshot on-fault variable PNGs and GIF animation in each `Q*` folder
3. `plotSlipAtPaleoSite` — accumulated slip at a paleoseismic station (50% along strike)
4. `post.process.eq.stats.py` — earthquake statistics CSV

All PNGs, GIFs, and the stats CSV are collected into `<case>/a.results/`.

Requirements: the case must have `Q*` cycle folders with `fault.*.nc` files, `fault.r.nc`, `user_defined_params.py`, and `lib.py`.

## Existing Case Scripts

Useful scripts already present inside the populated case directories include:

- `post.process.eq.info.py` or `post.process.eq.stats.py` for cycle-level event statistics
- `plotSlipAtPaleoSite` for stitching free-surface slip across `Q*` cycles at a chosen strike position
- `plot.rupture.length.vs.time.py` for rupture-length and rupture-evolution figures
- `plotProfile.py` and related plotting helpers for on-fault section views

These case-local scripts are the reference logic for how times and cumulative slips are pieced together across `Q0`, `Q1`, `Q2`, and so on.

## Notes

- Each `Q*` folder should be interpreted as one earthquake-cycle segment.
- `Q(n+1)` continues from the restart state of `Qn`, so catalogue time must be accumulated across cycles.
- Free-surface fault patches are on the maximum row index of the fault grid.
- Many publication-style PNG outputs already exist directly in the case directories.
- For the local git release workflow, the full case directories remain local-only. This README and the small top-level utility scripts are the tracked summary layer for `results/`.

## Inventory

Current top-level case inventory:

- `nz.bp5.qdc.dip50.2000.norm0`: `11` `Q*`, `36` root PNGs
- `nz.bp5.qdc.dip50.2000.norm_6mm_yr`: `24` `Q*`, `37` root PNGs
- `nz.bp5.qdc.dip90.4000`: `2` `Q*`, `1` root PNG
- `nz.bp5.qdc.noDipChange.2000.norm0.slowInitialLoad`: `28` `Q*`, `65` root PNGs
- `nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad`: `21` `Q*`, `23` root PNGs
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad`: `23` `Q*`, `30` root PNGs
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad.casedir.only`: `0` `Q*`, `0` root PNGs
- `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad`: `18` `Q*`, `0` root PNGs

Legacy utilities that were formerly under `post_utility_dev/` were intentionally reduced to a small reference archive at `archive/post_utility_dev_varyDip_legacy/`. The old scratch tree is no longer part of the active results workflow.
