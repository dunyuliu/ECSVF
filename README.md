# New Zealand Alpine Fault EQquasi Project

This repository is the working project for New Zealand Alpine Fault earthquake-cycle simulations run with `EQquasi`. The main body of the repo is your geometry preparation, case setup, simulation output, and post-processing workflow for Alpine Fault models. The `benchmark/` folder is the exception: it is reserved for cross-code benchmarking and intercomparison work.

## What This Repo Contains

- `geometry/`
  Geometry generation, STL processing, trace analysis, and conversion to `EQquasi` fault-grid inputs. This is where Alpine Fault surfaces are turned into `bFault_Rough_Geometry.txt`-style inputs and documented processed geometries.

- `results/`
  The main simulation archive. Each `results/<case>/` directory is a mostly self-contained `EQquasi` case directory with:
  - setup/configuration scripts such as `case.setup`, `create.newcase`, `user_defined_params.py`, and `defaultParameters.py`
  - launch scripts such as `run.sh`, `batch.hpc`, and `batch.cycle.eqquasi.hpc`
  - mesh and input files such as `fault.r.nc`, `disp.r.nc`, `eqquasi.mesh.*.nc`, `model.txt`, `stations.txt`, and `on_fault_vars_input.nc`
  - per-cycle outputs in `Q*` folders, typically containing `fault.*.nc`, `disp.*.nc`, `global.dat`, `cplot*`, `fltst*`, `srfst*`, and related post-processed products
  - analysis and plotting utilities such as `plotAccumulated`, `plotOnFaultVars`, `plotSlipAtPaleoSite`, `plot.rupture.length.vs.time.py`, and `post.process.eq.stats.py`

- `benchmark/`
  Cross-code benchmarking and intercomparison layer. Contains benchmark scripts, imported simulator result bundles, generated comparison outputs, reference docs, and archive material. Active codes in the comparison: `EQquasi`, `RSQSim`, `PyRSQSim`, `MCQsim`, `HBI`. Key utilities:
  - `benchmark/utils/process.eqquasi.py` — exports EQquasi Q-folder outputs into Avi-style benchmark CSVs; uses 0.5 m slip threshold and surface-row northing for rupture extents
  - `benchmark/utils/benchmark.comparison.analytics.py` — generates comparison figures for all codes (Mw histogram, recurrence CDF, site slip CDFs, rupture extents, long-term slip rate)
  - `benchmark/utils/plot.accumulated.slip.eqquasi.py` — standalone EQquasi post-processing: individual slip profiles, accumulated slip, and slip rate vs. fault strike for a single case

## Current Simulation Cases

Representative case directories under `results/` include:

- `nz.bp5.qdc.dip50.2000.norm0`
- `nz.bp5.qdc.dip50.2000.norm_6mm_yr`
- `nz.bp5.qdc.dip90.4000`
- `nz.bp5.qdc.noDipChange.2000.norm0.slowInitialLoad`
- `nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad`
- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad`
- `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad`

There are also partially staged or archive-like directories such as:

- `nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad.casedir.only`

The active varying-dip benchmark source is `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad`, which uses the newer mesh.

## Typical Workflow

### 1. Build or inspect Alpine Fault geometry

Geometry work lives in [`geometry/`](./geometry). The existing geometry README is the detailed reference for that part of the workflow.

Useful entry points:

- `geometry/utils/stl_to_eqquasi_fault_grid.py`
- `geometry/utils/fault_trace_analysis.py`
- `geometry/utils/earthquake_model_3d.py`
- `geometry/processed_geometry/`

Typical outputs are processed grids, summary reports, plots, and `bFault_Rough_Geometry.txt` files that can be used by `EQquasi`.

### 2. Stage geometry into an EQquasi case

Case directories under `results/` expect geometry and model inputs in the case root. The geometry-to-HPC handoff is represented by scripts such as:

- `geometry/transfer.geometry.to.ls6.sh`

and by the case setup scripts inside each results directory.

### 3. Configure a case

Within a case directory, the main setup logic is typically:

- `user_defined_params.py` for case parameters
- `case.setup` to generate:
  - `model.txt`
  - `stations.txt`
  - `on_fault_vars_input.nc`
  - `batch.hpc`
  - `batch.cycle.eqquasi.hpc`
  - `run.sh`

### 4. Run the simulation

Cases are launched either locally through MPI or on HPC. The repo includes both local and batch-style launch scripts.

Examples:

```bash
cd results/nz.bp5.qdc.noDipChange.2000.norm0.slowInitialLoad
python case.setup
bash run.sh
```

or on the cluster using generated batch scripts:

```bash
sbatch batch.hpc
sbatch batch.cycle.eqquasi.hpc
```

### 5. Archive outputs by cycle

`EQquasi` outputs are organized into `Q*` cycle folders. A typical `Q*` folder contains:

- `global.dat`
  Time history and scalar diagnostics used by multiple post-processing scripts.

- `fault.*.nc`
  On-fault field snapshots on the fault grid. Representative variables include:
  - `slip_rate`
  - `slips`
  - `slipd`
  - `slipn`
  - `shear_strike`
  - `shear_dip`
  - `effective_normal`

- `disp.*.nc`
  Off-fault displacement snapshots, typically stored as a `disp(node_id, disp_dim_id)` array.

- `fltst*`, `srfst*`, `tdyna*`, `cplot*`
  Additional simulator outputs and post-processing artifacts.

### 6. Sync results back from HPC

The root-level utility:

- `results/sync_models.py`

pulls selected case directories back from TACC/LS6 using `scp` or `rsync`.

For the varying-dip case running on UTIG Knox:

```bash
bash utig.to.local.sh
```

This rsyncs `nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad` from `dliu@knox.ig.utexas.edu` into `results/`. After syncing, re-run `process.eqquasi.sh` and the analytics to pick up any new Q folders.

## Post-Processing and Analysis

The repo already includes several analysis scripts in `results/` and duplicated inside case directories. Common uses include:

- event statistics: `post.process.eq.stats.py`
- rupture-length evolution: `plot.rupture.length.vs.time.py`
- accumulated slip or stress views: `plotAccumulated`, `plotOnFaultVars`, `plotProfiles`
- paleo-site slip diagnostics: `plotSlipAtPaleoSite`

These scripts work directly against the archived `Q*` folders and their NetCDF outputs.

## Dependencies

The Python tooling in this repo expects at least:

- Python 3
- `numpy`
- `matplotlib`
- `xarray`
- `netCDF4`

Operationally, the simulation workflow also assumes:

- the `eqquasi` executable is available in the execution environment
- MPI launch tools such as `mpirun.openmpi` or `ibrun`
- cluster modules such as `mumps` for the provided batch scripts

For plotting-heavy runs on shared systems, it is useful to set a writable `MPLCONFIGDIR`.

## Repository Notes

- This is a data-heavy working repository. `results/` is by far the largest component and contains many generated files.
- Case directories intentionally duplicate setup and plotting scripts so each case is relatively self-contained and reproducible.
- `benchmark/` should stay separate from the main `EQquasi` archive so intercomparison code does not get mixed into case directories.
- For the local git release workflow, large result archives stay local-only and are documented rather than committed wholesale.

## Suggested Top-Level Mental Model

Think of the repo as:

1. `geometry/` builds Alpine Fault inputs.
2. `results/` runs and archives `EQquasi` experiments.
3. `benchmark/` is the cross-simulator comparison layer.

## Related Documentation

- [`geometry/README.md`](./geometry/README.md)
- `geometry/processing_scenarios_summary.txt`
- the per-case setup and plotting scripts under `results/<case>/`
