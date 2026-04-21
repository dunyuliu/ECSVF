# Alpine Geometry - Earthquake Fault Modeling

3D earthquake fault geometry analysis and visualization toolkit for modeling fault surface topography, strain zones, and seismic rupture behavior.

## Overview

This project provides tools to:
- Load and visualize 3D fault surface geometry from STL and geometry text files
- Analyze fault traces (surface intersections)
- Define and visualize seismic zones (velocity-weakening, velocity-strengthening, creeping)
- Generate earthquake rupture models with variable fault dip angles
- Transform coordinates between geographic and fault-aligned systems

## Key Files

### Visualization & Modeling
- **`earthquake_model_3d.py`** - Main visualization script for 3D fault models with three cases:
  - Case 1: Simple planar fault geometry (50° dip)
  - Case 2: Variable fault geometry (no dip change)
  - Case 3: Variable dip with horizontal wing structure (HWS)

- **`fault_trace_analysis.py`** - Analyzes fault surface traces:
  - Extracts top trace (surface intersection) from STL files
  - Calculates strike direction and fault geometry
  - Rotates coordinates to align with strike direction
  - Transforms station coordinates for rupture modeling

### Geometry Processing
- **`stl_transform.py`** - STL file transformations and analysis
- **`stl_to_eqquasi_fault_grid.py`** - Converts STL fault surfaces to earthquake quasi-static grid format
- **`test_stl.py`** - Unit tests for STL processing

### Data Files
- **STL Files** (binary 3D models):
  - `alpine_planar500m.stl` - High-resolution planar fault model
  - `no_dip_change_1km.stl` - Variable geometry with constant dip
  - `variable_dip_simpler_geom1km_hws.stl` - Variable dip with HWS structure

- **Geometry Files** (fault surface data):
  - `*_Geometry.txt` - Processed fault geometry with coordinates and derivatives
  - Format: Grid points with x, z positions and dy/dx, dy/dz slope values

- **Results Directories**:
  - `.res2km` and `.res4km` - Output from earthquake simulation at different resolutions

### Visualization Outputs
- `earthquake_model_case*.png` - 3D fault diagrams with zone definitions (640×800 each)
- `fault_trace_analysis.png` - Original and rotated fault traces
- `station_coordinates_analysis.png` - Station locations in geographic and fault-aligned systems

## Usage

### Generate 3D Fault Models
```bash
# Planar fault example
python earthquake_model_3d.py case1

# Variable geometry models
python earthquake_model_3d.py case2
python earthquake_model_3d.py case3
```

### Analyze Fault Traces
```bash
python fault_trace_analysis.py
```
Analyzes STL geometry, finds fault surface trace, calculates strike angle, and transforms coordinates.

### Process STL Files
```bash
python stl_transform.py
python stl_to_eqquasi_fault_grid.py
```

## Fault Model Definitions

Models include three seismic zones:

- **VW (Velocity-Weakening)**: Locked seismogenic zone (depth: 2-20 km, 600 km strike length)
  - Nucleates earthquakes, high stress accumulation

- **VS (Velocity-Strengthening)**: Stable sliding zone around VW
  - 2 km depth extension above/below VW
  - 2 km lateral extension on ±x sides of VW
  - Transition zone with intermediate friction

- **Creeping**: Ductile flow zone
  - Shallow (above 2 km) and deep (below 20 km) regions
  - Continuous creep, minimal earthquake rupture

## Geometry Versions

**Current (December 2025):** `raw_geometry/post_workshop_geoms/`
- Updated planar and variable-dip meshes from Andy Howell's workshop (Nov 2025)
- GMSH `.geo` scripts for mesh generation/modification
- GNSS velocity data for model validation/loading
- Processed into 2km resolution grids: `202512_alpine_planar.res2km/`, `202512_variable_dip_1km.res2km/`
- Includes a case-matched planar crop for the `EQquasi` `dip50` benchmark run: `202603_alpine_planar_dip50_case.res2km/`
- See `raw_geometry/GEOMETRY_VERSIONS.md` for details

**Legacy (July-August 2025):** `raw_geometry/geoms_202508/`
- Original geometries used for initial simulations
- Processed into 2km and 4km resolution grids: `202508_*` folders in `processed_geometry/`
- Kept for reproducibility and reference

## Coordinate System

Models use a 3D coordinate system:
- **X-axis** (Strike direction): Horizontal, along fault trace (~E-W)
- **Y-axis** (Dip/Normal direction): Horizontal, perpendicular to strike (~N-S)
- **Z-axis** (Depth): Vertical, negative downward (0 = surface, -50 km = model base)

Geographic coordinates are rotated to align with fault strike for analysis.

## Dependencies

- NumPy - Numerical computations
- Matplotlib - 3D visualization
- Python 3.6+

## Project Structure

```
.
├── utils/                           # Operating scripts
│   ├── earthquake_model_3d.py       # Main visualization
│   ├── fault_trace_analysis.py      # Trace extraction & rotation
│   ├── stl_transform.py             # STL transformations
│   └── stl_to_eqquasi_fault_grid.py # STL to grid conversion
├── raw_geometry/                    # Geometry input files
│   ├── post_workshop_geoms/         # Latest (Nov 2025) - Andy Howell's meshes
│   │   ├── alpine_planar.stl
│   │   ├── alpine_planar.geo
│   │   ├── variable_dip_1km.stl
│   │   ├── variable_dip_renumbered.geo
│   │   ├── trace_points_alpine_planar.geo
│   │   └── GNSSNZ_comb_v10_sillremoved_19-Aug-2025.dat
│   ├── geoms_202508/                # Legacy (Jul-Aug 2025) - Original geometries
│   │   ├── alpine_planar500m.stl
│   │   ├── no_dip_change_1km.stl
│   │   └── variable_dip_simpler_geom1km_hws.stl
│   └── GEOMETRY_VERSIONS.md         # Detailed geometry documentation
├── processed_geometry/              # Processed geometry scenarios (dated prefixes)
│   ├── 202508_no_dip_change_1km.res2km/     # Aug 2025 (2km)
│   ├── 202508_no_dip_change_1km.res4km/     # Aug 2025 (4km)
│   ├── 202508_variable_dip_simpler_geom1km_hws.res2km/  # Aug 2025 (2km)
│   ├── 202512_alpine_planar.res2km/         # Dec 2025 planar (2km)
│   ├── 202512_variable_dip_1km.res2km/      # Dec 2025 variable dip (2km)
│   ├── no_dip_change_1km_Geometry.txt       # Original geometry file
│   ├── bFault_Rough_Geometry_noDip.txt
│   └── bFault_Rough_Geometry_varyDip.txt
├── *.res2km, *.res4km/              # Simulation results (2km & 4km resolution)
├── nz.bp5.qdc.dip90.4000/           # Original simulation archive
├── processing_scenarios_summary.txt # Processing documentation
└── README.md                        # This file
```

## References

Fault geometry models represent Alpine Fault (New Zealand) rupture scenarios with variable dip angles and seismic zone definitions for earthquake dynamic rupture simulations.
