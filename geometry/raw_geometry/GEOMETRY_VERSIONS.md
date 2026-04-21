# Geometry Versions Documentation

## Current Active Geometries
**Location:** `post_workshop_geoms/`
**Date:** November 2025 (post-workshop)
**Source:** Andy Howell's email (2025-11-17)
**Status:** Latest - Ready for new simulations

### Files:
- **`alpine_planar.stl`** (10M) - Updated planar fault geometry
- **`alpine_planar.geo`** (1.2K) - GMSH script for planar fault
- **`trace_points_alpine_planar.geo`** (563B) - Surface trace points definition
- **`variable_dip_1km.stl`** (9.9M) - Updated variable dip geometry
- **`variable_dip_renumbered.geo`** (17K) - GMSH script for variable dip
- **`GNSSNZ_comb_v10_sillremoved_19-Aug-2025.dat`** (140K) - GNSS velocity data (Australian Plate-fixed frame)

**Key Updates:**
- Cleaner meshes from Avi's improved .geo scripts
- Two .geo approaches for planar: 4-corners vs average strike/dip (average approach used - cleaner mesh)
- Coordinates in metres, NZTM projection (EPSG:2193)
- Suitable for trimming to shorter fault lengths if needed (maintain southern paleoseismic sites)

---

## Previous Geometries (Archive)
**Location:** `geoms_202508/`
**Date:** July-August 2025
**Status:** Legacy - Keep for reference/reproducibility

### Files:
- **`alpine_planar500m.stl`** (31M) - High-resolution planar model
- **`no_dip_change_1km.stl`** (10M) - Variable geometry, constant dip
- **`variable_dip_simpler_geom1km_hws.stl`** (9.8M) - Variable dip with HWS structure

**Notes:**
- These were used to generate processed_geometry/\*.txt outputs
- Results in .res2km and .res4km folders based on these geometries
- Keep for reproducibility and reference

---

## GNSS Data

**File:** `post_workshop_geoms/GNSSNZ_comb_v10_sillremoved_19-Aug-2025.dat`

- **Format:** Text file with CRLF line terminators
- **Columns:** lon lat ve vn vu se sn su name
  - lon, lat: Station coordinates (NZTM projection, metres)
  - ve, vn, vu: Velocity east, north, up (mm/yr)
  - se, sn, su: Velocity uncertainties
  - name: Station identifier

- **Reference Frame:** Australian Plate-fixed
- **Scale Reference:** Velocity at Christchurch = 34 mm/yr
- **Coverage:** All New Zealand
- **Status:** Can be subset for current exercise or used for model validation

---

## Recommendation

Use **`post_workshop_geoms/`** for new simulations and model development.
Keep **`geoms_202508/`** as archive for reference.
