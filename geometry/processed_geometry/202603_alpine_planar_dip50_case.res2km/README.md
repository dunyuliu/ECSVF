# Alpine Planar Dip50 Case Geometry

This folder contains the processed geometry artifact used by the benchmark planar `EQquasi` case:

- `results/nz.bp5.qdc.dip50.2000.norm_6mm_yr`

The grid file:

- `alpine_planar_dip50_case.res2km_grid.npz`

was created by cropping the centered `-320 km` to `320 km` strike window from:

- `geometry/processed_geometry/202512_alpine_planar.res2km/alpine_planar.res2km_grid.npz`

This keeps the same planar benchmark transform metadata and STL reference while matching the case-local `321 x 26` fault grid used by the `dip50` simulation archive.
