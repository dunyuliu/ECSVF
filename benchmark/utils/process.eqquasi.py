#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr


SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
EVENT_WINDOW_SIZE = 20
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = REPO_ROOT / "results/nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmark/Simulation_results/EQquasi_results"


@dataclass(frozen=True)
class Site:
    label: str
    file_suffix: str
    x_nztm_m: float
    y_nztm_m: float


DEFAULT_SITES = (
    Site("site_1", "site1_event_info.csv", 1207801.604, 5071777.726),
    Site("site_2", "site2_event_info.csv", 1382372.944, 5199926.971),
    Site("site_3", "site3_event_info.csv", 1536168.087, 5311181.708),
)


@dataclass
class GeometryModel:
    local_x_km: np.ndarray
    local_y_km: np.ndarray
    local_z_km: np.ndarray
    global_x_m: np.ndarray
    global_y_m: np.ndarray
    global_z_m: np.ndarray
    area_m2: np.ndarray
    rotation_matrix: np.ndarray
    center_xy_km: np.ndarray
    top_z_km: float

    def global_to_local(self, x_m: float, y_m: float) -> tuple[float, float]:
        xy_km = np.array([x_m, y_m], dtype=float) / 1000.0
        rotated_xy_km = xy_km @ self.rotation_matrix.T
        local_xy_km = rotated_xy_km - self.center_xy_km
        return float(local_xy_km[0]), float(local_xy_km[1])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export EQquasi results as Avi-style benchmark CSV files.")
    parser.add_argument("--input-dir", "--case-dir", dest="input_dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", "--outdir", dest="output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--geometry-npz", type=Path, default=None)
    parser.add_argument("--prefix", type=str, default=None)
    parser.add_argument("--min-mw", type=float, default=7.0)
    parser.add_argument("--slip-rate-threshold", type=float, default=0.1)
    parser.add_argument("--rupture-bound-slip-threshold", type=float, default=0.5)
    parser.add_argument("--slip-eps", type=float, default=0.0)
    parser.add_argument("--site-slip-threshold", type=float, default=0.0)
    return parser.parse_args()


def natural_sort_key(text: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def infer_prefix(case_name: str) -> str:
    lowered = case_name.lower()
    if "nodipchange" in lowered or "no_dip_change" in lowered:
        return "alpine_no_dip_change"
    if "varydip" in lowered or "varying_dip" in lowered:
        return "alpine_varying_dip"
    if "dip50" in lowered or "planar" in lowered:
        return "alpine_planar"
    return re.sub(r"[^a-zA-Z0-9]+", "_", case_name).strip("_")


def resolve_geometry_npz(case_dir: Path, geometry_npz: Path | None) -> Path:
    if geometry_npz is not None:
        return geometry_npz

    case_name = case_dir.name
    hints = (
        (
            "noDipChange",
            REPO_ROOT
            / "geometry/processed_geometry/202508_no_dip_change_1km.res2km/"
            / "no_dip_change_1km.res2km_grid.npz",
        ),
        (
            "varyDip20251202",
            REPO_ROOT
            / "geometry/processed_geometry/202512_variable_dip_1km.res2km/"
            / "variable_dip_1km.res2km_grid.npz",
        ),
        (
            "varyDip",
            REPO_ROOT
            / "geometry/processed_geometry/202508_variable_dip_simpler_geom1km_hws.res2km/"
            / "variable_dip_simpler_geom1km_hws.res2km_grid.npz",
        ),
        (
            "dip50",
            REPO_ROOT
            / "geometry/processed_geometry/202512_alpine_planar.res2km/"
            / "alpine_planar.res2km_grid.npz",
        ),
    )
    for token, npz_path in hints:
        if token in case_name:
            return npz_path
    raise FileNotFoundError(f"Could not infer geometry for {case_name}. Pass --geometry-npz.")


def load_geometry_model(npz_path: Path) -> GeometryModel:
    obj = np.load(npz_path, allow_pickle=True)
    metadata = obj["metadata"].item()
    nx, nz = map(int, metadata["grid_dimensions"])

    local_x_km = obj["x"].reshape(nz, nx)
    local_y_km = obj["y"].reshape(nz, nx)
    local_z_km = obj["z"].reshape(nz, nx)
    dy_dx = obj["dy_dx"].reshape(nz, nx)
    dy_dz = obj["dy_dz"].reshape(nz, nx)

    resolution_km = float(metadata["resolution"])
    area_m2 = np.sqrt(1.0 + dy_dx**2 + dy_dz**2) * (resolution_km * 1000.0) ** 2

    center_x_km, center_y_km, top_z_km = metadata["processing_info"]["center_offsets"]
    center_xy_km = np.array([float(center_x_km), float(center_y_km)], dtype=float)
    top_z_km = float(top_z_km)

    angle_rad = math.radians(float(metadata["rotation_angle"]))
    rotation_matrix = np.array(
        [
            [math.cos(angle_rad), -math.sin(angle_rad)],
            [math.sin(angle_rad), math.cos(angle_rad)],
        ]
    )

    rotated_xy_km = np.stack(
        [local_x_km + center_xy_km[0], local_y_km + center_xy_km[1]],
        axis=-1,
    )
    original_xy_km = rotated_xy_km @ rotation_matrix

    return GeometryModel(
        local_x_km=local_x_km,
        local_y_km=local_y_km,
        local_z_km=local_z_km,
        global_x_m=original_xy_km[..., 0] * 1000.0,
        global_y_m=original_xy_km[..., 1] * 1000.0,
        global_z_m=(local_z_km + top_z_km) * 1000.0,
        area_m2=area_m2,
        rotation_matrix=rotation_matrix,
        center_xy_km=center_xy_km,
        top_z_km=top_z_km,
    )


def validate_geometry_shape(case_dir: Path, geometry: GeometryModel) -> None:
    first_cycle = cycle_dirs(case_dir)
    if not first_cycle:
        raise FileNotFoundError(f"No Q* folders found in {case_dir}")
    with xr.open_dataset(first_cycle[0] / "fault.r.nc") as ds:
        grid_shape = tuple(int(v) for v in ds["shear_strike"].shape)
    if grid_shape != geometry.local_x_km.shape:
        raise ValueError(f"Geometry grid {geometry.local_x_km.shape} does not match {first_cycle[0] / 'fault.r.nc'} shape {grid_shape}")


def cycle_dirs(case_dir: Path) -> list[Path]:
    return sorted(
        [path for path in case_dir.iterdir() if path.is_dir() and path.name.startswith("Q")],
        key=lambda path: natural_sort_key(path.name),
    )


def q_index(q_dir: Path) -> int:
    return int(q_dir.name[1:])


def fault_snapshot_files(q_dir: Path) -> list[Path]:
    files = [path for path in q_dir.glob("fault.*.nc") if re.fullmatch(r"fault\.(\d+)\.nc", path.name)]
    return sorted(files, key=lambda path: int(re.fullmatch(r"fault\.(\d+)\.nc", path.name).group(1)))


def read_global_dat(q_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(q_dir / "global.dat")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]


def find_event_window(times_s: np.ndarray, peak_slip_rates: np.ndarray, slip_rate_threshold: float) -> tuple[int, int] | None:
    start_step = None
    end_step = None

    for i in range(len(peak_slip_rates) - EVENT_WINDOW_SIZE):
        if np.all(peak_slip_rates[i : i + EVENT_WINDOW_SIZE] > slip_rate_threshold):
            start_step = i + 1
            break

    if start_step is None:
        return None

    for i in range(start_step, len(peak_slip_rates) - EVENT_WINDOW_SIZE):
        if np.all(peak_slip_rates[i : i + EVENT_WINDOW_SIZE] < slip_rate_threshold):
            end_step = i + 1
            break

    if end_step is None:
        return None
    return start_step, end_step


def build_cycle_offsets(case_dir: Path) -> dict[int, float]:
    offsets: dict[int, float] = {}
    total_s = 0.0
    for q_dir in cycle_dirs(case_dir):
        times_s, _ = read_global_dat(q_dir)
        offsets[q_index(q_dir)] = total_s
        total_s += float(times_s[-1])
    return offsets


def load_material_properties(model_txt: Path) -> float:
    lines = [line.strip() for line in model_txt.read_text().splitlines() if line.strip()]
    if len(lines) < 8:
        raise ValueError(f"Unexpected model.txt format in {model_txt}")
    vp_m_s, vs_m_s, rho_kg_m3 = map(float, lines[7].split())
    del vp_m_s
    return rho_kg_m3 * vs_m_s * vs_m_s


def file_step_index(path: Path) -> int:
    match = re.fullmatch(r"fault\.(\d+)\.nc", path.name)
    if match is None:
        raise ValueError(f"Unexpected fault snapshot name: {path.name}")
    return int(match.group(1))


def compute_mw(moment_nm: float) -> float:
    if moment_nm <= 0.0:
        return float("nan")
    return (2.0 / 3.0) * math.log10(moment_nm * 1.0e7) - 10.7


def largest_contiguous_strike_segment(mask_2d: np.ndarray) -> np.ndarray:
    """Return a 2D mask restricted to the largest contiguous along-strike segment.

    A column is considered active if any row in that column is True.  The
    function finds all contiguous runs of active columns, selects the longest
    one, and returns a copy of mask_2d with all other columns zeroed out.
    """
    col_active = np.any(mask_2d, axis=0)  # (nid_strike,)
    n_cols = col_active.size

    best_start, best_len = 0, 0
    cur_start, cur_len = 0, 0
    for j in range(n_cols):
        if col_active[j]:
            if cur_len == 0:
                cur_start = j
            cur_len += 1
            if cur_len > best_len:
                best_len = cur_len
                best_start = cur_start
        else:
            cur_len = 0

    result = np.zeros_like(mask_2d)
    if best_len > 0:
        result[:, best_start : best_start + best_len] = mask_2d[:, best_start : best_start + best_len]
    return result


def surface_row_index(geometry: GeometryModel) -> int:
    row = int(np.argmax(np.mean(geometry.local_z_km, axis=1)))
    surface_z_km = geometry.local_z_km[row, :]
    if not np.allclose(surface_z_km, np.max(geometry.local_z_km), atol=1.0e-9):
        raise ValueError("Could not identify a unique free-surface row from the processed geometry")
    return row


def find_nearest_patch_indices(geometry: GeometryModel, sites: tuple[Site, ...]) -> dict[str, tuple[int, int]]:
    mapping: dict[str, tuple[int, int]] = {}
    surface_row = surface_row_index(geometry)
    for site in sites:
        site_x_km, site_y_km = geometry.global_to_local(site.x_nztm_m, site.y_nztm_m)
        distance_sq = (geometry.local_x_km[surface_row, :] - site_x_km) ** 2 + (geometry.local_y_km[surface_row, :] - site_y_km) ** 2
        strike_index = int(np.argmin(distance_sq))
        mapping[site.label] = (surface_row, strike_index)
    return mapping


def extract_cycle_details(
    q_dir: Path,
    geometry: GeometryModel,
    site_patch_indices: dict[str, tuple[int, int]],
    shear_modulus_pa: float,
    slip_rate_threshold: float,
    rupture_bound_slip_threshold: float,
    slip_eps: float,
    site_slip_threshold: float,
) -> dict | None:
    snapshots = fault_snapshot_files(q_dir)
    times_s, peak_slip_rates = read_global_dat(q_dir)
    event_window = find_event_window(times_s, peak_slip_rates, slip_rate_threshold)
    if event_window is None:
        return None
    start_step, end_step = event_window

    start_time_s = float(times_s[start_step - 1])
    end_time_s = float(times_s[end_step - 1])

    with xr.open_dataset(q_dir / "fault.r.nc") as ds:
        end_slips = ds["slips"].values
        end_slipd = ds["slipd"].values
    if end_slips.shape != geometry.local_x_km.shape or end_slipd.shape != geometry.local_x_km.shape:
        raise ValueError(f"Shape mismatch between {q_dir / 'fault.r.nc'} and geometry grid")
    end_slip_mag = np.hypot(end_slips, end_slipd)
    rupture_mask = end_slip_mag > slip_eps
    if not np.any(rupture_mask):
        return None

    bound_mask = largest_contiguous_strike_segment(end_slip_mag > rupture_bound_slip_threshold)
    surf_row = surface_row_index(geometry)
    bound_cols = np.any(bound_mask, axis=0)
    rupture_area_m2 = float(np.sum(geometry.area_m2[rupture_mask]))
    scalar_moment_nm = float(shear_modulus_pa * np.sum(geometry.area_m2[rupture_mask] * end_slip_mag[rupture_mask]))

    hypocenter_idx = None
    for snapshot in snapshots:
        step = file_step_index(snapshot)
        if step < start_step or step > end_step:
            continue
        with xr.open_dataset(snapshot) as ds:
            slip_rate = ds["slip_rate"].values
        if slip_rate.shape != geometry.local_x_km.shape:
            raise ValueError(f"Shape mismatch between {snapshot} and geometry grid")
        if np.any(slip_rate > slip_rate_threshold):
            hypocenter_idx = np.unravel_index(int(np.nanargmax(slip_rate)), slip_rate.shape)
            break
    site_rows: dict[str, dict] = {}
    for site in DEFAULT_SITES:
        patch_idx = site_patch_indices[site.label]
        slip_m = float(end_slip_mag[patch_idx])
        if slip_m < site_slip_threshold:
            continue
        rake_deg = (math.degrees(math.atan2(float(end_slipd[patch_idx]), float(end_slips[patch_idx]))) + 360.0) % 360.0
        site_rows[site.label] = {"Slip_m": slip_m, "Rake_deg": rake_deg}

    return {
        "start_time_s": start_time_s,
        "duration_s": end_time_s - start_time_s,
        "m0_Nm": scalar_moment_nm,
        "mw": compute_mw(scalar_moment_nm),
        "area_m2": rupture_area_m2,
        "max_slip_m": float(np.max(end_slip_mag[rupture_mask])),
        "x_NZTM_m": float(geometry.global_x_m[hypocenter_idx]) if hypocenter_idx is not None else float("nan"),
        "y_NZTM_m": float(geometry.global_y_m[hypocenter_idx]) if hypocenter_idx is not None else float("nan"),
        "z_NZTM_m": float(geometry.global_z_m[hypocenter_idx]) if hypocenter_idx is not None else float("nan"),
        # Surface row only: deep patches are displaced down-dip so their northing
        # diverges from the surface trace; bounds should reflect the surface extent.
        "sbound_NZTM_m": float(np.min(geometry.global_y_m[surf_row, bound_cols])) if np.any(bound_cols) else float("nan"),
        "nbound_NZTM_m": float(np.max(geometry.global_y_m[surf_row, bound_cols])) if np.any(bound_cols) else float("nan"),
        "site_rows": site_rows,
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_catalog_rows(
    case_dir: Path,
    geometry: GeometryModel,
    sites: tuple[Site, ...],
    min_mw: float,
    slip_rate_threshold: float,
    rupture_bound_slip_threshold: float,
    slip_eps: float,
    site_slip_threshold: float,
) -> tuple[dict[str, tuple[int, int]], list[dict], dict[str, list[dict]]]:
    validate_geometry_shape(case_dir, geometry)
    shear_modulus_pa = load_material_properties(case_dir / "model.txt")
    patch_indices = find_nearest_patch_indices(geometry, sites)
    cycle_offsets = build_cycle_offsets(case_dir)
    q_lookup = {q_index(q_dir): q_dir for q_dir in cycle_dirs(case_dir)}

    event_rows: list[dict] = []
    site_rows_by_label = {site.label: [] for site in sites}

    for event_id in sorted(q_lookup):
        q_dir = q_lookup[event_id]
        details = extract_cycle_details(
            q_dir=q_dir,
            geometry=geometry,
            site_patch_indices=patch_indices,
            shear_modulus_pa=shear_modulus_pa,
            slip_rate_threshold=slip_rate_threshold,
            rupture_bound_slip_threshold=rupture_bound_slip_threshold,
            slip_eps=slip_eps,
            site_slip_threshold=site_slip_threshold,
        )
        if details is None:
            continue

        t0_s = cycle_offsets[event_id] + details["start_time_s"]
        if details["mw"] > min_mw:
            if (
                math.isnan(details["x_NZTM_m"])
                or math.isnan(details["y_NZTM_m"])
                or math.isnan(details["z_NZTM_m"])
            ):
                raise RuntimeError(f"No saved fault snapshot brackets the event in {q_dir}; cannot export hypocenter for File 1")
            if math.isnan(details["sbound_NZTM_m"]) or math.isnan(details["nbound_NZTM_m"]):
                raise RuntimeError(f"No patches exceeded slip threshold {rupture_bound_slip_threshold} m in {q_dir}")
            event_rows.append(
                {
                    "event_id": event_id,
                    "t0_s": t0_s,
                    "t0_year": t0_s / SECONDS_PER_YEAR,
                    "m0_Nm": details["m0_Nm"],
                    "mw": details["mw"],
                    "area_m2": details["area_m2"],
                    "dt_s": details["duration_s"],
                    "max_slip_m": details["max_slip_m"],
                    "x_NZTM_m": details["x_NZTM_m"],
                    "y_NZTM_m": details["y_NZTM_m"],
                    "z_NZTM_m": details["z_NZTM_m"],
                    "sbound_NZTM_m": details["sbound_NZTM_m"],
                    "nbound_NZTM_m": details["nbound_NZTM_m"],
                }
            )

        for site in sites:
            site_row = details["site_rows"].get(site.label)
            if site_row is None:
                continue
            site_rows_by_label[site.label].append(
                {
                    "event_id": event_id,
                    "t0_s": t0_s,
                    "Mw": details["mw"],
                    "Slip_m": site_row["Slip_m"],
                    "Rake_deg": site_row["Rake_deg"],
                }
            )

    if not event_rows:
        raise RuntimeError("No events met the export criteria")
    return patch_indices, event_rows, site_rows_by_label


def main() -> None:
    args = parse_args()

    case_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    geometry_npz = resolve_geometry_npz(case_dir, args.geometry_npz)
    prefix = args.prefix or infer_prefix(case_dir.name)

    geometry = load_geometry_model(geometry_npz)
    _, event_rows, site_rows_by_label = build_catalog_rows(
        case_dir=case_dir,
        geometry=geometry,
        sites=DEFAULT_SITES,
        min_mw=args.min_mw,
        slip_rate_threshold=args.slip_rate_threshold,
        rupture_bound_slip_threshold=args.rupture_bound_slip_threshold,
        slip_eps=args.slip_eps,
        site_slip_threshold=args.site_slip_threshold,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    event_path = output_dir / f"{prefix}_event_info.csv"
    write_csv(
        event_path,
        event_rows,
        [
            "event_id",
            "t0_s",
            "t0_year",
            "m0_Nm",
            "mw",
            "area_m2",
            "dt_s",
            "max_slip_m",
            "x_NZTM_m",
            "y_NZTM_m",
            "z_NZTM_m",
            "sbound_NZTM_m",
            "nbound_NZTM_m",
        ],
    )

    print(event_path)
    for site in DEFAULT_SITES:
        site_path = output_dir / f"{prefix}_{site.file_suffix}"
        write_csv(
            site_path,
            site_rows_by_label[site.label],
            ["event_id", "t0_s", "Mw", "Slip_m", "Rake_deg"],
        )
        print(site_path)


if __name__ == "__main__":
    main()
