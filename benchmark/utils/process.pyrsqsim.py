#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


SECONDS_PER_YEAR = 3.15e7
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = REPO_ROOT.parent / "PyRSQSim.dev/examples/alpine/results/variable_dip_1km_5km"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmark/Simulation_results/PyRSQSim_results"
DEFAULT_RAKE_DEG = 160.0
DEFAULT_SHEAR_MODULUS_PA = 30.0e9
CATALOG_SLIP_EPS = 1.0e-6

SITES = (
    ("site_1", "site1_event_info.csv", 1207801.604, 5071777.726),
    ("site_2", "site2_event_info.csv", 1382372.944, 5199926.971),
    ("site_3", "site3_event_info.csv", 1536168.087, 5311181.708),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export PyRSQSim Alpine example results as Avi-style benchmark CSV files.")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--stl", type=Path, default=None)
    parser.add_argument("--summary-npz", type=Path, default=None)
    parser.add_argument("--history-npz", type=Path, default=None)
    parser.add_argument("--catalog-csv", type=Path, default=None)
    parser.add_argument("--output-dir", "--outdir", dest="output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--prefix", type=str, default=None)
    parser.add_argument("--min-mw", type=float, default=7.0)
    parser.add_argument("--site-slip-threshold", type=float, default=0.0)
    parser.add_argument("--rupture-bound-slip-threshold", type=float, default=0.1)
    parser.add_argument("--rake-deg", type=float, default=DEFAULT_RAKE_DEG)
    parser.add_argument("--shear-modulus-pa", type=float, default=DEFAULT_SHEAR_MODULUS_PA)
    parser.add_argument("--hypocenter-policy", choices=("strict", "onset_centroid"), default="strict")
    return parser.parse_args()


def parse_stl(stl_path: Path) -> list[dict[str, np.ndarray]]:
    facets: list[dict[str, np.ndarray]] = []
    current_normal: np.ndarray | None = None
    current_vertices: list[list[float]] = []
    with stl_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("facet normal"):
                parts = line.split()
                current_normal = np.array([float(parts[2]), float(parts[3]), float(parts[4])], dtype=float)
                current_vertices = []
            elif line.startswith("vertex"):
                parts = line.split()
                current_vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith("endfacet"):
                if current_normal is None or len(current_vertices) != 3:
                    raise ValueError(f"Malformed STL facet in {stl_path}")
                facets.append({"normal": current_normal, "vertices": np.asarray(current_vertices, dtype=float)})
                current_normal = None
                current_vertices = []
    if not facets:
        raise ValueError(f"No facets parsed from {stl_path}")
    return facets


def triangle_properties(vertices: np.ndarray, normal: np.ndarray) -> tuple[np.ndarray, float]:
    v0, v1, v2 = vertices
    centroid = (v0 + v1 + v2) / 3.0
    edge1 = v1 - v0
    edge2 = v2 - v0
    cross = np.cross(edge1, edge2)
    area = 0.5 * np.linalg.norm(cross)
    if area <= 0.0:
        raise ValueError("Degenerate STL triangle encountered")
    if np.linalg.norm(normal) <= 0.0:
        raise ValueError("Zero STL normal encountered")
    return centroid, float(area)


def infer_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    results_dir = args.results_dir.resolve()
    if not results_dir.exists():
        raise FileNotFoundError(results_dir)
    stl = (args.stl or (results_dir.parents[1] / f"{results_dir.name}.stl")).resolve()
    summary = (args.summary_npz or (results_dir / "summary_data.npz")).resolve()
    history = (args.history_npz or (results_dir / "rsfsim_history.npz")).resolve()
    catalog = (args.catalog_csv or (results_dir / "earthquakes.csv")).resolve()
    for path in (stl, summary, history, catalog):
        if not path.exists():
            raise FileNotFoundError(path)
    return stl, summary, history, catalog


def infer_prefix(results_dir: Path) -> str:
    name = results_dir.name.lower()
    if "variable_dip" in name or "varying_dip" in name or "varydip" in name:
        return "alpine_varying_dip_5km"
    if "planar" in name or "dip50" in name:
        return "alpine_planar"
    return re.sub(r"[^a-zA-Z0-9]+", "_", results_dir.name).strip("_")


def load_catalog(catalog_csv: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with catalog_csv.open() as handle:
        reader = csv.DictReader(handle)
        for event_id, row in enumerate(reader):
            rows.append(
                {
                    "event_id": event_id,
                    "time_years": float(row["time_years"]),
                    "magnitude": float(row["magnitude"]),
                    "total_slip_m": float(row["total_slip_m"]),
                    "total_area_m2": float(row["total_area_m2"]),
                    "duration_s": float(row["duration_s"]),
                    "element_ids": np.fromstring(row["element_ids"], sep=" ", dtype=int),
                }
            )
    if not rows:
        raise RuntimeError(f"No catalog rows in {catalog_csv}")
    return rows


def build_segments(state: np.ndarray) -> list[tuple[int, int]]:
    rupturing = np.any(state == 2, axis=1)
    segments: list[tuple[int, int]] = []
    start_idx: int | None = None
    for i, is_rupturing in enumerate(rupturing):
        if is_rupturing and start_idx is None:
            start_idx = i
        elif (not is_rupturing) and start_idx is not None:
            segments.append((start_idx, i))
            start_idx = None
    if start_idx is not None:
        raise RuntimeError("History ends with an active rupture segment")
    return segments


def load_elements(stl_path: Path, summary_npz: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    facets = parse_stl(stl_path)
    centroids = np.zeros((len(facets), 3), dtype=float)
    areas_m2 = np.zeros(len(facets), dtype=float)
    for i, facet in enumerate(facets):
        centroid, area = triangle_properties(facet["vertices"], facet["normal"])
        centroids[i] = centroid
        areas_m2[i] = area

    summary = np.load(summary_npz, allow_pickle=True)
    summary_xyz = np.column_stack([summary["x_km"], summary["y_km"], summary["z_km"]]) * 1000.0
    if summary_xyz.shape != centroids.shape:
        raise ValueError(f"STL element count {centroids.shape[0]} does not match summary_data count {summary_xyz.shape[0]}")
    max_abs_diff = float(np.max(np.abs(summary_xyz - centroids)))
    if max_abs_diff > 1.0e-6:
        raise ValueError(f"STL centroid ordering does not match summary_data.npz (max abs diff {max_abs_diff})")
    return centroids[:, 0], centroids[:, 1], centroids[:, 2], areas_m2


def select_site_indices(x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray) -> dict[str, int]:
    surface_mask = np.isclose(z_m, np.max(z_m))
    if not np.any(surface_mask):
        raise RuntimeError("No shallowest PyRSQSim elements identified for site mapping")
    surface_ids = np.where(surface_mask)[0]
    mapping: dict[str, int] = {}
    for label, _, x_site_m, y_site_m in SITES:
        distance = np.sqrt((x_m[surface_ids] - x_site_m) ** 2 + (y_m[surface_ids] - y_site_m) ** 2 + z_m[surface_ids] ** 2)
        mapping[label] = int(surface_ids[int(np.argmin(distance))])
    return mapping


def compute_moment(area_m2: np.ndarray, slip_m: np.ndarray, shear_modulus_pa: float) -> float:
    return float(shear_modulus_pa * np.sum(area_m2 * slip_m))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def hypocenter_coordinates(
    first_rupture_ids: np.ndarray,
    x_m: np.ndarray,
    y_m: np.ndarray,
    z_m: np.ndarray,
    policy: str,
) -> tuple[float, float, float]:
    if first_rupture_ids.size == 1:
        idx = int(first_rupture_ids[0])
        return float(x_m[idx]), float(y_m[idx]), float(z_m[idx])
    if policy == "strict":
        raise RuntimeError(f"Ambiguous hypocenter onset set with {first_rupture_ids.size} elements")
    return (
        float(np.mean(x_m[first_rupture_ids])),
        float(np.mean(y_m[first_rupture_ids])),
        float(np.mean(z_m[first_rupture_ids])),
    )


def main() -> None:
    args = parse_args()
    stl_path, summary_npz, history_npz, catalog_csv = infer_paths(args)
    output_dir = args.output_dir.resolve()
    prefix = args.prefix or infer_prefix(args.results_dir.resolve())

    x_m, y_m, z_m, area_m2 = load_elements(stl_path, summary_npz)
    catalog_rows = load_catalog(catalog_csv)
    history = np.load(history_npz, allow_pickle=True)
    time_s = history["time"]
    state = history["state"]
    slip = history["slip"]
    slip_rate = history["slip_rate"]
    if slip.shape != state.shape or slip_rate.shape != state.shape:
        raise ValueError("History array shapes are inconsistent")
    if slip.shape[1] != x_m.shape[0]:
        raise ValueError("History element count does not match STL element count")

    segments = build_segments(state)
    if len(segments) != len(catalog_rows):
        raise RuntimeError(f"History rupture segments {len(segments)} do not match catalog events {len(catalog_rows)}")

    site_indices = select_site_indices(x_m, y_m, z_m)
    event_rows: list[dict[str, object]] = []
    site_rows = {label: [] for label, *_ in SITES}

    for catalog_row, (start_idx, end_idx) in zip(catalog_rows, segments):
        event_id = int(catalog_row["event_id"])
        prev_state = state[start_idx - 1] if start_idx > 0 else np.zeros(state.shape[1], dtype=state.dtype)
        if start_idx > 0:
            pre_start_slip = slip[start_idx - 1]
            if float(np.max(np.abs(slip[start_idx] - pre_start_slip))) > 1.0e-9:
                raise RuntimeError(f"Slip changed before recorded rupture onset for event {event_id}")
        event_slip = slip[end_idx] - slip[start_idx]
        if np.min(event_slip) < -1.0e-9:
            raise RuntimeError(f"Negative event slip encountered for event {event_id}")
        event_slip = np.maximum(event_slip, 0.0)

        rupture_ids = np.asarray(catalog_row["element_ids"], dtype=int)
        reconstructed_ids = np.where(event_slip > CATALOG_SLIP_EPS)[0]
        if not np.array_equal(reconstructed_ids, rupture_ids):
            raise RuntimeError(f"Catalog element_ids do not match reconstructed slipped elements for event {event_id}")
        if not math.isclose(float(time_s[end_idx] / SECONDS_PER_YEAR), float(catalog_row["time_years"]), rel_tol=0.0, abs_tol=1.0e-9):
            raise RuntimeError(f"Catalog end time does not match history segment end for event {event_id}")

        t0_s = float(time_s[start_idx])
        duration_s = float(time_s[end_idx] - time_s[start_idx])
        moment_nm = compute_moment(area_m2[rupture_ids], event_slip[rupture_ids], args.shear_modulus_pa)
        mw = float(catalog_row["magnitude"])

        first_rupture_ids = rupture_ids[(state[start_idx, rupture_ids] == 2) & (prev_state[rupture_ids] != 2)]
        if first_rupture_ids.size == 0:
            raise RuntimeError(f"No first-rupturing elements identified for event {event_id}")
        hypocenter_x_m, hypocenter_y_m, hypocenter_z_m = hypocenter_coordinates(
            first_rupture_ids,
            x_m,
            y_m,
            z_m,
            args.hypocenter_policy,
        )

        bound_ids = rupture_ids[event_slip[rupture_ids] > args.rupture_bound_slip_threshold]

        for label, _, _, _ in SITES:
            site_idx = site_indices[label]
            site_slip = float(event_slip[site_idx])
            if site_slip > args.site_slip_threshold:
                site_rows[label].append(
                    {
                        "event_id": event_id,
                        "t0_s": t0_s,
                        "Mw": mw,
                        "Slip_m": site_slip,
                        "Rake_deg": args.rake_deg,
                    }
                )

        if mw > args.min_mw:
            if bound_ids.size == 0:
                raise RuntimeError(f"No bound patches above {args.rupture_bound_slip_threshold} m for M>{args.min_mw} event {event_id}")
            event_rows.append(
                {
                    "event_id": event_id,
                    "t0_s": t0_s,
                    "t0_year": t0_s / SECONDS_PER_YEAR,
                    "m0_Nm": moment_nm,
                    "mw": mw,
                    "area_m2": float(np.sum(area_m2[rupture_ids])),
                    "dt_s": duration_s,
                    "max_slip_m": float(np.max(event_slip[rupture_ids])),
                    "x_NZTM_m": hypocenter_x_m,
                    "y_NZTM_m": hypocenter_y_m,
                    "z_NZTM_m": hypocenter_z_m,
                    "sbound_NZTM_m": float(np.min(y_m[bound_ids])),
                    "nbound_NZTM_m": float(np.max(y_m[bound_ids])),
                }
            )

    if not event_rows:
        raise RuntimeError("No events met the export criteria")

    write_csv(
        output_dir / f"{prefix}_event_info.csv",
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
    for label, file_suffix, _, _ in SITES:
        write_csv(
            output_dir / f"{prefix}_{file_suffix}",
            site_rows[label],
            ["event_id", "t0_s", "Mw", "Slip_m", "Rake_deg"],
        )

    print(output_dir / f"{prefix}_event_info.csv")
    for _, file_suffix, _, _ in SITES:
        print(output_dir / f"{prefix}_{file_suffix}")


if __name__ == "__main__":
    main()
