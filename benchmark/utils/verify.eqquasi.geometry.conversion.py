#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(BENCHMARK_ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PUBLICATION_DPI = 300
GEOMETRY_FIGSIZE = (7.4, 11.0)


def configure_matplotlib() -> None:
    matplotlib.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 16,
            "axes.titleweight": "bold",
            "axes.labelsize": 14,
            "axes.labelweight": "bold",
            "axes.linewidth": 1.8,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "xtick.major.width": 1.8,
            "ytick.major.width": 1.8,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "legend.fontsize": 11,
            "savefig.dpi": PUBLICATION_DPI,
        }
    )


def style_axis(ax) -> None:
    ax.tick_params(axis="both", which="major", width=1.8, length=6, pad=6)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        label.set_fontweight("bold")


configure_matplotlib()


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CASE_DIR = REPO_ROOT / "results/nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad"
DEFAULT_EXPORT_DIR = REPO_ROOT / "benchmark/Simulation_results/EQquasi_results"


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str


def load_process_module():
    module_path = Path(__file__).with_name("process.eqquasi.py")
    spec = importlib.util.spec_from_file_location("process_eqquasi", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify EQquasi geometry conversion and export files.")
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--geometry-npz", type=Path, default=None)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT_DIR)
    parser.add_argument("--prefix", type=str, default=None)
    parser.add_argument("--min-mw", type=float, default=7.0)
    parser.add_argument("--slip-rate-threshold", type=float, default=0.1)
    parser.add_argument("--rupture-bound-slip-threshold", type=float, default=0.5)
    parser.add_argument("--slip-eps", type=float, default=0.0)
    parser.add_argument("--site-slip-threshold", type=float, default=0.0)
    parser.add_argument("--comparison-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--figure-path", type=Path, default=None)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def load_geometry_metadata(npz_path: Path) -> dict:
    obj = np.load(npz_path, allow_pickle=True)
    return obj["metadata"].item()


def resolve_stl_path(stl_filename: str) -> Path:
    rel_path = Path(stl_filename)
    candidates = [
        rel_path,
        REPO_ROOT / rel_path,
        REPO_ROOT / "geometry" / rel_path,
        REPO_ROOT / "geometry" / "raw_geometry" / rel_path,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate

    matches = sorted((REPO_ROOT / "geometry/raw_geometry").rglob(rel_path.name))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise FileNotFoundError(f"Ambiguous STL basename {rel_path.name}: {matches}")
    raise FileNotFoundError(f"Could not find original STL {stl_filename}")


def load_stl_vertices(stl_path: Path) -> np.ndarray:
    with stl_path.open("rb") as handle:
        header = handle.read(80)
        count_bytes = handle.read(4)
        if len(count_bytes) == 4:
            triangle_count = int(np.frombuffer(count_bytes, dtype=np.uint32)[0])
            expected_size = 84 + triangle_count * 50
            if stl_path.stat().st_size == expected_size:
                data = np.fromfile(stl_path, dtype=np.dtype([("chunk", "f4", 12), ("attr", "u2")]), offset=84, count=triangle_count)
                return data["chunk"][:, 3:].reshape(-1, 3)
    vertices: list[list[float]] = []
    with stl_path.open() as handle:
        for line in handle:
            text = line.strip().split()
            if len(text) == 4 and text[0].lower() == "vertex":
                vertices.append([float(text[1]), float(text[2]), float(text[3])])
    if not vertices:
        raise ValueError(f"Could not parse STL vertices from {stl_path}")
    return np.asarray(vertices, dtype=float)


def subsample_points(points: np.ndarray, max_points: int) -> np.ndarray:
    if points.shape[0] <= max_points:
        return points
    step = int(np.ceil(points.shape[0] / max_points))
    return points[::step]


def write_geometry_figure(process_module, geometry, geometry_npz: Path, mapping: dict[str, tuple[int, int]], figure_path: Path) -> Path:
    metadata = load_geometry_metadata(geometry_npz)
    stl_path = resolve_stl_path(str(metadata["stl_filename"]))
    stl_vertices = subsample_points(load_stl_vertices(stl_path), 30000) / 1000.0

    local_station_xy = []
    global_station_xy = []
    patch_global_xy = []
    patch_local_xy = []
    for site in process_module.DEFAULT_SITES:
        site_x_km, site_y_km = geometry.global_to_local(site.x_nztm_m, site.y_nztm_m)
        local_station_xy.append((site_x_km, site_y_km))
        global_station_xy.append((site.x_nztm_m / 1000.0, site.y_nztm_m / 1000.0))
        patch_idx = mapping[site.label]
        patch_global_xy.append((geometry.global_x_m[patch_idx] / 1000.0, geometry.global_y_m[patch_idx] / 1000.0))
        patch_local_xy.append((geometry.local_x_km[patch_idx], geometry.local_y_km[patch_idx]))

    fig, axes = plt.subplots(2, 1, figsize=GEOMETRY_FIGSIZE)

    axes[0].scatter(stl_vertices[:, 0], stl_vertices[:, 1], s=0.4, c="0.75", linewidths=0, rasterized=True)
    for site, station_xy, patch_xy in zip(process_module.DEFAULT_SITES, global_station_xy, patch_global_xy):
        axes[0].scatter(station_xy[0], station_xy[1], marker="*", s=340, c="tab:red", edgecolors="black", linewidths=1.2, zorder=3)
        axes[0].scatter(patch_xy[0], patch_xy[1], marker="o", s=130, facecolors="none", edgecolors="black", linewidths=1.6, zorder=4)
        axes[0].text(station_xy[0], station_xy[1], f" {site.label}", fontsize=11, fontweight="bold", va="bottom")
    axes[0].set_title(f"Original STL in UTM\n{stl_path.name}")
    axes[0].set_xlabel("Easting (km)")
    axes[0].set_ylabel("Northing (km)")
    axes[0].set_aspect("equal")
    style_axis(axes[0])

    processed_xy = np.column_stack([geometry.local_x_km.ravel(), geometry.local_y_km.ravel()])
    axes[1].scatter(processed_xy[:, 0], processed_xy[:, 1], s=4, c="0.75", linewidths=0, rasterized=True)
    for site, station_xy, patch_xy in zip(process_module.DEFAULT_SITES, local_station_xy, patch_local_xy):
        axes[1].scatter(station_xy[0], station_xy[1], marker="*", s=340, c="tab:blue", edgecolors="black", linewidths=1.2, zorder=3)
        axes[1].scatter(patch_xy[0], patch_xy[1], marker="o", s=130, facecolors="none", edgecolors="black", linewidths=1.6, zorder=4)
        axes[1].plot([station_xy[0], patch_xy[0]], [station_xy[1], patch_xy[1]], c="black", lw=1.8, alpha=0.65, zorder=2)
        axes[1].text(station_xy[0], station_xy[1], f" {site.label}", fontsize=11, fontweight="bold", va="bottom")
    axes[1].set_title("Processed EQquasi Geometry")
    axes[1].set_xlabel("Local strike x (km)")
    axes[1].set_ylabel("Local fault-normal y (km)")
    axes[1].set_aspect("equal")
    style_axis(axes[1])

    plt.tight_layout()
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(figure_path, dpi=PUBLICATION_DPI, bbox_inches="tight")
    plt.close(fig)
    return figure_path


def check_round_trip(process_module, geometry) -> CheckResult:
    xy_km = np.column_stack([geometry.global_x_m.ravel(), geometry.global_y_m.ravel()]) / 1000.0
    rotated_xy_km = xy_km @ geometry.rotation_matrix.T
    local_xy_km = rotated_xy_km - geometry.center_xy_km
    x_local = local_xy_km[:, 0]
    y_local = local_xy_km[:, 1]
    z_local = geometry.global_z_m.ravel() / 1000.0 - geometry.top_z_km
    err_x = float(np.max(np.abs(x_local.reshape(geometry.local_x_km.shape) - geometry.local_x_km)))
    err_y = float(np.max(np.abs(y_local.reshape(geometry.local_y_km.shape) - geometry.local_y_km)))
    err_z = float(np.max(np.abs(z_local.reshape(geometry.local_z_km.shape) - geometry.local_z_km)))
    ok = err_x < 1.0e-9 and err_y < 1.0e-9 and err_z < 1.0e-9
    return CheckResult("geometry_round_trip", ok, f"max_abs_error_km: x={err_x:.3e}, y={err_y:.3e}, z={err_z:.3e}")


def check_site_mapping(process_module, geometry) -> tuple[CheckResult, dict[str, tuple[int, int]]]:
    mapping = process_module.find_nearest_patch_indices(geometry, process_module.DEFAULT_SITES)
    details = []
    ok = True
    nx = geometry.local_x_km.shape[1]
    surface_row = process_module.surface_row_index(geometry)
    for site in process_module.DEFAULT_SITES:
        site_x_km, site_y_km = geometry.global_to_local(site.x_nztm_m, site.y_nztm_m)
        iz, ix = mapping[site.label]
        dx = float(geometry.local_x_km[iz, ix] - site_x_km)
        dy = float(geometry.local_y_km[iz, ix] - site_y_km)
        dz = float(geometry.local_z_km[iz, ix])
        distance_km = float(np.sqrt(dx**2 + dy**2 + dz**2))
        pct = 100.0 * ix / (nx - 1)
        details.append(f"{site.label}: idx=({iz},{ix}) pct={pct:.2f} distance_km={distance_km:.3f}")
        ok = ok and distance_km <= 2.0 * np.sqrt(2.0) and iz == surface_row
    return CheckResult("station_patch_mapping", ok, " | ".join(details)), mapping


def build_expected_rows(process_module, case_dir: Path, geometry, min_mw: float, slip_rate_threshold: float, rupture_bound_slip_threshold: float, slip_eps: float, site_slip_threshold: float) -> tuple[list[dict], dict[str, list[dict]]]:
    _, event_rows, site_rows_by_label = process_module.build_catalog_rows(
        case_dir=case_dir,
        geometry=geometry,
        sites=process_module.DEFAULT_SITES,
        min_mw=min_mw,
        slip_rate_threshold=slip_rate_threshold,
        rupture_bound_slip_threshold=rupture_bound_slip_threshold,
        slip_eps=slip_eps,
        site_slip_threshold=site_slip_threshold,
    )
    return event_rows, site_rows_by_label


def compare_keyed_rows(
    expected_rows: list[dict],
    exported_rows: list[dict[str, str]],
    numeric_fields: tuple[str, ...],
    tolerance: float,
) -> tuple[bool, str]:
    expected = {int(row["event_id"]): row for row in expected_rows}
    exported = {int(row["event_id"]): row for row in exported_rows}
    missing_ids = sorted(set(expected) - set(exported))
    extra_ids = sorted(set(exported) - set(expected))
    if missing_ids or extra_ids:
        return False, f"missing={missing_ids} extra={extra_ids}"

    max_delta = 0.0
    max_field = ""
    for event_id, expected_row in expected.items():
        exported_row = exported[event_id]
        for field in numeric_fields:
            delta = abs(float(expected_row[field]) - float(exported_row[field]))
            if delta > max_delta:
                max_delta = delta
                max_field = f"{event_id}:{field}"
            if delta > tolerance:
                return False, f"event_id={event_id} field={field} delta={delta:.3e}"
    return True, f"rows={len(expected_rows)} max_delta={max_delta:.3e} at {max_field or 'n/a'}"


def check_site_export_consistency(process_module, export_dir: Path, prefix: str, expected_site_rows: dict[str, list[dict]], tolerance: float) -> CheckResult:
    details = []
    ok = True
    for site in process_module.DEFAULT_SITES:
        site_path = export_dir / f"{prefix}_{site.file_suffix}"
        if not site_path.exists():
            return CheckResult("site_export_consistency", False, f"missing export file: {site_path}")
        match_ok, detail = compare_keyed_rows(
            expected_rows=expected_site_rows[site.label],
            exported_rows=read_csv_rows(site_path),
            numeric_fields=("event_id", "t0_s", "Mw", "Slip_m", "Rake_deg"),
            tolerance=tolerance,
        )
        ok = ok and match_ok
        details.append(f"{site.label}: {detail}")
    return CheckResult("site_export_consistency", ok, " | ".join(details))


def check_event_summary_consistency(export_dir: Path, prefix: str, expected_event_rows: list[dict], tolerance: float) -> CheckResult:
    export_path = export_dir / f"{prefix}_event_info.csv"
    if not export_path.exists():
        return CheckResult("event_summary_consistency", False, f"missing export file: {export_path}")
    ok, detail = compare_keyed_rows(
        expected_rows=expected_event_rows,
        exported_rows=read_csv_rows(export_path),
        numeric_fields=(
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
        ),
        tolerance=tolerance,
    )
    return CheckResult("event_summary_consistency", ok, detail)


def main() -> int:
    args = parse_args()
    process_module = load_process_module()
    case_dir = args.case_dir.resolve()
    geometry_npz = process_module.resolve_geometry_npz(case_dir, args.geometry_npz).resolve()
    geometry = process_module.load_geometry_model(geometry_npz)
    export_dir = args.export_dir.resolve()
    prefix = args.prefix or process_module.infer_prefix(case_dir.name)
    figure_path = (args.figure_path.resolve() if args.figure_path is not None else export_dir / f"{prefix}_geometry_conversion_check.png")
    expected_event_rows, expected_site_rows = build_expected_rows(
        process_module,
        case_dir,
        geometry,
        args.min_mw,
        args.slip_rate_threshold,
        args.rupture_bound_slip_threshold,
        args.slip_eps,
        args.site_slip_threshold,
    )

    results: list[CheckResult] = []
    mapping_check, mapping = check_site_mapping(process_module, geometry)
    results.append(check_round_trip(process_module, geometry))
    results.append(mapping_check)
    results.append(
        check_site_export_consistency(
            process_module,
            export_dir,
            prefix,
            expected_site_rows,
            args.comparison_tolerance,
        )
    )
    results.append(
        check_event_summary_consistency(
            export_dir,
            prefix,
            expected_event_rows,
            args.comparison_tolerance,
        )
    )
    written_figure = write_geometry_figure(process_module, geometry, geometry_npz, mapping, figure_path)
    results.append(CheckResult("geometry_visualization", True, str(written_figure)))

    failed = False
    for result in results:
        print(f"[{'PASS' if result.ok else 'FAIL'}] {result.name}: {result.detail}")
        failed = failed or (not result.ok)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
