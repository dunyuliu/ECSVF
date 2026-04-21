#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import importlib.util
import itertools
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = BENCHMARK_ROOT.parent
os.environ.setdefault("MPLCONFIGDIR", str(BENCHMARK_ROOT / ".mplconfig"))


def _load_local_module(stem: str):
    """Import a module from the same directory as this script by filename stem.

    The stem may contain dots (e.g. 'plot.accumulated.slip.eqquasi'); they are
    replaced with underscores for the internal module name so that importlib
    does not treat them as package separators.
    """
    path = Path(__file__).with_name(stem + ".py")
    mod_name = stem.replace(".", "_")
    spec = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod  # register before exec so @dataclass __module__ lookup works
    spec.loader.exec_module(mod)
    return mod


def _find_eqquasi_case_dir(scenario: str) -> Path | None:
    """Search REPO_ROOT/results for the EQquasi case directory matching scenario."""
    results_dir = REPO_ROOT / "results"
    if not results_dir.exists():
        return None
    # More-specific tokens listed first so they win when multiple dirs match.
    # Order mirrors the hints in process.eqquasi.py::resolve_geometry_npz.
    tokens = {
        "planar": ["dip50", "planar"],
        "varying_dip": ["varyDip20251202", "varyDip", "varying_dip"],
        "no_dip_change": ["noDipChange", "no_dip_change"],
    }
    for tok in tokens.get(scenario, []):
        for candidate in sorted(results_dir.iterdir()):
            if candidate.is_dir() and tok in candidate.name:
                return candidate
    return None

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PUBLICATION_DPI = 300
FIGURE_SIZE = (7.2, 5.2)
LINE_WIDTH = 2.8


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
            "xtick.minor.width": 1.2,
            "ytick.minor.width": 1.2,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "legend.fontsize": 11,
            "lines.linewidth": LINE_WIDTH,
            "savefig.dpi": PUBLICATION_DPI,
        }
    )


def style_axis(ax) -> None:
    ax.tick_params(axis="both", which="major", width=1.8, length=6, pad=6)
    ax.tick_params(axis="both", which="minor", width=1.2, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        label.set_fontweight("bold")


def style_legend(ax, loc: str = "best"):
    legend = ax.legend(frameon=False, loc=loc)
    if legend is None:
        return None
    for text in legend.get_texts():
        text.set_fontweight("bold")
    return legend


configure_matplotlib()


EVENT_NUMERIC_FIELDS = (
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
)
SITE_NUMERIC_FIELDS = ("event_id", "t0_s", "Mw", "Slip_m", "Rake_deg")


@dataclass(frozen=True)
class Bundle:
    source: str
    bundle_id: str
    scenario: str
    event_file: Path
    site_files: dict[int, Path]

    @property
    def label(self) -> str:
        return f"{self.source}:{self.bundle_id}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare complete Avi-style benchmark bundles.")
    parser.add_argument("--simulation-root", type=Path, default=BENCHMARK_ROOT / "Simulation_results")
    parser.add_argument("--outdir", type=Path, default=BENCHMARK_ROOT / "analytics")
    return parser.parse_args()


def infer_scenario(name: str) -> str:
    lowered = name.lower()
    if "planar" in lowered:
        return "planar"
    if "varying_dip" in lowered or "varydip" in lowered:
        return "varying_dip"
    if "no_dip_change" in lowered or "nodipchange" in lowered:
        return "no_dip_change"
    return "unknown"


def discover_bundles(simulation_root: Path) -> list[Bundle]:
    if not simulation_root.exists():
        raise FileNotFoundError(simulation_root)

    bundles: list[Bundle] = []
    for result_dir in sorted(simulation_root.glob("*_results")):
        source = result_dir.name.replace("_results", "")
        event_files = [
            path
            for path in sorted(result_dir.glob("*_event_info.csv"))
            if re.fullmatch(r".+_site[123]_event_info\.csv", path.name) is None
        ]
        for event_file in event_files:
            bundle_id = event_file.stem[: -len("_event_info")]
            scenario = infer_scenario(bundle_id)
            if scenario == "unknown":
                raise ValueError(f"Unrecognized benchmark scenario name: {bundle_id}")
            site_files = {
                1: result_dir / f"{bundle_id}_site1_event_info.csv",
                2: result_dir / f"{bundle_id}_site2_event_info.csv",
                3: result_dir / f"{bundle_id}_site3_event_info.csv",
            }
            missing = [str(path) for path in site_files.values() if not path.exists()]
            if missing:
                raise FileNotFoundError(f"Incomplete bundle for {event_file}: {missing}")
            bundles.append(
                Bundle(
                    source=source,
                    bundle_id=bundle_id,
                    scenario=scenario,
                    event_file=event_file,
                    site_files=site_files,
                )
            )
    if not bundles:
        raise RuntimeError(f"No complete bundles found under {simulation_root}")
    return bundles


def load_rows(path: Path, numeric_fields: tuple[str, ...]) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parsed: dict[str, float] = {}
            for key, value in row.items():
                parsed[key] = float(value) if key in numeric_fields else value
            rows.append(parsed)
    return rows


def safe_mean(values: np.ndarray) -> float:
    return float(np.mean(values)) if values.size else float("nan")


def safe_std(values: np.ndarray) -> float:
    return float(np.std(values, ddof=1)) if values.size >= 2 else float("nan")


def safe_quantile(values: np.ndarray, q: float) -> float:
    return float(np.quantile(values, q)) if values.size else float("nan")


def empirical_cdf_distance(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 or b.size == 0:
        return float("nan")
    values = np.unique(np.concatenate([a, b]))
    cdf_a = np.searchsorted(np.sort(a), values, side="right") / a.size
    cdf_b = np.searchsorted(np.sort(b), values, side="right") / b.size
    return float(np.max(np.abs(cdf_a - cdf_b)))


def circular_mean_deg(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    angles = np.deg2rad(values)
    mean_vector = np.mean(np.exp(1j * angles))
    if np.isclose(abs(mean_vector), 0.0):
        return float("nan")
    return float((np.rad2deg(np.angle(mean_vector)) + 360.0) % 360.0)


def circular_delta_deg(a_deg: float, b_deg: float) -> float:
    if math.isnan(a_deg) or math.isnan(b_deg):
        return float("nan")
    return float(abs((a_deg - b_deg + 180.0) % 360.0 - 180.0))


def bundle_summary_rows(bundles: list[Bundle]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for bundle in bundles:
        events = load_rows(bundle.event_file, EVENT_NUMERIC_FIELDS)
        event_years = np.array([row["t0_year"] for row in events], dtype=float)
        recurrence = np.diff(np.sort(event_years)) if event_years.size >= 2 else np.array([], dtype=float)
        row: dict[str, object] = {
            "source": bundle.source,
            "bundle_id": bundle.bundle_id,
            "scenario": bundle.scenario,
            "n_events": len(events),
            "mw_mean": safe_mean(np.array([row["mw"] for row in events], dtype=float)),
            "mw_max": float(np.max([row["mw"] for row in events])) if events else float("nan"),
            "duration_mean_s": safe_mean(np.array([row["dt_s"] for row in events], dtype=float)),
            "area_mean_m2": safe_mean(np.array([row["area_m2"] for row in events], dtype=float)),
            "recurrence_mean_years": safe_mean(recurrence),
            "recurrence_std_years": safe_std(recurrence),
        }
        for site in (1, 2, 3):
            site_rows = load_rows(bundle.site_files[site], SITE_NUMERIC_FIELDS)
            slips = np.array([item["Slip_m"] for item in site_rows], dtype=float)
            row[f"site{site}_count"] = len(site_rows)
            row[f"site{site}_slip_mean_m"] = safe_mean(slips)
            row[f"site{site}_slip_p95_m"] = safe_quantile(slips, 0.95)
        rows.append(row)
    return rows


def pairwise_event_rows(bundles: list[Bundle]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    grouped: dict[str, list[Bundle]] = {}
    for bundle in bundles:
        grouped.setdefault(bundle.scenario, []).append(bundle)
    for scenario_bundles in grouped.values():
        for lhs, rhs in itertools.combinations(sorted(scenario_bundles, key=lambda item: item.label), 2):
            lhs_rows = load_rows(lhs.event_file, EVENT_NUMERIC_FIELDS)
            rhs_rows = load_rows(rhs.event_file, EVENT_NUMERIC_FIELDS)
            lhs_mw = np.array([row["mw"] for row in lhs_rows], dtype=float)
            rhs_mw = np.array([row["mw"] for row in rhs_rows], dtype=float)
            lhs_recur = np.diff(np.sort(np.array([row["t0_year"] for row in lhs_rows], dtype=float))) if len(lhs_rows) >= 2 else np.array([], dtype=float)
            rhs_recur = np.diff(np.sort(np.array([row["t0_year"] for row in rhs_rows], dtype=float))) if len(rhs_rows) >= 2 else np.array([], dtype=float)
            rows.append(
                {
                    "scenario": lhs.scenario,
                    "lhs": lhs.label,
                    "rhs": rhs.label,
                    "lhs_n_events": len(lhs_rows),
                    "rhs_n_events": len(rhs_rows),
                    "delta_n_events": len(lhs_rows) - len(rhs_rows),
                    "lhs_mw_mean": safe_mean(lhs_mw),
                    "rhs_mw_mean": safe_mean(rhs_mw),
                    "delta_mw_mean": safe_mean(lhs_mw) - safe_mean(rhs_mw),
                    "mw_cdf_distance": empirical_cdf_distance(lhs_mw, rhs_mw),
                    "lhs_recurrence_mean_years": safe_mean(lhs_recur),
                    "rhs_recurrence_mean_years": safe_mean(rhs_recur),
                    "delta_recurrence_mean_years": safe_mean(lhs_recur) - safe_mean(rhs_recur),
                    "recurrence_cdf_distance": empirical_cdf_distance(lhs_recur, rhs_recur),
                }
            )
    return rows


def pairwise_site_rows(bundles: list[Bundle]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    grouped: dict[str, list[Bundle]] = {}
    for bundle in bundles:
        grouped.setdefault(bundle.scenario, []).append(bundle)
    for scenario_bundles in grouped.values():
        for lhs, rhs in itertools.combinations(sorted(scenario_bundles, key=lambda item: item.label), 2):
            for site in (1, 2, 3):
                lhs_rows = load_rows(lhs.site_files[site], SITE_NUMERIC_FIELDS)
                rhs_rows = load_rows(rhs.site_files[site], SITE_NUMERIC_FIELDS)
                lhs_slip = np.array([row["Slip_m"] for row in lhs_rows], dtype=float)
                rhs_slip = np.array([row["Slip_m"] for row in rhs_rows], dtype=float)
                lhs_rake = np.array([row["Rake_deg"] for row in lhs_rows], dtype=float)
                rhs_rake = np.array([row["Rake_deg"] for row in rhs_rows], dtype=float)
                lhs_rake_mean = circular_mean_deg(lhs_rake)
                rhs_rake_mean = circular_mean_deg(rhs_rake)
                rows.append(
                    {
                        "scenario": lhs.scenario,
                        "site": site,
                        "lhs": lhs.label,
                        "rhs": rhs.label,
                        "lhs_count": len(lhs_rows),
                        "rhs_count": len(rhs_rows),
                        "lhs_slip_mean_m": safe_mean(lhs_slip),
                        "rhs_slip_mean_m": safe_mean(rhs_slip),
                        "delta_slip_mean_m": safe_mean(lhs_slip) - safe_mean(rhs_slip),
                        "lhs_slip_p95_m": safe_quantile(lhs_slip, 0.95),
                        "rhs_slip_p95_m": safe_quantile(rhs_slip, 0.95),
                        "slip_cdf_distance": empirical_cdf_distance(lhs_slip, rhs_slip),
                        "lhs_rake_mean_deg": lhs_rake_mean,
                        "rhs_rake_mean_deg": rhs_rake_mean,
                        "rake_circular_delta_deg": circular_delta_deg(lhs_rake_mean, rhs_rake_mean),
                    }
                )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.unlink(missing_ok=True)
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def remove_if_unwritten(path: Path, written: bool) -> None:
    if not written:
        path.unlink(missing_ok=True)


def plot_event_mw_histogram(bundles: list[Bundle], outdir: Path) -> None:
    output_path = outdir / "event_mw_histogram.png"
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plotted = False
    bundle_series: list[tuple[Bundle, np.ndarray]] = []
    for bundle in sorted(bundles, key=lambda item: item.label):
        rows = load_rows(bundle.event_file, EVENT_NUMERIC_FIELDS)
        values = np.array([row["mw"] for row in rows], dtype=float)
        if values.size == 0:
            continue
        bundle_series.append((bundle, values))

    if bundle_series:
        all_values = np.concatenate([values for _, values in bundle_series])
        bin_min = math.floor(all_values.min() * 10.0) / 10.0
        bin_max = math.ceil(all_values.max() * 10.0) / 10.0
        if math.isclose(bin_min, bin_max):
            bin_max = bin_min + 0.1
        bins = np.arange(bin_min, bin_max + 0.11, 0.1)
        line_styles = ["-", "--", "-.", ":"]

        for idx, (bundle, values) in enumerate(bundle_series):
            density, edges = np.histogram(values, bins=bins, density=True)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.plot(
                centers,
                density,
                linestyle=line_styles[idx % len(line_styles)],
                linewidth=LINE_WIDTH,
                label=f"{bundle.label} (n={values.size})",
            )
            plotted = True

        ax.set_title("Event Magnitude Distribution")
        ax.set_xlabel("Mw")
        ax.set_ylabel("Density")
        ax.grid(True, alpha=0.3)
        style_axis(ax)
        style_legend(ax)
    if plotted:
        fig.tight_layout()
        fig.savefig(output_path, dpi=PUBLICATION_DPI)
    plt.close(fig)
    remove_if_unwritten(output_path, plotted)


def plot_recurrence_cdf(bundles: list[Bundle], outdir: Path) -> None:
    output_path = outdir / "recurrence_cdf.png"
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plotted = False
    for bundle in bundles:
        rows = load_rows(bundle.event_file, EVENT_NUMERIC_FIELDS)
        years = np.sort(np.array([row["t0_year"] for row in rows], dtype=float))
        recurrence = np.diff(years) if years.size >= 2 else np.array([], dtype=float)
        if recurrence.size == 0:
            continue
        x = np.sort(recurrence)
        y = np.arange(1, len(x) + 1) / len(x)
        ax.step(x, y, where="post", label=bundle.label, linewidth=LINE_WIDTH)
        plotted = True
    if plotted:
        ax.set_xlabel("Recurrence Interval (years)")
        ax.set_ylabel("Empirical CDF")
        ax.set_title("Recurrence Interval Comparison")
        ax.grid(True, alpha=0.3)
        style_axis(ax)
        style_legend(ax)
        fig.tight_layout()
        fig.savefig(output_path, dpi=PUBLICATION_DPI)
    plt.close(fig)
    remove_if_unwritten(output_path, plotted)


def plot_site_slip_cdfs(bundles: list[Bundle], outdir: Path, label_to_color: dict[str, str]) -> None:
    for site in (1, 2, 3):
        output_path = outdir / f"site{site}_slip_cdf.png"
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        plotted = False
        for bundle in bundles:
            rows = load_rows(bundle.site_files[site], SITE_NUMERIC_FIELDS)
            slips = np.sort(np.array([row["Slip_m"] for row in rows], dtype=float))
            if slips.size == 0:
                continue
            y = np.arange(1, len(slips) + 1) / len(slips)
            ax.step(slips, y, where="post", label=bundle.label, linewidth=LINE_WIDTH,
                    color=label_to_color[bundle.label])
            plotted = True
        if plotted:
            ax.set_xlabel("Slip (m)")
            ax.set_ylabel("Empirical CDF")
            ax.set_title(f"Site {site} Slip Comparison")
            ax.grid(True, alpha=0.3)
            style_axis(ax)
            style_legend(ax)
            fig.tight_layout()
            fig.savefig(output_path, dpi=PUBLICATION_DPI)
        plt.close(fig)
        remove_if_unwritten(output_path, plotted)


def plot_rupture_extents(bundles: list[Bundle], outdir: Path, label_to_color: dict[str, str]) -> None:
    output_path = outdir / "rupture_extents.png"
    SOURCES = {"EQquasi", "RSQSim", "MCQsim", "HBI"}
    WINDOW_YR = 3000.0

    # One bundle per source: prefer varying_dip, else first available
    source_bundle: dict[str, Bundle] = {}
    for bundle in bundles:
        if bundle.source not in SOURCES:
            continue
        existing = source_bundle.get(bundle.source)
        if existing is None or (bundle.scenario == "varying_dip" and existing.scenario != "varying_dip"):
            source_bundle[bundle.source] = bundle

    if not source_bundle:
        return

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    plotted = False

    for source in sorted(source_bundle):
        bundle = source_bundle[source]
        rows = load_rows(bundle.event_file, EVENT_NUMERIC_FIELDS)
        sbounds = np.array([row["sbound_NZTM_m"] for row in rows], dtype=float)
        nbounds = np.array([row["nbound_NZTM_m"] for row in rows], dtype=float)
        years = np.array([row["t0_year"] for row in rows], dtype=float)
        valid = np.isfinite(sbounds) & np.isfinite(nbounds) & np.isfinite(years)
        if not valid.any():
            continue
        sbounds = sbounds[valid]
        nbounds = nbounds[valid]
        years = years[valid]

        t_start = years.min()
        window = years <= t_start + WINDOW_YR
        sbounds_km = sbounds[window] / 1e3
        nbounds_km = nbounds[window] / 1e3
        years_rel = years[window] - t_start

        color = label_to_color[bundle.label]
        n_shown = int(window.sum())
        alpha = max(0.25, min(0.85, 30.0 / max(n_shown, 1)))
        lw = max(0.6, min(LINE_WIDTH, 60.0 / max(n_shown, 1)))
        label_added = False
        for s, n, t in zip(sbounds_km, nbounds_km, years_rel):
            ax.plot(
                [s, n], [t, t],
                color=color, linewidth=lw, alpha=alpha,
                label=f"{bundle.label} (n={n_shown})" if not label_added else "_nolegend_",
                solid_capstyle="butt",
            )
            label_added = True
        plotted = True

    if plotted:
        ax.set_xlabel("Along-Strike Northing (km NZTM)")
        ax.set_ylabel("Time (years)")
        ax.set_title("Rupture Extents — varying dip (3000-yr window)")
        ax.grid(True, alpha=0.3)
        style_axis(ax)
        style_legend(ax, loc="upper left")
        fig.tight_layout()
        fig.savefig(output_path, dpi=PUBLICATION_DPI)
    plt.close(fig)
    remove_if_unwritten(output_path, plotted)


def plot_long_term_slip_rate(bundles: list[Bundle], outdir: Path, label_to_color: dict[str, str]) -> None:
    output_path = outdir / "long_term_slip_rate.png"
    BIN_WIDTH_M = 10e3  # 10 km bins along strike

    fig, ax = plt.subplots(figsize=(8, 6))
    plotted = False

    for bundle in bundles:
        color = label_to_color[bundle.label]

        # For EQquasi: use actual spatial slip profiles via plot.accumulated.slip.eqquasi functions
        if bundle.source == "EQquasi":
            case_dir = _find_eqquasi_case_dir(bundle.scenario)
            plot_mod = _load_local_module("plot.accumulated.slip.eqquasi")

            q_dirs = plot_mod.find_q_dirs(str(case_dir))
            end_times_yr = plot_mod.cumulative_end_times_yr(q_dirs)

            # Accumulate profiles; determine ncols from first valid profile
            accumulated = None
            total_yr = 0.0
            for qd, t_yr in zip(q_dirs, end_times_yr):
                if t_yr is None:
                    continue
                profile = plot_mod.load_slip_profile(qd)
                if profile is None:
                    continue
                if accumulated is None:
                    accumulated = np.zeros(profile.size)
                accumulated += profile
                total_yr = t_yr

            if accumulated is None or total_yr <= 0:
                continue

            # Map columns linearly to NZTM northing using event CSV fault extent
            event_rows_eq = load_rows(bundle.event_file, EVENT_NUMERIC_FIELDS)
            sbounds_eq = [r["sbound_NZTM_m"] for r in event_rows_eq if np.isfinite(r["sbound_NZTM_m"])]
            nbounds_eq = [r["nbound_NZTM_m"] for r in event_rows_eq if np.isfinite(r["nbound_NZTM_m"])]
            northing_m = np.linspace(min(sbounds_eq), max(nbounds_eq), accumulated.size)

            slip_rate_mm_yr = accumulated / total_yr * 1e3
            ax.plot(northing_m / 1e3, slip_rate_mm_yr, linewidth=LINE_WIDTH,
                    color=color, label=bundle.label)
            plotted = True
            continue

        # All other models: uniform 0.5 * max_slip over rupture extent
        rows = load_rows(bundle.event_file, EVENT_NUMERIC_FIELDS)
        sbounds = np.array([row["sbound_NZTM_m"] for row in rows], dtype=float)
        nbounds = np.array([row["nbound_NZTM_m"] for row in rows], dtype=float)
        max_slips = np.array([row["max_slip_m"] for row in rows], dtype=float)
        years = np.array([row["t0_year"] for row in rows], dtype=float)
        valid = np.isfinite(sbounds) & np.isfinite(nbounds) & np.isfinite(max_slips) & np.isfinite(years)
        if not valid.any():
            continue
        sbounds, nbounds, max_slips, years = sbounds[valid], nbounds[valid], max_slips[valid], years[valid]
        time_span_yr = years.max() - years.min()
        if time_span_yr <= 0:
            continue
        bin_edges = np.arange(sbounds.min(), nbounds.max() + BIN_WIDTH_M, BIN_WIDTH_M)
        if bin_edges.size < 2:
            continue
        bin_centers_km = 0.5 * (bin_edges[:-1] + bin_edges[1:]) / 1e3
        slip_accum = np.zeros(len(bin_centers_km))
        for s, n, slip in zip(sbounds, nbounds, max_slips):
            in_rupture = (bin_edges[1:] > s) & (bin_edges[:-1] < n)
            slip_accum[in_rupture] += 0.5 * slip
        slip_rate_mm_yr = slip_accum / time_span_yr * 1e3
        ax.plot(bin_centers_km, slip_rate_mm_yr, linewidth=LINE_WIDTH,
                color=color, label=bundle.label)
        plotted = True

    if plotted:
        ax.set_xlabel("Along-Strike Northing (km NZTM)")
        ax.set_ylabel("Slip Rate (mm/yr)")
        ax.set_title("Long-Term Coseismic Slip Rate")
        ax.grid(True, alpha=0.3)
        style_axis(ax)
        style_legend(ax)
        fig.tight_layout()
        fig.savefig(output_path, dpi=PUBLICATION_DPI)
    plt.close(fig)
    remove_if_unwritten(output_path, plotted)


def write_overview(path: Path, bundles: list[Bundle]) -> None:
    lines = ["Benchmark comparison analytics", ""]
    for bundle in bundles:
        lines.append(f"- {bundle.label} | scenario={bundle.scenario} | event_file={bundle.event_file}")
    lines.append("")
    lines.append("Generated files:")
    generated = [
        "bundle_summary.csv",
        "pairwise_event_comparison.csv",
        "pairwise_site_comparison.csv",
        "event_mw_histogram.png",
        "recurrence_cdf.png",
        "site1_slip_cdf.png",
        "site2_slip_cdf.png",
        "site3_slip_cdf.png",
        "rupture_extents.png",
        "long_term_slip_rate.png",
    ]
    for name in generated:
        if (path.parent / name).exists():
            lines.append(f"- {name}")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    bundles = discover_bundles(args.simulation_root.resolve())
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    label_to_color = {bundle.label: palette[idx % len(palette)] for idx, bundle in enumerate(bundles)}

    summary_rows = bundle_summary_rows(bundles)
    event_rows = pairwise_event_rows(bundles)
    site_rows = pairwise_site_rows(bundles)

    write_csv(outdir / "bundle_summary.csv", summary_rows)
    write_csv(outdir / "pairwise_event_comparison.csv", event_rows)
    write_csv(outdir / "pairwise_site_comparison.csv", site_rows)
    plot_event_mw_histogram(bundles, outdir)
    plot_recurrence_cdf(bundles, outdir)
    plot_site_slip_cdfs(bundles, outdir, label_to_color)
    plot_rupture_extents(bundles, outdir, label_to_color)
    plot_long_term_slip_rate(bundles, outdir, label_to_color)
    write_overview(outdir / "README.txt", bundles)

    print(f"{len(bundles)} bundles")
    print(outdir)


if __name__ == "__main__":
    main()
