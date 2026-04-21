#!/usr/bin/env python3
"""
Plot accumulated coseismic slip and long-term slip rate for an EQquasi case.

Produces three separate figures:
  1. Individual earthquake slip profiles (colored by time)
  2. Running accumulated coseismic slip
  3. Long-term coseismic slip rate along fault strike

Each Q folder's fault.r.nc contains the coseismic slip distribution for that
earthquake. Slip is extracted at DEPTH_PCT (88%) depth, consistent with other
post-processing scripts in this project.

Usage:
    python benchmark/utils/plot.accumulated.slip.eqquasi.py --case-dir <path/to/case>

Example:
    python benchmark/utils/plot.accumulated.slip.eqquasi.py \
        --case-dir results/nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad
"""

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

FIGURE_SIZE = (8, 6)
PUBLICATION_DPI = 300
LINE_WIDTH = 2.5

matplotlib.rcParams.update({
    "font.size": 14,
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
    "legend.fontsize": 12,
    "lines.linewidth": LINE_WIDTH,
    "savefig.dpi": PUBLICATION_DPI,
    "font.weight": "bold",
})

GRID_SPACING_M = 2000.0
SEC_PER_YEAR = 365.25 * 24 * 3600
DEPTH_PCT = 88.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot accumulated EQquasi coseismic slip along fault strike.")
    parser.add_argument("--case-dir", type=str, required=True,
                        help="Path to EQquasi case directory containing Q* folders")
    parser.add_argument("--outdir", type=str, default=None,
                        help="Output directory for figures (default: <case-dir>)")
    return parser.parse_args()


def natural_sort_key(s: str) -> list:
    return [int(c) if c.isdigit() else c for c in re.split(r"(\d+)", s)]


def style_axis(ax) -> None:
    ax.tick_params(axis="both", which="major", width=1.8, length=6, pad=6)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        label.set_fontweight("bold")


def find_q_dirs(case_dir: str) -> list[str]:
    dirs = sorted(
        [d for d in os.listdir(case_dir) if re.fullmatch(r"Q\d+", d)],
        key=natural_sort_key,
    )
    return [os.path.join(case_dir, d) for d in dirs]


def cumulative_end_times_yr(q_dir_list: list[str]) -> list[float | None]:
    total_s = 0.0
    times: list[float | None] = []
    for qd in q_dir_list:
        gdat = os.path.join(qd, "global.dat")
        if not os.path.exists(gdat):
            times.append(None)
            continue
        data = np.loadtxt(gdat)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        total_s += float(data[-1, 0])
        times.append(total_s / SEC_PER_YEAR)
    return times


SLIP_THRESHOLD_M = 0.5  # consistent with process.eqquasi.py --rupture-bound-slip-threshold


def load_slip_profile(q_dir: str) -> np.ndarray | None:
    rfile = os.path.join(q_dir, "fault.r.nc")
    if not os.path.exists(rfile):
        return None
    with xr.open_dataset(rfile) as ds:
        slips = ds["slips"].values.astype(float)
        slipd = ds["slipd"].values.astype(float)
    nrows = slips.shape[0]
    row = int(round((DEPTH_PCT / 100.0) * (nrows - 1)))
    row = max(0, min(row, nrows - 1))
    profile = np.hypot(slips, slipd)[row, :]
    # Zero out sub-threshold columns so small-slip noise doesn't accumulate
    profile[profile < SLIP_THRESHOLD_M] = 0.0
    return profile


def main() -> None:
    args = parse_args()
    case_dir = os.path.abspath(args.case_dir)
    case_name = os.path.basename(case_dir)
    outdir = os.path.abspath(args.outdir) if args.outdir else case_dir

    qdirs = find_q_dirs(case_dir)
    if not qdirs:
        raise RuntimeError(f"No Q* directories found in {case_dir}")

    end_times_yr = cumulative_end_times_yr(qdirs)
    profiles: list[np.ndarray] = []
    times_yr: list[float] = []

    for qd, t_yr in zip(qdirs, end_times_yr):
        if t_yr is None:
            print(f"Skipping {os.path.basename(qd)}: no global.dat")
            continue
        profile = load_slip_profile(qd)
        if profile is None:
            print(f"Skipping {os.path.basename(qd)}: no fault.r.nc")
            continue
        profiles.append(profile)
        times_yr.append(t_yr)
        print(f"{os.path.basename(qd)}: t={t_yr:.1f} yr  peak_slip={profile.max():.2f} m")

    if not profiles:
        raise RuntimeError("No slip profiles loaded")

    n_cols = profiles[0].size
    distance_km = np.arange(n_cols) * GRID_SPACING_M / 1e3

    cumulative = np.zeros(n_cols)
    running: list[np.ndarray] = []
    for p in profiles:
        cumulative = cumulative + p
        running.append(cumulative.copy())

    # Trim to active rupture zone (0.5 m threshold consistent with process.eqquasi.py)
    active = np.where(running[-1] > 0.5)[0]
    col_lo, col_hi = (int(active[0]), int(active[-1])) if active.size else (0, n_cols - 1)
    distance_km = distance_km[col_lo : col_hi + 1]
    profiles = [p[col_lo : col_hi + 1] for p in profiles]
    running = [r[col_lo : col_hi + 1] for r in running]

    t_arr = np.array(times_yr)
    cmap = plt.cm.viridis
    t_norm = matplotlib.colors.Normalize(vmin=t_arr.min(), vmax=t_arr.max())
    depth_row = int(round((DEPTH_PCT / 100.0) * 25))
    subtitle = f"{len(profiles)} earthquakes  |  {DEPTH_PCT:.0f}% depth  |  {case_name}"

    # --- Figure 1: individual earthquake slip profiles ---
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for profile, t_yr in zip(profiles, times_yr):
        ax.plot(distance_km, profile, color=cmap(t_norm(t_yr)), linewidth=LINE_WIDTH * 0.6, alpha=0.8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=t_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Elapsed Time (yr)", fontweight="bold")
    cbar.ax.tick_params(labelsize=11)
    ax.set_xlabel("Along-Strike Distance (km)")
    ax.set_ylabel("Coseismic Slip (m)")
    ax.set_title("Individual Earthquake Slip Profiles")
    ax.grid(True, alpha=0.3)
    style_axis(ax)
    fig.suptitle(subtitle, fontsize=11, fontweight="bold")
    fig.tight_layout()
    out1 = os.path.join(outdir, "eqquasi_slip_profiles.png")
    fig.savefig(out1)
    plt.close(fig)
    print(f"Saved: {out1}")

    # --- Figure 2: running accumulated slip ---
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for run_profile, t_yr in zip(running, times_yr):
        ax.plot(distance_km, run_profile, color=cmap(t_norm(t_yr)), linewidth=LINE_WIDTH * 0.6, alpha=0.8)
    ax.plot(distance_km, running[-1], color="black", linewidth=LINE_WIDTH,
            label=f"Total ({t_arr.max():.0f} yr)")
    sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=t_norm)
    sm2.set_array([])
    cbar2 = fig.colorbar(sm2, ax=ax, pad=0.02)
    cbar2.set_label("Elapsed Time (yr)", fontweight="bold")
    cbar2.ax.tick_params(labelsize=11)
    ax.set_xlabel("Along-Strike Distance (km)")
    ax.set_ylabel("Accumulated Slip (m)")
    ax.set_title("Running Accumulated Coseismic Slip")
    ax.legend(loc="upper left", frameon=False)
    ax.grid(True, alpha=0.3)
    style_axis(ax)
    fig.suptitle(subtitle, fontsize=11, fontweight="bold")
    fig.tight_layout()
    out2 = os.path.join(outdir, "eqquasi_accumulated_slip.png")
    fig.savefig(out2)
    plt.close(fig)
    print(f"Saved: {out2}")

    # --- Figure 3: long-term slip rate ---
    slip_rate_mm_yr = running[-1] / t_arr.max() * 1e3
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.plot(distance_km, slip_rate_mm_yr, color="black", linewidth=LINE_WIDTH)
    ax.set_xlabel("Along-Strike Distance (km)")
    ax.set_ylabel("Slip Rate (mm/yr)")
    ax.set_title(f"Long-Term Coseismic Slip Rate ({t_arr.max():.0f} yr)")
    ax.grid(True, alpha=0.3)
    style_axis(ax)
    fig.suptitle(subtitle, fontsize=11, fontweight="bold")
    fig.tight_layout()
    out3 = os.path.join(outdir, "eqquasi_slip_rate.png")
    fig.savefig(out3)
    plt.close(fig)
    print(f"Saved: {out3}")


if __name__ == "__main__":
    main()
