#!/usr/bin/env python3
"""
Convert raw HBI catalogue CSVs to benchmark bundle format.

Reads all job entries from Alpine_jobid_readme.txt, finds the matching
full-catalogue and station CSV files, remaps column names to the benchmark
standard, and writes *_event_info.csv / *_site{n}_event_info.csv bundles.

Output filenames encode the job ID and scenario, e.g.:
    alpine_hbi_81103019_planar_event_info.csv
    alpine_hbi_81102019_varying_dip_site1_event_info.csv

Usage:
    python benchmark/utils/process.hbi.py
    python benchmark/utils/process.hbi.py --input-dir <path> --outdir <path>
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = REPO_ROOT / "benchmark/Simulation_results/HBI_results"
DEFAULT_OUTPUT_DIR = DEFAULT_INPUT_DIR

JOBID_README = "Alpine_jobid_readme.txt"

GEOMETRY_TO_SCENARIO = {
    "planar": "planar",
    "non-planar": "varying_dip",
}

EVENT_COL_MAP = {
    "event_ID":       "event_id",
    "start_time_sec": "t0_s",
    "start_time_yr":  "t0_year",
    "magnitude_M0":   "m0_Nm",
    "magnitude_Mw":   "mw",
    "rupture_area_m2": "area_m2",
    "duration_sec":   "dt_s",
    "max_slip_m":     "max_slip_m",
    "hypo_x_NZTM":    "x_NZTM_m",
    "hypo_y_NZTM":    "y_NZTM_m",
    "hypo_z_NZTM":    "z_NZTM_m",
    "s_bound_NZTM":   "sbound_NZTM_m",
    "n_bound_NZTM":   "nbound_NZTM_m",
}

SITE_COL_MAP = {
    "event_ID":       "event_id",
    "event_time_sec": "t0_s",
    "magnitude_Mw":   "Mw",
    "event_slip_m":   "Slip_m",
    "rake":           "Rake_deg",
}

EVENT_FIELDNAMES = [
    "event_id", "t0_s", "t0_year", "m0_Nm", "mw", "area_m2",
    "dt_s", "max_slip_m", "x_NZTM_m", "y_NZTM_m", "z_NZTM_m",
    "sbound_NZTM_m", "nbound_NZTM_m",
]

SITE_FIELDNAMES = ["event_id", "t0_s", "Mw", "Slip_m", "Rake_deg"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert raw HBI catalogue CSVs to benchmark bundle format."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--outdir", "--output-dir", dest="output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def parse_jobid_readme(path: Path) -> list[dict[str, str]]:
    jobs: list[dict[str, str]] = []
    with path.open() as fh:
        header: list[str] | None = None
        for line in fh:
            cols = [c.strip() for c in re.split(r"\t+", line.rstrip("\n"))]
            cols = [c for c in cols if c]
            if not cols:
                continue
            if header is None:
                header = cols
                continue
            jobs.append({header[i]: cols[i] for i in range(min(len(header), len(cols)))})
    return jobs


def find_file(input_dir: Path, glob: str) -> Path:
    matches = sorted(input_dir.glob(glob))
    if not matches:
        raise FileNotFoundError(f"No file matching {glob!r} in {input_dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple files match {glob!r}: {[str(m) for m in matches]}")
    return matches[0]


def remap_rows(src: Path, col_map: dict[str, str], fieldnames: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with src.open() as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            out: dict[str, str] = {}
            for src_col, dst_col in col_map.items():
                if src_col in raw:
                    v = raw[src_col].strip()
                    out[dst_col] = v if v else "nan"
            for f in fieldnames:
                if f not in out:
                    out[f] = "nan"
            rows.append({f: out[f] for f in fieldnames})
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path.name}  ({len(rows)} rows)")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    readme = input_dir / JOBID_README
    if not readme.exists():
        raise FileNotFoundError(readme)

    jobs = parse_jobid_readme(readme)
    print(f"Found {len(jobs)} jobs in {JOBID_README}")

    for job in jobs:
        jobid = job.get("jobid", "").strip()
        geometry = job.get("geometry", "").strip()
        scenario = GEOMETRY_TO_SCENARIO.get(geometry)
        if scenario is None:
            print(f"  Skipping jobid {jobid!r}: unrecognised geometry {geometry!r}")
            continue

        prefix = f"alpine_hbi_{jobid}_{scenario}"
        print(f"\n  jobid={jobid}  geometry={geometry}  scenario={scenario}  prefix={prefix}")

        event_src = find_file(input_dir, f"Alpine_catalogue_summary_full_{jobid}_*.csv")
        event_rows = remap_rows(event_src, EVENT_COL_MAP, EVENT_FIELDNAMES)
        write_csv(output_dir / f"{prefix}_event_info.csv", event_rows, EVENT_FIELDNAMES)

        for n in (1, 2, 3):
            site_src = find_file(input_dir, f"Alpine_catalogue_summary_station{n}_{jobid}_*.csv")
            site_rows = remap_rows(site_src, SITE_COL_MAP, SITE_FIELDNAMES)
            write_csv(output_dir / f"{prefix}_site{n}_event_info.csv", site_rows, SITE_FIELDNAMES)


if __name__ == "__main__":
    main()
