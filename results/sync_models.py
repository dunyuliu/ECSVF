#!/usr/bin/env python3

import subprocess
import sys

model_names = [
    "nz.bp5.qdc.dip50.2000.norm_6mm_yr",
    "nz.bp5.qdc.dip50.2000.norm0",
    "nz.bp5.qdc.noDipChange.2000.norm0.slowInitialLoad",
    "nz.bp5.qdc.noDipChange.2000.norm_6mm_yr.slowInitialLoad",
    "nz.bp5.qdc.varyDip.2000.norm_6mm_yr.slowInitialLoad",
    "nz.bp5.qdc.varyDip20251202.2000.norm_6mm_yr.slowInitialLoad"
]

def show_usage():
    print("Usage: ./sync_models.py <index> [--rsync]")
    print("Available models:")
    for i, name in enumerate(model_names):
        print(f"  {i}: {name}")
    print("\nOptions:")
    print("  --rsync    Use rsync instead of scp (recommended for resumable transfers)")

if len(sys.argv) < 2 or len(sys.argv) > 3:
    show_usage()
    sys.exit(1)

use_rsync = "--rsync" in sys.argv
index_arg = sys.argv[1] if sys.argv[1] != "--rsync" else sys.argv[2]

try:
    index = int(index_arg)
    if index < 0 or index >= len(model_names):
        raise ValueError
except ValueError:
    print(f"Invalid index. Please use 0-{len(model_names)-1}")
    sys.exit(1)

model_name = model_names[index]
remote_path = f"dunyuliu@ls6.tacc.utexas.edu:/scratch/07931/dunyuliu/0.Dunyu/EQquasi/nz.work/{model_name}"

if use_rsync:
    sync_command = f"rsync -avz --progress {remote_path}/ {model_name}/"
    print(f"Syncing model with rsync: {model_name}")
else:
    sync_command = f"scp -r {remote_path} ."
    print(f"Syncing model with scp: {model_name}")

subprocess.run(sync_command, shell=True)
