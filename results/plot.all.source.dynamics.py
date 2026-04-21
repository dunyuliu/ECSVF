#! /usr/bin/env python3

import os
import glob

def clean_files(folder):
    """Remove all .png and .gif files in the given folder."""
    for file_pattern in ["*.png", "*.gif"]:
        for file_path in glob.glob(os.path.join(folder, file_pattern)):
            os.remove(file_path)
            print(f"Removed: {file_path}")

def process_folders():
    """Loop over folders Q0, Q1, Q2, ..., QN and run ../plotOnFaultVars_mod."""
    base_dir = os.getcwd()
    q_folders = sorted([d for d in os.listdir(base_dir) if d.startswith("Q") and d[1:].isdigit()],
                       key=lambda x: int(x[1:]))

    for folder in q_folders:
        folder_path = os.path.join(base_dir, folder)
        if os.path.isdir(folder_path):
            print(f"Processing folder: {folder}")
            os.chdir(folder_path)
            
            # Clean .png and .gif files
            clean_files(folder_path)
            
            # Run the command
            os.system("../plotOnFaultVars --cores 8")
            
            # Return to the base directory
            os.chdir(base_dir)

if __name__ == "__main__":
    process_folders()
