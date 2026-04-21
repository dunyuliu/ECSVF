#!/usr/bin/env python3
"""
Script to plot rupture length evolution using plotAccumulated functionality:
- Gets accumulated slip data for horizontal profile at 100% depth (surface)
- Calculates rupture length and nucleation location
- Plots vertical lines for rupture length vs accumulated time from global.dat
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os, glob, re
from lib import boxcar_average

# Set matplotlib defaults for better readability
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'font.weight': 'bold',
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold'
})

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    return [tryint(c) for c in re.split('([0-9]+)', s)]

def sort_nicely(l):
    l.sort(key=alphanum_key, reverse=False)
    return l

def read_global_dat_times(model_path):
    """
    Read time data from global.dat file
    Returns array of time values (first column)
    """
    global_file = os.path.join(model_path, 'global.dat')
    if not os.path.exists(global_file):
        print(f"Warning: {global_file} not found")
        return np.array([])
    
    try:
        # Read the first column (time) from global.dat
        data = np.loadtxt(global_file, usecols=[0])
        return data
    except Exception as e:
        print(f"Error reading {global_file}: {e}")
        return np.array([])

def generate_timeseries_profile_all(model_path, mode, var, ninterval, nwindow, depth_pct=None, horizontal_pct=None):
    """
    Modified from plotAccumulated - extract time series information from model folders
    """
    filenames = glob.glob(model_path+'/fault.?????.nc')
    filenames = sort_nicely(filenames)
    num_files = len(filenames)

    # fault.r.nc must exist - extract dimension information from restart file
    a = xr.open_dataset(model_path+'/fault.r.nc')
    shear_strike = a.shear_strike
    num_rows = shear_strike.shape[0]
    num_cols = shear_strike.shape[1]
    
    # For horizontal profile at surface (100% depth)
    if mode == 2:  # Horizontal profile
        if depth_pct is not None:
            row_index = int(round((depth_pct / 100.0) * (num_rows - 1)))
            row_index = max(0, min(row_index, num_rows - 1))
            actual_depth_pct = (row_index / (num_rows - 1)) * 100
        else:
            row_index = num_rows - 1  # Surface (100% depth)
            actual_depth_pct = 100.0
        col_index = None
        actual_horizontal_pct = None
        print(f'Horizontal profile at {actual_depth_pct:.1f}% depth (row {row_index})')

    print(f'Grid dimensions: {num_files} files, {num_rows} rows, {num_cols} cols')

    # Plot horizontal profile, count columns
    data = np.zeros((num_cols, num_files))
    data_ave = np.zeros((num_cols-nwindow+1, num_files))
    
    itag = 0
    for file1 in filenames:
        b = xr.open_dataset(file1)
        if var == 'slips':
            res = b.slips
        elif var == 'slipd':
            res = b.slipd
            
        itag = itag + 1
        for i in range(num_cols):
            data[i,itag-1] = res[row_index,i]
        
        data_ave[:,itag-1] = boxcar_average(data[:,itag-1],nwindow) 
        
    data0 = data_ave[:,::ninterval]
    return data0, actual_depth_pct, actual_horizontal_pct

def detect_rupture_and_calculate_extent(model_path, slip_data, slip_rate_threshold=0.1, slip_threshold=0.1, grid_spacing=2000.0):
    """
    Detect rupture using slip rate threshold and calculate extent using slip threshold
    
    Parameters:
    - model_path: path to Q cycle folder
    - slip_data: accumulated slip data [num_strike_nodes, num_time_steps] 
    - slip_rate_threshold: minimum slip rate to consider as active rupture (m/s)
    - slip_threshold: minimum slip to determine rupture extent (m)
    - grid_spacing: spatial resolution in meters
    
    Returns:
    - rupture_start: start position of rupture (km), 0 if no rupture detected
    - rupture_end: end position of rupture (km), 0 if no rupture detected  
    - nucleation_location: where slip rate first exceeds threshold (km), 0 if no rupture
    """
    # Get fault files to check slip rates and find nucleation
    filenames = glob.glob(model_path+'/fault.?????.nc')
    filenames = sort_nicely(filenames)
    
    if len(filenames) == 0:
        return 0.0, 0.0, 0.0
    
    # Find nucleation location (where slip rate first exceeds threshold)
    nucleation_location = 0.0
    rupture_detected = False
    
    for filename in filenames:
        with xr.open_dataset(filename) as ds:
            if 'slip_rate' in ds:
                slip_rate = ds.slip_rate.values
                # Check surface row (same as slip data)
                surface_slip_rate = slip_rate[-1, :]  # Last row = surface
                high_rate_indices = np.where(surface_slip_rate > slip_rate_threshold)[0]
                
                if len(high_rate_indices) > 0:
                    # First occurrence of high slip rate is nucleation
                    nucleation_idx = high_rate_indices[0]
                    nucleation_location = nucleation_idx * grid_spacing / 1000.0  # convert to km
                    rupture_detected = True
                    break
    
    if not rupture_detected:
        print(f"  No rupture detected (max slip rate < {slip_rate_threshold} m/s)")
        return 0.0, 0.0, 0.0
    
    # If rupture detected, calculate extent using final slip distribution
    final_slip_profile = slip_data[:, -1]  # Last time step
    
    # Find rupture extent using slip threshold
    ruptured_mask = final_slip_profile > slip_threshold
    if np.any(ruptured_mask):
        ruptured_indices = np.where(ruptured_mask)[0]
        
        # Rupture start and end positions
        rupture_start = ruptured_indices[0] * grid_spacing / 1000.0  # convert to km
        rupture_end = ruptured_indices[-1] * grid_spacing / 1000.0  # convert to km
        rupture_length = rupture_end - rupture_start
        
        print(f"  Rupture detected: {rupture_start:.1f} to {rupture_end:.1f} km (length = {rupture_length:.1f} km), nucleation at {nucleation_location:.1f} km")
    else:
        print(f"  Rupture detected but no slip > {slip_threshold} m for extent calculation")
        rupture_start = 0.0
        rupture_end = 0.0
        nucleation_location = 0.0
    
    return rupture_start, rupture_end, nucleation_location

def plot_rupture_length_evolution(ncyc, nwindow, var_name='slips', 
                                slip_rate_threshold=0.1, slip_threshold=0.1, grid_spacing=2000.0,
                                output_file='rupture_length_evolution.png'):
    """
    Plot ONE vertical line per Q cycle showing final rupture length
    """
    
    print(f'Generating rupture length evolution from accumulated {var_name}')
    
    ninterval = 1  # Use every time step
    depth_pct = 100.0  # Surface
    
    cycle_times = []
    cycle_rupture_starts = []
    cycle_rupture_ends = []
    cycle_nucleation_locations = []
    cumulative_time_offset = 0.0
    
    for i in range(ncyc):
        model_path = 'Q'+str(i) 
        print(f'Processing cycle {i}: {model_path}')
        
        try:
            # Read time data from global.dat for this cycle
            global_times = read_global_dat_times(model_path)
            if len(global_times) == 0:
                print(f"No time data found for {model_path}")
                continue
            
            # Get accumulated slip data for this cycle
            data, actual_depth_pct, actual_horizontal_pct = generate_timeseries_profile_all(
                model_path, mode=2, var=var_name, ninterval=ninterval, nwindow=nwindow, 
                depth_pct=depth_pct, horizontal_pct=None)
            
            # Detect rupture and calculate extent using proper thresholds
            rupture_start, rupture_end, nucleation_location = detect_rupture_and_calculate_extent(
                model_path, data, slip_rate_threshold, slip_threshold, grid_spacing)
            
            # Use the final time from global.dat for this cycle, accumulated across cycles
            final_time = global_times[-1] + cumulative_time_offset
            
            # Store one point per cycle
            cycle_times.append(final_time)
            cycle_rupture_starts.append(rupture_start)
            cycle_rupture_ends.append(rupture_end)
            cycle_nucleation_locations.append(nucleation_location)
            
            # Update cumulative time offset for next cycle
            cumulative_time_offset = final_time
            
            rupture_length = rupture_end - rupture_start if rupture_end > rupture_start else 0.0
            print(f"Cycle {i}: final time = {final_time:.2e}, rupture length = {rupture_length:.2f} km")
            
        except Exception as e:
            print(f'Error processing {model_path}: {e}')
            continue
    
    if not cycle_times:
        print("No data to plot")
        return
    
    # Create the plot - ONE vertical line per Q cycle
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
    
    # Convert times from seconds to years (1 year = 365.25 * 24 * 3600 seconds)
    seconds_per_year = 365.25 * 24 * 3600
    cycle_times_years = [t / seconds_per_year for t in cycle_times]
    
    # Plot: Vertical lines showing rupture extent (from start to end position)
    rupture_count = 0
    for i in range(len(cycle_times_years)):
        time_years = cycle_times_years[i]
        rupture_start = cycle_rupture_starts[i]
        rupture_end = cycle_rupture_ends[i]
        nucleation_location = cycle_nucleation_locations[i]
        
        # Only plot if rupture detected
        if rupture_end > rupture_start:
            # Rupture extent line
            ax1.plot([time_years, time_years], [rupture_start, rupture_end], 'b-', linewidth=3.0, alpha=0.8, label=f'Q{i}')
            # Nucleation point
            if nucleation_location > 0:
                ax1.plot(time_years, nucleation_location, 'ro', markersize=6, alpha=0.9)
            rupture_count += 1
    
    # Labels and formatting
    slip_type = 'Strike Slip' if var_name == 'slips' else 'Dip Slip'
    
    ax1.set_ylabel('Position along fault (km)', fontweight='bold')
    ax1.set_xlabel('Accumulated Time (years)', fontweight='bold')
    ax1.set_title(f'Rupture Extent per Q Cycle - {slip_type} at Surface', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)
    
    # Add legend for nucleation points if any ruptures detected
    if rupture_count > 0:
        ax1.plot([], [], 'ro', markersize=6, label='Nucleation')
        ax1.legend()
    
    # Format x-axis for scientific notation if needed
    ax1.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f'Saved: {output_file}')
    plt.show()

def main():
    """Main function with user interaction"""
    print("=== Rupture Length Evolution Plotter ===")
    print("This script calculates and plots rupture length evolution:")
    print("- Uses slip rate > threshold to detect rupture occurrence")
    print("- Uses slip > threshold to determine rupture extent")
    print("- Plots vertical lines vs accumulated time from global.dat")
    print()
    
    # Get number of cycles
    while True:
        try:
            ncyc = int(input("Number of cycles to process: "))
            if ncyc > 0:
                break
            print("Number of cycles must be positive")
        except ValueError:
            print("Please enter a positive integer")
    
    # Check if cycle directories exist
    existing_cycles = []
    for i in range(ncyc):
        model_dir = f'Q{i}'
        if os.path.exists(model_dir):
            existing_cycles.append(i)
            print(f"Found: {model_dir}")
        else:
            print(f"Warning: Directory not found: {model_dir}")
    
    if not existing_cycles:
        print("No valid model directories found. Exiting.")
        return
    
    # Get slip variable
    while True:
        var_name = input("Slip variable (slips/slipd) [default: slips]: ").strip().lower()
        if var_name == "":
            var_name = "slips"
        if var_name in ['slips', 'slipd']:
            break
        print("Please enter 'slips' or 'slipd'")
    
    # Get window size for averaging
    while True:
        try:
            nwindow = int(input("Window size for averaging [default: 1]: ") or "1")
            if nwindow > 0:
                break
            print("Window size must be positive")
        except ValueError:
            print("Please enter a positive integer")
    
    # Get parameters
    while True:
        try:
            slip_rate_threshold = float(input("Slip rate threshold for rupture detection (m/s) [default: 0.1]: ") or "0.1")
            if slip_rate_threshold > 0:
                break
            print("Slip rate threshold must be positive")
        except ValueError:
            print("Please enter a valid number")
    
    while True:
        try:
            slip_threshold = float(input("Slip threshold for rupture extent (m) [default: 0.1]: ") or "0.1")
            if slip_threshold > 0:
                break
            print("Slip threshold must be positive")
        except ValueError:
            print("Please enter a valid number")
    
    while True:
        try:
            grid_spacing = float(input("Grid spacing (m) [default: 2000]: ") or "2000")
            if grid_spacing > 0:
                break
            print("Grid spacing must be positive")
        except ValueError:
            print("Please enter a valid number")
    
    # Get output filename
    output_file = input("Output filename [default: rupture_length_evolution.png]: ").strip()
    if not output_file:
        output_file = "rupture_length_evolution.png"
    
    # Create the plot
    plot_rupture_length_evolution(ncyc, nwindow, var_name, slip_rate_threshold, slip_threshold, grid_spacing, output_file)

if __name__ == "__main__":
    main()