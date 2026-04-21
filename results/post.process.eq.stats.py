#! /usr/bin/env python3

import os
import numpy as np
import xarray as xr
from glob import glob
import re

def natural_sort_key(s):
    """Sort strings containing numbers in natural order."""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def process_fault_files(base_path, variable_name, slip_rate_threshold=0.01):
    # Find all Q folders
    q_folders = sorted(glob(os.path.join(base_path, 'Q*')), key=natural_sort_key)
    print(f"Found Q folders: {q_folders}")
    if not q_folders:
        print("No Q folders found.")
        return []

    all_values = []
    
    for q_folder in q_folders:
        # Read global.dat file to get times
        global_dat_path = os.path.join(q_folder, 'global.dat')
        try:
            times = np.loadtxt(global_dat_path)[:, 0]  # First column contains times
            peak_slip_rates = np.loadtxt(global_dat_path)[:, 1]  # Second column contains peak slip rates
            
            # Find start and end of events
            threshold = slip_rate_threshold  # You can adjust this threshold
            window_size = 20
            start_step = None
            end_step = None
            
            # Find start step
            for i in range(len(peak_slip_rates) - window_size):
                if all(peak_slip_rates[i:i+window_size] > threshold):
                    start_step = i
                    break
            
            # Find first end step after start_step
            if start_step is not None:
                for i in range(start_step, len(peak_slip_rates) - window_size):
                    if all(peak_slip_rates[i:i+window_size] < threshold):
                        end_step = i
                        break
            
            #if start_step is not None and end_step is not None:
            #    print(f"Event starts at step: {start_step}")
            #    print(f"Event ends at step: {end_step}")
            #    print(f"Earthquake lasts: {times[end_step] - times[start_step]} seconds")
        except Exception as e:
            print(f"Error reading {global_dat_path}: {str(e)}")
            continue
        
        # Find all fault files in each Q folder, excluding fault.r.nc
        fault_files = [f for f in sorted(glob(os.path.join(q_folder, 'fault.*.nc')), 
                           key=natural_sort_key) if not f.endswith('fault.r.nc')]
        
        for fault_file in fault_files:
            try:
                # Extract time step number from filename
                time_step = int(fault_file.split('.')[-2])  # Assumes format fault.XXXX.nc
                if time_step > len(times):
                    print(f"Warning: time step {time_step} exceeds available times in global.dat {len(times)}")
                    continue
                
                # Open the NetCDF file
                with xr.open_dataset(fault_file) as ds:
                    # Extract the specified variable
                    value = ds[variable_name].values
                    
                    # Store data along with metadata and time
                    all_values.append({
                        'folder': os.path.basename(q_folder),
                        'file': os.path.basename(fault_file),
                        'time': times[time_step-1],
                        'peak_slip_rate': peak_slip_rates[time_step-1],
                        'start_step': start_step,
                        'end_step': end_step,
                        variable_name: value
                    })
                    
            except Exception as e:
                print(f"Error processing {fault_file}: {str(e)}")
    
    return all_values

def calculate_shear_modulus(vs, rho):
    return (vs**2) * rho

def eq_analyzer(slip_rates, slip_strike, slip_dip, 
                slip_rate_threshold=1e-3,
                area_per_node=1.0, shear_modulus=3e10):
    """Analyze earthquake data to find slip rates, moments, and magnitudes."""
    # Group data by Q folder
    folder_data = {}
    for sr, ss, sd in zip(slip_rates, slip_strike, slip_dip):
        folder = sr['folder']
        if folder not in folder_data:
            folder_data[folder] = {'slip_rates': [], 'slip_strike': [], 'slip_dip': []}
        folder_data[folder]['slip_rates'].append(sr)
        folder_data[folder]['slip_strike'].append(ss)
        folder_data[folder]['slip_dip'].append(sd)
       
    results = {}
    for folder, data in folder_data.items():
        print(f"\nProcessing folder: {folder}")
        
        # Sort by time
        sorted_indices = sorted(range(len(data['slip_rates'])), 
                              key=lambda i: data['slip_rates'][i]['time'])
        times = [data['slip_rates'][i]['time'] for i in sorted_indices]
        rates = [data['slip_rates'][i]['slip_rate'] for i in sorted_indices]
        slip_strike_values = [data['slip_strike'][i]['slips'] for i in sorted_indices]
        slip_dip_values = [data['slip_dip'][i]['slipd'] for i in sorted_indices]
        
        # Find when slip rate exceeds threshold
        start_idx = next((i for i, rate in enumerate(rates) if np.any(rate > slip_rate_threshold)), None)
        if start_idx is None:
            print(f"  No slip rate exceeding threshold found in {folder}")
            continue
            
        # Find when slip rate drops below threshold after start
        end_idx = next((i for i in range(start_idx + 1, len(rates)) 
                       if np.all(rates[i] < slip_rate_threshold)), None)
        if end_idx is None:
            print(f"  No slip rate drop below threshold found in {folder}")
            continue
            
        # Find nodes that exceed slip rate threshold at any point
        slip_rate_mask = np.any(np.array(rates) > slip_rate_threshold, axis=0)  # Using slip_threshold instead of 0.01
        
        # Calculate total slip and count nodes where slip rate exceeded threshold
        total_slip = np.sqrt(slip_strike_values[end_idx]**2 + slip_dip_values[end_idx]**2)
        nodes_above_threshold = np.sum(slip_rate_mask)
        
        # Calculate average slip for nodes above threshold
        average_slip = np.mean(total_slip[slip_rate_mask]) if nodes_above_threshold > 0 else 0
        max_slip = np.max(total_slip[slip_rate_mask]) if nodes_above_threshold > 0 else 0
        
        # Calculate rupture area
        rupture_area = nodes_above_threshold * area_per_node
        
        # Calculate moment only for nodes that had slip rate above threshold
        moment = shear_modulus * rupture_area * average_slip

        # Calculate magnitude
        magnitude = 2/3 * np.log10(moment*1e7) - 10.7  # Using different constant
        
        print(f"  Start time: {times[start_idx]:.2f}")
        print(f"  End time: {times[end_idx]:.2f}")
        print(f"  Duration: {times[end_idx] - times[start_idx]:.2f}")
        print(f"  Nodes above threshold: {nodes_above_threshold}")
        print(f"  Average slip: {average_slip:.2f} m")
        print(f"  Moment: {moment:.2e}")
        print(f"  Magnitude: {magnitude:.2f}")
        print(f"  Rupture area: {rupture_area:.2f} m^2")
        
        results[folder] = {
            'start_time': times[start_idx],
            'end_time': times[end_idx],
            'duration': times[end_idx] - times[start_idx],
            'nodes_above_threshold': nodes_above_threshold,
            'average_slip': average_slip,
            'max_slip': max_slip,
            'moment': moment,
            'magnitude': magnitude,
            'rupture_area': rupture_area
        }
    
    return results

def write_stats_to_csv(statistics, output_file='earthquake_stats.csv'):
    """
    Write earthquake statistics to a CSV file.
    Args:
        statistics: Dictionary containing earthquake statistics
        output_file: Path to output CSV file
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write('event_id,t0_s,t0_year,m0_Nm,mw,area_m2,dt_s,max_slip_m,x_NZTM_m,y_NZTM_m,z_NZTM_m,sbound_NZTM_m,nbound_NZTM_m\n')
        
        # Process each event
        for folder, stats in statistics.items():
            event_id = int(re.search(r'Q(\d+)', folder).group(1))
            t0_s = stats['start_time']
            t0_year = stats['start_time']/3600/24/365  # Default year if not specified
            m0_Nm = stats['moment']
            mw = stats['magnitude']
            area_m2 = stats['rupture_area']
            dt_s = stats['duration']
            max_slip_m = stats['max_slip']  # Using average slip as max slip
            # Default values for spatial coordinates
            x_NZTM_m = y_NZTM_m = z_NZTM_m = sbound_NZTM_m = nbound_NZTM_m = 0
            
            f.write(f'{event_id},{t0_s},{t0_year},{m0_Nm},{mw},{area_m2},{dt_s},{max_slip_m},'
                    f'{x_NZTM_m},{y_NZTM_m},{z_NZTM_m},{sbound_NZTM_m},{nbound_NZTM_m}\n')
                
def main():
    # Set your base path
    base_path = "."
    slip_rate_threshold = 0.01  # Adjust as needed
    # Process all fault files

    slip_rates = process_fault_files(base_path, variable_name='slip_rate', slip_rate_threshold=slip_rate_threshold)
    slip_strike = process_fault_files(base_path, variable_name='slips', slip_rate_threshold=slip_rate_threshold)  # If slip is also needed
    slip_dip = process_fault_files(base_path, variable_name='slipd', slip_rate_threshold=slip_rate_threshold)  # If slip is also needed

    rho = 2.67e3
    vs = 3.464e3
    shear_modulus = calculate_shear_modulus(vs, rho)
    # Perform statistical analysis
    statistics = eq_analyzer(slip_rates, slip_strike, slip_dip,
                             slip_rate_threshold=slip_rate_threshold,
                             area_per_node=2e3*2e3, shear_modulus=shear_modulus)
    


    write_stats_to_csv(statistics)
    # Print results
    print("\nStatistical Analysis Results:")
    for key, value in statistics.items():
        print(f"{key}: {value}")

if __name__ == "__main__":
    main()