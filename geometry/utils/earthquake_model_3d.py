import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_geometry_file(geometry_filename):
    """
    Load geometry data from processed .txt file
    Format:
    Line 1: nx nz 0
    Line 2: resolution min_x min_z
    Line 3+: y(m) dy/dx dy/dz
    """
    print(f"Loading geometry file: {geometry_filename}")
    
    with open(geometry_filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header lines
    header1 = lines[0].split()
    nx, nz = int(header1[0]), int(header1[1])
    
    header2 = lines[1].split()
    resolution, x_min, z_min = float(header2[0]), float(header2[1]), float(header2[2])
    
    print(f"Grid parameters: {nx} × {nz}, resolution: {resolution} km")
    print(f"Grid origin: x_min={x_min} km, z_min={z_min} km")
    
    # Load data starting from line 3 (index 2)
    data_lines = [line.split() for line in lines[2:]]
    data = np.array([[float(val) for val in line] for line in data_lines])
    
    # Data format: y(m) dy/dx dy/dz
    y_values = data[:, 0] / 1000.0  # Convert from meters to km
    derivatives = data[:, 1:]       # dy/dx, dy/dz (dimensionless)
    
    # Reconstruct x,z coordinates using header information
    x_max = x_min + (nx - 1) * resolution
    z_max = z_min + (nz - 1) * resolution
    
    x_grid = np.arange(x_min, x_max + resolution, resolution)
    z_grid = np.arange(z_max, z_min - resolution, -resolution)  # Top to bottom
    
    print(f"Reconstructed grid ranges:")
    print(f"X: {x_min} to {x_max} km ({len(x_grid)} points)")
    print(f"Z: {z_max} to {z_min} km ({len(z_grid)} points)")
    
    # The data is ordered: bottom z first, then -x to +x
    vertices = []
    for i, y in enumerate(y_values):
        # Calculate z_idx and x_idx from linear index
        z_idx = i // nx
        x_idx = i % nx
        
        # Convert to actual coordinates
        # Note: z_grid goes from top to bottom, but data is ordered bottom first
        # So we reverse the z index
        actual_z_idx = nz - 1 - z_idx
        
        if actual_z_idx < len(z_grid) and x_idx < len(x_grid):
            x = x_grid[x_idx]
            z = z_grid[actual_z_idx]
            vertices.append([x, z, y])
    
    vertices = np.array(vertices)
    
    print(f"Loaded {len(vertices)} vertices")
    print(f"X range: {np.min(vertices[:, 0]):.1f} to {np.max(vertices[:, 0]):.1f} km")
    print(f"Y range: {np.min(vertices[:, 2]):.1f} to {np.max(vertices[:, 2]):.1f} km")
    print(f"Z range: {np.min(vertices[:, 1]):.1f} to {np.max(vertices[:, 1]):.1f} km")
    
    return vertices, derivatives

def create_earthquake_model_case1():
    """
    Case 1: Simple planar fault
    """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Model dimensions
    x_min, x_max = -320, 320
    y_pos, y_neg = 30, -60
    z_depth = 50
    
    # Draw box edges
    x = [x_min, x_max, x_max, x_min, x_min]
    y_top = [y_pos, y_pos, y_pos, y_pos, y_pos]
    y_bottom = [y_neg, y_neg, y_neg, y_neg, y_neg]
    z_top = [0, 0, 0, 0, 0]
    z_bottom = [-z_depth, -z_depth, -z_depth, -z_depth, -z_depth]
    
    ax.plot(x, y_top, z_top, 'k-', linewidth=1, alpha=0.6)
    ax.plot(x, y_bottom, z_top, 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_min, x_min], [y_neg, y_pos], [0, 0], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_neg, y_pos], [0, 0], 'k-', linewidth=1, alpha=0.6)
    
    ax.plot(x, y_top, z_bottom, 'k-', linewidth=1, alpha=0.6)
    ax.plot(x, y_bottom, z_bottom, 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_min, x_min], [y_neg, y_pos], [-z_depth, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_neg, y_pos], [-z_depth, -z_depth], 'k-', linewidth=1, alpha=0.6)
    
    ax.plot([x_min, x_min], [y_pos, y_pos], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_pos, y_pos], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_min, x_min], [y_neg, y_neg], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_neg, y_neg], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    
    # Create planar fault
    dip_angle = 50
    dip_rad = np.radians(dip_angle)
    
    x_fault = np.linspace(x_min, x_max, 200)
    z_fault = np.linspace(0, -z_depth, 100)
    X_fault, Z_fault = np.meshgrid(x_fault, z_fault)
    
    horizontal_offset = Z_fault / np.tan(dip_rad)
    Y_fault = horizontal_offset
    
    # Define zones with updated VS definition
    vw_mask = (Z_fault >= -20) & (Z_fault <= -2) & (np.abs(X_fault) <= 300)
    
    # VS zone: 2 km around VW in depth AND 2 km on ±x sides of VW
    vs_mask_depth = (Z_fault >= -22) & (Z_fault <= 0) & (np.abs(X_fault) <= 300)  # 2 km above/below VW
    vs_mask_sides = (Z_fault >= -20) & (Z_fault <= -2) & (np.abs(X_fault) > 300) & (np.abs(X_fault) <= 302)  # 2 km on sides
    vs_mask = (vs_mask_depth | vs_mask_sides) & ~vw_mask
    
    creep_mask = ~(vw_mask | vs_mask)
    
    # Plot zones
    if np.any(vw_mask):
        X_vw = np.where(vw_mask, X_fault, np.nan)
        Y_vw = np.where(vw_mask, Y_fault, np.nan)
        Z_vw = np.where(vw_mask, Z_fault, np.nan)
        ax.plot_surface(X_vw, Y_vw, Z_vw, alpha=0.8, color='yellow')
    
    if np.any(vs_mask):
        X_vs = np.where(vs_mask, X_fault, np.nan)
        Y_vs = np.where(vs_mask, Y_fault, np.nan)
        Z_vs = np.where(vs_mask, Z_fault, np.nan)
        ax.plot_surface(X_vs, Y_vs, Z_vs, alpha=0.8, color='green')
    
    if np.any(creep_mask):
        X_creep = np.where(creep_mask, X_fault, np.nan)
        Y_creep = np.where(creep_mask, Y_fault, np.nan)
        Z_creep = np.where(creep_mask, Z_fault, np.nan)
        ax.plot_surface(X_creep, Y_creep, Z_creep, alpha=0.6, color='gray')
    
    # Fault trace
    x_trace = np.linspace(x_min, x_max, 100)
    ax.plot(x_trace, np.zeros_like(x_trace), np.zeros_like(x_trace), 'r-', linewidth=3)
    
    # Annotations
    ax.text(0, y_pos + 10, 10, 'Strike Direction (X)', fontsize=10)
    ax.text(x_min - 50, 0, 0, 'Fault Normal\nDirection (Y)', fontsize=10)
    ax.text(0, 0, 20, 'Depth (Z)', fontsize=10)
    ax.text(0, -20, -25, f'Fault Dip: {dip_angle}°', fontsize=10, color='red')
    
    ax.text(-200, -15, -10, 'VW', fontsize=12, color='orange', weight='bold')
    ax.text(-200, -15, -1, 'VS', fontsize=12, color='darkgreen', weight='bold')
    ax.text(-200, -5, -30, 'Creeping', fontsize=12, color='darkgray', weight='bold')
    
    ax.set_xlabel('Strike Direction (km)')
    ax.set_ylabel('Fault Normal Direction (km)')
    ax.set_zlabel('Depth (km)')
    
    # Static dimensions for Case 1
    x_width = x_max - x_min
    y_width = y_pos - y_neg
    ax.set_title(f'planar\n{x_width:.0f}km × {y_width:.0f}km × {z_depth}km')
    
    ax.set_xlim(x_min - 50, x_max + 50)
    ax.set_ylim(y_neg - 20, y_pos + 20)
    ax.set_zlim(-z_depth - 10, 10)
    
    ax.view_init(elev=20, azim=-45)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, ax

def create_earthquake_model_from_geometry(geometry_filename):
    """
    Create earthquake model from processed geometry file
    """
    # Load geometry data
    vertices, derivatives = load_geometry_file(geometry_filename)
    
    # Extract coordinates (vertices format: [x, z, y])
    x_coords = vertices[:, 0]  # Strike direction
    z_coords = vertices[:, 1]  # Depth direction  
    y_coords = vertices[:, 2]  # Fault normal direction
    
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Dynamic model box dimensions based on actual fault extent
    fault_x_min, fault_x_max = np.min(x_coords), np.max(x_coords)
    fault_y_min, fault_y_max = np.min(y_coords), np.max(y_coords)
    fault_z_min, fault_z_max = np.min(z_coords), np.max(z_coords)
    
    # Add buffer around fault (50 km on each side for x, 30 km for y)
    x_buffer = 50
    y_buffer = 30
    
    x_min = fault_x_min - x_buffer
    x_max = fault_x_max + x_buffer
    y_pos = max(30, fault_y_max + y_buffer)
    y_neg = fault_y_min - y_buffer
    z_depth = abs(fault_z_min)
    
    print(f"Model box dimensions:")
    print(f"X: {x_min:.1f} to {x_max:.1f} km (fault: {fault_x_min:.1f} to {fault_x_max:.1f} km)")
    print(f"Y: {y_neg:.1f} to {y_pos:.1f} km (fault: {fault_y_min:.1f} to {fault_y_max:.1f} km)")
    print(f"Z: -{z_depth:.1f} to 0 km")
    
    # Draw box edges
    x = [x_min, x_max, x_max, x_min, x_min]
    y_top = [y_pos, y_pos, y_pos, y_pos, y_pos]
    y_bottom = [y_neg, y_neg, y_neg, y_neg, y_neg]
    z_top = [0, 0, 0, 0, 0]
    z_bottom = [-z_depth, -z_depth, -z_depth, -z_depth, -z_depth]
    
    ax.plot(x, y_top, z_top, 'k-', linewidth=1, alpha=0.6)
    ax.plot(x, y_bottom, z_top, 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_min, x_min], [y_neg, y_pos], [0, 0], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_neg, y_pos], [0, 0], 'k-', linewidth=1, alpha=0.6)
    
    ax.plot(x, y_top, z_bottom, 'k-', linewidth=1, alpha=0.6)
    ax.plot(x, y_bottom, z_bottom, 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_min, x_min], [y_neg, y_pos], [-z_depth, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_neg, y_pos], [-z_depth, -z_depth], 'k-', linewidth=1, alpha=0.6)
    
    ax.plot([x_min, x_min], [y_pos, y_pos], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_pos, y_pos], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_min, x_min], [y_neg, y_neg], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    ax.plot([x_max, x_max], [y_neg, y_neg], [0, -z_depth], 'k-', linewidth=1, alpha=0.6)
    
    # Color vertices by zone with updated VS definition
    colors = []
    for i in range(len(vertices)):
        x, z, y = vertices[i]
        
        # VW zone: 2-20 km depth, |x| <= 300 km
        if -20 <= z <= -2 and abs(x) <= 300:
            colors.append('yellow')  # VW zone
        
        # VS zone: 2 km around VW in depth AND 2 km on ±x sides of VW
        elif ((-22 <= z <= 0 and abs(x) <= 300) or  # 2 km above/below VW in depth
              (-20 <= z <= -2 and 300 < abs(x) <= 302)):  # 2 km on ±x sides of VW
            colors.append('green')   # VS zone
        
        # Everything else: creeping
        else:
            colors.append('gray')    # Creeping zone
    
    # Plot fault vertices
    ax.scatter(x_coords, y_coords, z_coords, c=colors, s=0.5, alpha=0.7)
    
    # Add fault trace at surface
    surface_verts = vertices[z_coords > -2]  # Within 2km of surface
    if len(surface_verts) > 1:
        sorted_indices = np.argsort(surface_verts[:, 0])
        sorted_surface = surface_verts[sorted_indices]
        ax.plot(sorted_surface[:, 0], sorted_surface[:, 2], sorted_surface[:, 1], 
                'r-', linewidth=3, alpha=0.8, label='Fault Trace')
    
    # Add annotations
    ax.text(0, y_pos + 10, 10, 'Strike Direction (X)', fontsize=10)
    ax.text(x_min - 50, 0, 0, 'Fault Normal\nDirection (Y)', fontsize=10)
    ax.text(0, 0, 20, 'Depth (Z)', fontsize=10)
    
    # Add zone labels
    ax.text(-200, -15, -10, 'VW', fontsize=12, color='orange', weight='bold')
    ax.text(-200, -15, -1, 'VS', fontsize=12, color='darkgreen', weight='bold')
    ax.text(-200, -5, -30, 'Creeping', fontsize=12, color='darkgray', weight='bold')
    
    # Set labels and title
    ax.set_xlabel('Strike Direction (km)')
    ax.set_ylabel('Fault Normal Direction (km)')
    ax.set_zlabel('Depth (km)')
    
    # Dynamic title with geometry filename and actual dimensions
    x_width = x_max - x_min
    y_width = y_pos - y_neg
    geometry_name = geometry_filename.replace('_Geometry.txt', '').replace('.txt', '')
    ax.set_title(f'{geometry_name}\n{x_width:.0f}km × {y_width:.0f}km × {z_depth:.0f}km')
    
    # Set axis limits
    ax.set_xlim(x_min - 50, x_max + 50)
    ax.set_ylim(y_neg - 20, y_pos + 20)
    ax.set_zlim(-z_depth - 10, 10)
    
    # Set view angle: -x left, +x right, y- near viewer
    ax.view_init(elev=20, azim=-45)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, ax

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        case = sys.argv[1].lower()
        
        if case == "case1":
            # Case 1: Planar fault
            print("Creating Case 1: Planar fault geometry...")
            fig, ax = create_earthquake_model_case1()
            
            plt.savefig('earthquake_model_case1.png', dpi=300, bbox_inches='tight')
            print("Case 1 diagram saved as 'earthquake_model_case1.png'")
            plt.close(fig)
            
        elif case == "case2":
            # Case 2: no_dip_change_1km geometry
            geometry_file = "no_dip_change_1km_Geometry.txt"
            print("Creating Case 2: no_dip_change_1km geometry...")
            
            try:
                fig, ax = create_earthquake_model_from_geometry(geometry_file)
                
                plt.savefig('earthquake_model_case2.png', dpi=300, bbox_inches='tight')
                print("Case 2 diagram saved as 'earthquake_model_case2.png'")
                plt.close(fig)
                
            except Exception as e:
                print(f"Error processing Case 2: {e}")
                print(f"Make sure {geometry_file} exists in the current directory")
                
        elif case == "case3":
            # Case 3: variable_dip_simpler_geom1km_hws geometry
            geometry_file = "variable_dip_simpler_geom1km_hws_Geometry.txt"
            print("Creating Case 3: variable_dip_simpler_geom1km_hws geometry...")
            
            try:
                fig, ax = create_earthquake_model_from_geometry(geometry_file)
                
                plt.savefig('earthquake_model_case3.png', dpi=300, bbox_inches='tight')
                print("Case 3 diagram saved as 'earthquake_model_case3.png'")
                plt.close(fig)
                
            except Exception as e:
                print(f"Error processing Case 3: {e}")
                print(f"Make sure {geometry_file} exists in the current directory")
        
        else:
            print("Usage: python earthquake_model_3d.py [case1|case2|case3]")
            print("  case1: Planar fault")
            print("  case2: no_dip_change_1km geometry")
            print("  case3: variable_dip_simpler_geom1km_hws geometry")
            
    else:
        print("Usage: python earthquake_model_3d.py [case1|case2|case3]")
        print("  case1: Planar fault")
        print("  case2: no_dip_change_1km geometry")
        print("  case3: variable_dip_simpler_geom1km_hws geometry")