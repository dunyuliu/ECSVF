import numpy as np
import re
from scipy.interpolate import griddata
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RBFInterpolator, interp2d
from scipy.optimize import minimize
import os

def read_ascii_stl_vertices(filename, max_vertices=50000):
    """
    Read ASCII STL file and return unique vertices
    """
    vertices_set = set()
    vertex_count = 0
    
    print(f"Reading ASCII STL file: {filename}")
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num % 100000 == 0 and line_num > 0:
                print(f"Processing line {line_num}, found {len(vertices_set)} unique vertices")
            
            # Look for vertex lines
            if 'vertex' in line:
                # Extract coordinates using regex
                match = re.search(r'vertex\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if match:
                    x, y, z = map(float, match.groups())
                    vertices_set.add((x, y, z))
                    vertex_count += 1
                    
                    # Limit vertices for performance
                    if len(vertices_set) >= max_vertices:
                        print(f"Reached maximum vertex limit ({max_vertices}), stopping")
                        break
    
    vertices = np.array(list(vertices_set))
    print(f"Loaded {len(vertices)} unique vertices from {vertex_count} total vertices")
    return vertices

def calculate_strike_angle_from_coords(vertices_km):
    """
    Calculate strike angle from coordinate distribution using PCA
    """
    # Use PCA on surface vertices to find principal direction
    surface_vertices = vertices_km[vertices_km[:, 2] > np.max(vertices_km[:, 2]) - 5]  # Top 5km
    
    if len(surface_vertices) < 10:
        surface_vertices = vertices_km  # Use all vertices if not enough surface vertices
    
    # Center the data
    centered = surface_vertices[:, :2] - np.mean(surface_vertices[:, :2], axis=0)  # Only x,y
    
    # Compute covariance matrix
    cov_matrix = np.cov(centered.T)
    
    # Get eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    
    # Strike is along the principal (longest) direction
    principal_vector = eigenvectors[:, np.argmax(eigenvalues)]
    
    # Calculate angle from x-axis
    strike_angle = np.degrees(np.arctan2(principal_vector[1], principal_vector[0]))
    
    return strike_angle

def calculate_planar_reference(vertices, method='best_fit_plane'):
    """
    Calculate a planar reference surface for smoothness comparison
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    method : str
        Method for calculating reference plane:
        - 'best_fit_plane': Least squares fit to all vertices
        - 'top_surface': Best fit to surface vertices only
        - 'pca_plane': Principal component analysis based plane
        
    Returns:
    --------
    plane_params : dict
        Dictionary containing plane parameters and reference surface info
    """
    vertices = np.array(vertices)
    
    if method == 'best_fit_plane':
        # Fit plane ax + by + cz + d = 0 to all vertices using SVD
        centroid = np.mean(vertices, axis=0)
        centered = vertices - centroid
        
        # SVD to find normal vector (smallest singular value direction)
        _, _, vh = np.linalg.svd(centered)
        normal = vh[-1, :]  # Last row is normal to best-fit plane
        
        # Plane equation: normal · (point - centroid) = 0
        # Or: ax + by + cz = d where (a,b,c) = normal, d = normal · centroid
        d = np.dot(normal, centroid)
        
        plane_params = {
            'method': method,
            'normal': normal,
            'centroid': centroid,
            'd': d,
            'equation': f"{normal[0]:.6f}x + {normal[1]:.6f}y + {normal[2]:.6f}z = {d:.6f}"
        }
        
    elif method == 'top_surface':
        # Use only surface vertices (top 10% by z-coordinate)
        z_threshold = np.percentile(vertices[:, 2], 90)
        surface_vertices = vertices[vertices[:, 2] >= z_threshold]
        
        if len(surface_vertices) < 10:
            surface_vertices = vertices  # Fallback to all vertices
            
        return calculate_planar_reference(surface_vertices, 'best_fit_plane')
        
    elif method == 'pca_plane':
        # Use PCA to find the plane of maximum variance
        pca = PCA(n_components=3)
        pca.fit(vertices)
        
        # The plane normal is the direction of minimum variance (3rd component)
        normal = pca.components_[2]  # Third principal component
        centroid = np.mean(vertices, axis=0)
        d = np.dot(normal, centroid)
        
        plane_params = {
            'method': method,
            'normal': normal,
            'centroid': centroid,
            'd': d,
            'equation': f"{normal[0]:.6f}x + {normal[1]:.6f}y + {normal[2]:.6f}z = {d:.6f}",
            'explained_variance_ratio': pca.explained_variance_ratio_
        }
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return plane_params

def calculate_distance_to_plane(vertices, plane_params):
    """
    Calculate perpendicular distance from vertices to reference plane
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    plane_params : dict
        Plane parameters from calculate_planar_reference
        
    Returns:
    --------
    distances : array_like, shape (n,)
        Perpendicular distances to plane (signed)
    """
    vertices = np.array(vertices)
    normal = plane_params['normal']
    d = plane_params['d']
    
    # Distance = |ax + by + cz - d| / sqrt(a² + b² + c²)
    # Since normal is unit vector, denominator = 1
    distances = np.dot(vertices, normal) - d
    
    return distances

def calculate_surface_roughness(vertices, neighborhood_size=5):
    """
    Calculate surface roughness metrics
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    neighborhood_size : int
        Number of nearest neighbors to consider for local roughness
        
    Returns:
    --------
    roughness_metrics : dict
        Dictionary containing various roughness measurements
    """
    vertices = np.array(vertices)
    n_points = len(vertices)
    
    if n_points < neighborhood_size + 1:
        neighborhood_size = max(1, n_points - 1)
    
    # Calculate distances between all pairs of points
    distances = cdist(vertices, vertices)
    
    # For each point, find its k nearest neighbors
    local_roughness = np.zeros(n_points)
    local_std = np.zeros(n_points)
    
    for i in range(n_points):
        # Get indices of k nearest neighbors (excluding self)
        neighbor_indices = np.argsort(distances[i])[1:neighborhood_size+1]
        neighbor_points = vertices[neighbor_indices]
        
        # Calculate local plane through neighbors
        if len(neighbor_points) >= 3:
            try:
                # Fit plane to neighbors
                centroid = np.mean(neighbor_points, axis=0)
                centered = neighbor_points - centroid
                _, _, vh = np.linalg.svd(centered)
                local_normal = vh[-1, :]
                
                # Distance from current point to local plane
                point_to_centroid = vertices[i] - centroid
                local_roughness[i] = abs(np.dot(point_to_centroid, local_normal))
                
                # Standard deviation of neighbor z-coordinates
                local_std[i] = np.std(neighbor_points[:, 2])
                
            except:
                # Fallback: use z-coordinate standard deviation
                local_roughness[i] = np.std(neighbor_points[:, 2])
                local_std[i] = local_roughness[i]
        else:
            # Not enough neighbors, use simple z variation
            if len(neighbor_points) > 0:
                local_roughness[i] = np.std(neighbor_points[:, 2])
                local_std[i] = local_roughness[i]
    
    roughness_metrics = {
        'mean_local_roughness': np.mean(local_roughness),
        'std_local_roughness': np.std(local_roughness),
        'max_local_roughness': np.max(local_roughness),
        'rms_roughness': np.sqrt(np.mean(local_roughness**2)),
        'mean_local_std': np.mean(local_std),
        'roughness_distribution': local_roughness,
        'neighborhood_size': neighborhood_size
    }
    
    return roughness_metrics

def calculate_curvature_metrics(vertices, grid_shape=None, resolution=1.0):
    """
    Calculate surface curvature metrics if vertices form a regular grid
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    grid_shape : tuple, optional
        Shape of the grid (nz, nx) if vertices form a regular grid
    resolution : float
        Grid resolution for derivative calculations
        
    Returns:
    --------
    curvature_metrics : dict
        Dictionary containing curvature measurements
    """
    vertices = np.array(vertices)
    
    if grid_shape is None:
        # Estimate if this could be a regular grid
        unique_x = len(np.unique(np.round(vertices[:, 0], 3)))
        unique_z = len(np.unique(np.round(vertices[:, 2], 3)))
        
        if unique_x * unique_z == len(vertices):
            grid_shape = (unique_z, unique_x)
        else:
            return {'error': 'Cannot determine regular grid structure for curvature calculation'}
    
    nz, nx = grid_shape
    
    # Reshape z-coordinates into grid
    try:
        # Sort vertices by z (descending) then x (ascending) to match grid structure
        sorted_indices = np.lexsort((vertices[:, 0], -vertices[:, 2]))
        sorted_vertices = vertices[sorted_indices]
        
        z_grid = sorted_vertices[:, 2].reshape(nz, nx)
        
        # Calculate second derivatives for curvature
        # d²z/dx² (curvature in x-direction)
        d2z_dx2 = np.zeros_like(z_grid)
        d2z_dx2[:, 1:-1] = (z_grid[:, 2:] - 2*z_grid[:, 1:-1] + z_grid[:, :-2]) / (resolution**2)
        
        # d²z/dz² (curvature in z-direction) 
        d2z_dz2 = np.zeros_like(z_grid)
        d2z_dz2[1:-1, :] = (z_grid[2:, :] - 2*z_grid[1:-1, :] + z_grid[:-2, :]) / (resolution**2)
        
        # Mixed derivative d²z/dxdz
        d2z_dxdz = np.zeros_like(z_grid)
        if nz > 2 and nx > 2:
            d2z_dxdz[1:-1, 1:-1] = ((z_grid[2:, 2:] - z_grid[2:, :-2]) - 
                                   (z_grid[:-2, 2:] - z_grid[:-2, :-2])) / (4 * resolution**2)
        
        # Mean curvature H = 0.5 * (κ1 + κ2)
        # Gaussian curvature K = κ1 * κ2
        # For surface z = f(x,y): H ≈ 0.5 * (fxx + fyy) in low-slope regions
        mean_curvature = 0.5 * (d2z_dx2 + d2z_dz2)
        gaussian_curvature = d2z_dx2 * d2z_dz2 - d2z_dxdz**2
        
        # Calculate curvature statistics (excluding boundary points)
        interior_mask = np.ones_like(mean_curvature, dtype=bool)
        interior_mask[0, :] = interior_mask[-1, :] = False
        interior_mask[:, 0] = interior_mask[:, -1] = False
        
        valid_mean_curv = mean_curvature[interior_mask]
        valid_gauss_curv = gaussian_curvature[interior_mask]
        
        curvature_metrics = {
            'mean_curvature_stats': {
                'mean': np.mean(valid_mean_curv),
                'std': np.std(valid_mean_curv),
                'rms': np.sqrt(np.mean(valid_mean_curv**2)),
                'max_abs': np.max(np.abs(valid_mean_curv))
            },
            'gaussian_curvature_stats': {
                'mean': np.mean(valid_gauss_curv),
                'std': np.std(valid_gauss_curv),
                'rms': np.sqrt(np.mean(valid_gauss_curv**2)),
                'max_abs': np.max(np.abs(valid_gauss_curv))
            },
            'grid_shape': grid_shape,
            'resolution': resolution,
            'interior_points': np.sum(interior_mask)
        }
        
    except Exception as e:
        curvature_metrics = {'error': f'Curvature calculation failed: {str(e)}'}
    
    return curvature_metrics

def calculate_smoothness_metrics(vertices, reference_method='best_fit_plane', 
                                neighborhood_size=5, grid_shape=None, resolution=1.0):
    """
    Calculate comprehensive smoothness metrics for a fault surface
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    reference_method : str
        Method for calculating reference plane ('best_fit_plane', 'top_surface', 'pca_plane')
    neighborhood_size : int
        Size of neighborhood for roughness calculation
    grid_shape : tuple, optional
        Grid shape for curvature analysis
    resolution : float
        Grid resolution
        
    Returns:
    --------
    smoothness_data : dict
        Comprehensive smoothness analysis results
    """
    vertices = np.array(vertices)
    
    print(f"Calculating smoothness metrics for {len(vertices)} vertices...")
    
    # 1. Calculate planar reference
    plane_params = calculate_planar_reference(vertices, method=reference_method)
    distances_to_plane = calculate_distance_to_plane(vertices, plane_params)
    
    # 2. Planar deviation metrics
    planar_metrics = {
        'rms_deviation': np.sqrt(np.mean(distances_to_plane**2)),
        'max_deviation': np.max(np.abs(distances_to_plane)),
        'mean_abs_deviation': np.mean(np.abs(distances_to_plane)),
        'std_deviation': np.std(distances_to_plane),
        'deviation_range': np.max(distances_to_plane) - np.min(distances_to_plane)
    }
    
    # 3. Surface roughness metrics
    roughness_metrics = calculate_surface_roughness(vertices, neighborhood_size)
    
    # 4. Curvature metrics (if applicable)
    curvature_metrics = calculate_curvature_metrics(vertices, grid_shape, resolution)
    
    # 5. Overall smoothness score (lower = smoother)
    # Combine multiple metrics into a single score
    smoothness_score = (planar_metrics['rms_deviation'] + 
                       roughness_metrics['rms_roughness'] + 
                       (curvature_metrics.get('mean_curvature_stats', {}).get('rms', 0) * 100))
    
    smoothness_data = {
        'planar_reference': plane_params,
        'planar_deviation_metrics': planar_metrics,
        'surface_roughness_metrics': roughness_metrics,
        'curvature_metrics': curvature_metrics,
        'distances_to_plane': distances_to_plane,
        'overall_smoothness_score': smoothness_score,
        'analysis_parameters': {
            'reference_method': reference_method,
            'neighborhood_size': neighborhood_size,
            'grid_shape': grid_shape,
            'resolution': resolution,
            'n_vertices': len(vertices)
        }
    }
    
    return smoothness_data

def compare_smoothness_preservation(original_vertices, processed_vertices, 
                                   original_grid_shape=None, processed_grid_shape=None, 
                                   resolution=1.0):
    """
    Compare smoothness between original STL and processed mesh
    
    Parameters:
    -----------
    original_vertices : array_like, shape (n, 3)
        Original STL vertices
    processed_vertices : array_like, shape (m, 3)  
        Processed mesh vertices
    original_grid_shape : tuple, optional
        Grid shape for original vertices
    processed_grid_shape : tuple, optional
        Grid shape for processed vertices
    resolution : float
        Grid resolution
        
    Returns:
    --------
    comparison_data : dict
        Detailed comparison of smoothness metrics
    """
    print("Comparing smoothness preservation between original and processed surfaces...")
    
    # Calculate smoothness for both surfaces
    original_smoothness = calculate_smoothness_metrics(
        original_vertices, grid_shape=original_grid_shape, resolution=resolution)
    
    processed_smoothness = calculate_smoothness_metrics(
        processed_vertices, grid_shape=processed_grid_shape, resolution=resolution)
    
    # Calculate preservation ratios
    orig_planar = original_smoothness['planar_deviation_metrics']
    proc_planar = processed_smoothness['planar_deviation_metrics']
    
    orig_rough = original_smoothness['surface_roughness_metrics']
    proc_rough = processed_smoothness['surface_roughness_metrics']
    
    preservation_ratios = {
        'rms_deviation_ratio': proc_planar['rms_deviation'] / orig_planar['rms_deviation'],
        'max_deviation_ratio': proc_planar['max_deviation'] / orig_planar['max_deviation'],
        'roughness_ratio': proc_rough['rms_roughness'] / orig_rough['rms_roughness'],
        'smoothness_score_ratio': (processed_smoothness['overall_smoothness_score'] / 
                                  original_smoothness['overall_smoothness_score'])
    }
    
    # Overall preservation assessment
    avg_preservation = np.mean(list(preservation_ratios.values()))
    preservation_quality = "excellent" if avg_preservation < 1.1 else \
                          "good" if avg_preservation < 1.5 else \
                          "moderate" if avg_preservation < 2.0 else "poor"
    
    comparison_data = {
        'original_smoothness': original_smoothness,
        'processed_smoothness': processed_smoothness,
        'preservation_ratios': preservation_ratios,
        'average_preservation_ratio': avg_preservation,
        'preservation_quality': preservation_quality,
        'summary': {
            'original_score': original_smoothness['overall_smoothness_score'],
            'processed_score': processed_smoothness['overall_smoothness_score'],
            'score_change': processed_smoothness['overall_smoothness_score'] - 
                           original_smoothness['overall_smoothness_score'],
            'relative_change_percent': ((processed_smoothness['overall_smoothness_score'] / 
                                        original_smoothness['overall_smoothness_score']) - 1) * 100
        }
    }
    
    return comparison_data

def apply_gaussian_smoothing(vertices, grid_shape, sigma=1.0, preserve_boundaries=True):
    """
    Apply Gaussian smoothing to surface vertices
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    grid_shape : tuple
        Shape of the grid (nz, nx)
    sigma : float
        Standard deviation for Gaussian kernel
    preserve_boundaries : bool
        Whether to preserve boundary values
        
    Returns:
    --------
    smoothed_vertices : array_like, shape (n, 3)
        Smoothed vertex coordinates
    smoothing_info : dict
        Information about the smoothing operation
    """
    vertices = np.array(vertices)
    nz, nx = grid_shape
    
    try:
        # Sort vertices by z (descending) then x (ascending) to match grid structure
        sorted_indices = np.lexsort((vertices[:, 0], -vertices[:, 2]))
        sorted_vertices = vertices[sorted_indices]
        
        # Extract coordinates and reshape into grids
        x_grid = sorted_vertices[:, 0].reshape(nz, nx)
        y_grid = sorted_vertices[:, 1].reshape(nz, nx)
        z_grid = sorted_vertices[:, 2].reshape(nz, nx)
        
        # Apply Gaussian smoothing to y-coordinates (fault surface position)
        if preserve_boundaries:
            # Create mask for interior points
            interior_mask = np.ones((nz, nx), dtype=bool)
            interior_mask[0, :] = interior_mask[-1, :] = False  # Top and bottom boundaries
            interior_mask[:, 0] = interior_mask[:, -1] = False  # Left and right boundaries
            
            # Smooth only interior points
            y_smoothed = y_grid.copy()
            interior_smoothed = gaussian_filter(y_grid[interior_mask], sigma=sigma)
            y_smoothed[interior_mask] = interior_smoothed
        else:
            # Smooth entire grid
            y_smoothed = gaussian_filter(y_grid, sigma=sigma)
        
        # Reconstruct smoothed vertices
        smoothed_vertices = np.column_stack([
            x_grid.ravel(), 
            y_smoothed.ravel(), 
            z_grid.ravel()
        ])
        
        # Restore original order
        restore_indices = np.argsort(sorted_indices)
        smoothed_vertices = smoothed_vertices[restore_indices]
        
        # Calculate smoothing metrics
        original_roughness = np.std(y_grid)
        smoothed_roughness = np.std(y_smoothed)
        smoothing_factor = original_roughness / smoothed_roughness if smoothed_roughness > 0 else np.inf
        
        smoothing_info = {
            'method': 'gaussian_smoothing',
            'sigma': sigma,
            'preserve_boundaries': preserve_boundaries,
            'original_roughness': original_roughness,
            'smoothed_roughness': smoothed_roughness,
            'smoothing_factor': smoothing_factor,
            'grid_shape': grid_shape
        }
        
    except Exception as e:
        print(f"Gaussian smoothing failed: {str(e)}")
        smoothed_vertices = vertices.copy()
        smoothing_info = {'error': str(e)}
    
    return smoothed_vertices, smoothing_info

def apply_adaptive_smoothing(vertices, target_roughness=0.1, max_iterations=10, 
                           grid_shape=None, neighborhood_size=5):
    """
    Apply adaptive smoothing to achieve target roughness level
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    target_roughness : float
        Target RMS roughness level
    max_iterations : int
        Maximum number of smoothing iterations
    grid_shape : tuple, optional
        Grid shape for structured smoothing
    neighborhood_size : int
        Size of neighborhood for roughness calculation
        
    Returns:
    --------
    smoothed_vertices : array_like, shape (n, 3)
        Adaptively smoothed vertices
    smoothing_info : dict
        Detailed information about the adaptive smoothing process
    """
    vertices = np.array(vertices)
    current_vertices = vertices.copy()
    
    smoothing_history = []
    
    for iteration in range(max_iterations):
        # Calculate current roughness
        roughness_metrics = calculate_surface_roughness(current_vertices, neighborhood_size)
        current_roughness = roughness_metrics['rms_roughness']
        
        print(f"Iteration {iteration + 1}: Current RMS roughness = {current_roughness:.6f}")
        
        smoothing_history.append({
            'iteration': iteration + 1,
            'roughness': current_roughness,
            'target_roughness': target_roughness
        })
        
        # Check if target is reached
        if current_roughness <= target_roughness:
            print(f"Target roughness {target_roughness:.6f} achieved in {iteration + 1} iterations")
            break
        
        # Determine smoothing parameters based on current roughness
        roughness_ratio = current_roughness / target_roughness
        if roughness_ratio > 10:
            sigma = 2.0  # Strong smoothing
        elif roughness_ratio > 5:
            sigma = 1.5  # Moderate smoothing
        elif roughness_ratio > 2:
            sigma = 1.0  # Mild smoothing
        else:
            sigma = 0.5  # Light smoothing
        
        # Apply smoothing
        if grid_shape is not None:
            current_vertices, step_info = apply_gaussian_smoothing(
                current_vertices, grid_shape, sigma=sigma, preserve_boundaries=True)
        else:
            # Use RBF-based smoothing for unstructured data
            current_vertices, step_info = apply_rbf_smoothing(
                current_vertices, sigma=sigma)
        
        smoothing_history[-1]['sigma'] = sigma
        smoothing_history[-1]['step_info'] = step_info
        
        # Prevent over-smoothing
        if iteration > 0 and abs(smoothing_history[-1]['roughness'] - 
                                smoothing_history[-2]['roughness']) < 1e-6:
            print("Convergence reached (minimal roughness change)")
            break
    
    final_roughness = calculate_surface_roughness(current_vertices, neighborhood_size)['rms_roughness']
    
    smoothing_info = {
        'method': 'adaptive_smoothing',
        'target_roughness': target_roughness,
        'final_roughness': final_roughness,
        'iterations_used': len(smoothing_history),
        'max_iterations': max_iterations,
        'target_achieved': final_roughness <= target_roughness,
        'smoothing_history': smoothing_history,
        'roughness_reduction_factor': vertices.std() / current_vertices.std() if current_vertices.std() > 0 else 1.0
    }
    
    return current_vertices, smoothing_info

def apply_rbf_smoothing(vertices, sigma=1.0, epsilon=None, function='thin_plate_spline'):
    """
    Apply RBF-based smoothing for unstructured vertex data
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    sigma : float
        Smoothing parameter
    epsilon : float, optional
        Shape parameter for RBF
    function : str
        RBF function type
        
    Returns:
    --------
    smoothed_vertices : array_like, shape (n, 3)
        RBF-smoothed vertices
    smoothing_info : dict
        Information about RBF smoothing
    """
    vertices = np.array(vertices)
    
    try:
        # Use subset of points for RBF fitting to improve performance
        n_points = len(vertices)
        if n_points > 1000:
            # Sample points for RBF fitting
            sample_indices = np.random.choice(n_points, size=min(1000, n_points), replace=False)
            sample_vertices = vertices[sample_indices]
        else:
            sample_vertices = vertices
        
        # Create RBF interpolator for y-coordinates (surface position)
        if epsilon is None:
            # Auto-determine epsilon based on point spacing
            distances = cdist(sample_vertices[:, [0, 2]], sample_vertices[:, [0, 2]])
            epsilon = np.median(distances[distances > 0]) * sigma
        
        # Fit RBF to sampled data
        rbf = RBFInterpolator(
            sample_vertices[:, [0, 2]],  # x, z coordinates
            sample_vertices[:, 1],        # y coordinates
            kernel=function,
            epsilon=epsilon,
            smoothing=sigma**2  # Smoothing parameter
        )
        
        # Apply RBF to all vertices
        y_smoothed = rbf(vertices[:, [0, 2]])
        
        # Create smoothed vertices
        smoothed_vertices = vertices.copy()
        smoothed_vertices[:, 1] = y_smoothed
        
        # Calculate smoothing metrics
        original_roughness = np.std(vertices[:, 1])
        smoothed_roughness = np.std(y_smoothed)
        smoothing_factor = original_roughness / smoothed_roughness if smoothed_roughness > 0 else np.inf
        
        smoothing_info = {
            'method': 'rbf_smoothing',
            'function': function,
            'sigma': sigma,
            'epsilon': epsilon,
            'n_sample_points': len(sample_vertices),
            'original_roughness': original_roughness,
            'smoothed_roughness': smoothed_roughness,
            'smoothing_factor': smoothing_factor
        }
        
    except Exception as e:
        print(f"RBF smoothing failed: {str(e)}")
        smoothed_vertices = vertices.copy()
        smoothing_info = {'error': str(e)}
    
    return smoothed_vertices, smoothing_info

def apply_constrained_smoothing(vertices, constraints=None, smoothing_weight=1.0, 
                              constraint_weight=10.0, grid_shape=None):
    """
    Apply smoothing while preserving specific geometric constraints
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    constraints : dict, optional
        Dictionary specifying constraints:
        - 'fixed_points': indices of vertices to keep fixed
        - 'target_values': target y-values for specific vertices
        - 'gradient_constraints': target gradients at specific locations
    smoothing_weight : float
        Weight for smoothness term in optimization
    constraint_weight : float
        Weight for constraint terms in optimization
    grid_shape : tuple, optional
        Grid shape for structured processing
        
    Returns:
    --------
    smoothed_vertices : array_like, shape (n, 3)
        Constraint-preserving smoothed vertices
    smoothing_info : dict
        Information about constrained smoothing
    """
    vertices = np.array(vertices)
    
    if constraints is None:
        constraints = {}
    
    # Extract y-coordinates for optimization
    y_original = vertices[:, 1].copy()
    y_optimized = y_original.copy()
    
    # Define objective function
    def objective(y_vals):
        total_cost = 0.0
        
        # Smoothness term: minimize second derivatives
        if grid_shape is not None:
            nz, nx = grid_shape
            try:
                y_grid = y_vals.reshape(nz, nx)
                
                # Second derivatives in x and z directions
                d2y_dx2 = np.zeros_like(y_grid)
                d2y_dz2 = np.zeros_like(y_grid)
                
                # Calculate second derivatives
                d2y_dx2[:, 1:-1] = y_grid[:, 2:] - 2*y_grid[:, 1:-1] + y_grid[:, :-2]
                d2y_dz2[1:-1, :] = y_grid[2:, :] - 2*y_grid[1:-1, :] + y_grid[:-2, :]
                
                # Smoothness cost
                smoothness_cost = np.sum(d2y_dx2**2) + np.sum(d2y_dz2**2)
                total_cost += smoothing_weight * smoothness_cost
                
            except:
                # Fallback: minimize differences between adjacent points
                smoothness_cost = np.sum(np.diff(y_vals)**2)
                total_cost += smoothing_weight * smoothness_cost
        else:
            # For unstructured data, minimize local variations
            smoothness_cost = np.var(y_vals)
            total_cost += smoothing_weight * smoothness_cost
        
        # Constraint terms
        if 'fixed_points' in constraints:
            fixed_indices = constraints['fixed_points']
            fixed_cost = np.sum((y_vals[fixed_indices] - y_original[fixed_indices])**2)
            total_cost += constraint_weight * fixed_cost
        
        if 'target_values' in constraints:
            target_indices = constraints['target_values']['indices']
            target_values = constraints['target_values']['values']
            target_cost = np.sum((y_vals[target_indices] - target_values)**2)
            total_cost += constraint_weight * target_cost
        
        return total_cost
    
    try:
        # Optimize
        result = minimize(objective, y_optimized, method='L-BFGS-B')
        
        if result.success:
            y_smoothed = result.x
            
            # Create smoothed vertices
            smoothed_vertices = vertices.copy()
            smoothed_vertices[:, 1] = y_smoothed
            
            # Calculate metrics
            original_roughness = np.std(y_original)
            smoothed_roughness = np.std(y_smoothed)
            
            smoothing_info = {
                'method': 'constrained_smoothing',
                'optimization_success': True,
                'iterations': result.nit,
                'final_cost': result.fun,
                'smoothing_weight': smoothing_weight,
                'constraint_weight': constraint_weight,
                'original_roughness': original_roughness,
                'smoothed_roughness': smoothed_roughness,
                'constraints_applied': list(constraints.keys()),
                'smoothing_factor': original_roughness / smoothed_roughness if smoothed_roughness > 0 else np.inf
            }
        else:
            print(f"Constrained optimization failed: {result.message}")
            smoothed_vertices = vertices.copy()
            smoothing_info = {'error': f'Optimization failed: {result.message}'}
            
    except Exception as e:
        print(f"Constrained smoothing failed: {str(e)}")
        smoothed_vertices = vertices.copy()
        smoothing_info = {'error': str(e)}
    
    return smoothed_vertices, smoothing_info

def apply_smoothing_control(vertices, smoothing_method='gaussian', 
                          smoothing_params=None, grid_shape=None):
    """
    Apply specified smoothing method with given parameters
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices
    smoothing_method : str
        Smoothing method: 'gaussian', 'adaptive', 'rbf', 'constrained', 'none'
    smoothing_params : dict
        Parameters specific to each smoothing method
    grid_shape : tuple, optional
        Grid shape for structured methods
        
    Returns:
    --------
    smoothed_vertices : array_like, shape (n, 3)
        Smoothed vertices
    smoothing_info : dict
        Detailed smoothing information
    """
    if smoothing_params is None:
        smoothing_params = {}
        
    print(f"Applying {smoothing_method} smoothing...")
    
    if smoothing_method == 'gaussian':
        sigma = smoothing_params.get('sigma', 1.0)
        preserve_boundaries = smoothing_params.get('preserve_boundaries', True)
        
        if grid_shape is not None:
            return apply_gaussian_smoothing(vertices, grid_shape, sigma, preserve_boundaries)
        else:
            return apply_rbf_smoothing(vertices, sigma=sigma)
            
    elif smoothing_method == 'adaptive':
        target_roughness = smoothing_params.get('target_roughness', 0.1)
        max_iterations = smoothing_params.get('max_iterations', 10)
        neighborhood_size = smoothing_params.get('neighborhood_size', 5)
        
        return apply_adaptive_smoothing(vertices, target_roughness, max_iterations, 
                                      grid_shape, neighborhood_size)
                                      
    elif smoothing_method == 'rbf':
        sigma = smoothing_params.get('sigma', 1.0)
        epsilon = smoothing_params.get('epsilon', None)
        function = smoothing_params.get('function', 'thin_plate_spline')
        
        return apply_rbf_smoothing(vertices, sigma, epsilon, function)
        
    elif smoothing_method == 'constrained':
        constraints = smoothing_params.get('constraints', {})
        smoothing_weight = smoothing_params.get('smoothing_weight', 1.0)
        constraint_weight = smoothing_params.get('constraint_weight', 10.0)
        
        return apply_constrained_smoothing(vertices, constraints, smoothing_weight, 
                                        constraint_weight, grid_shape)
                                        
    elif smoothing_method == 'none':
        return vertices.copy(), {'method': 'none', 'no_smoothing_applied': True}
        
    else:
        raise ValueError(f"Unknown smoothing method: {smoothing_method}")

def plot_smoothing_effects(original_vertices, smoothed_vertices, smoothing_info, output_prefix):
    """
    Create plots showing the effects of smoothing
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(18, 12))
    
    # Plot 1: Original surface
    ax1 = fig.add_subplot(231, projection='3d')
    ax1.scatter(original_vertices[:, 0], original_vertices[:, 1], original_vertices[:, 2],
                c=original_vertices[:, 1], cmap='viridis', s=1, alpha=0.7)
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_zlabel('Z (km)')
    ax1.set_title('Original Surface')
    ax1.set_box_aspect([1,1,1])
    
    # Plot 2: Smoothed surface
    ax2 = fig.add_subplot(232, projection='3d')
    ax2.scatter(smoothed_vertices[:, 0], smoothed_vertices[:, 1], smoothed_vertices[:, 2],
                c=smoothed_vertices[:, 1], cmap='viridis', s=1, alpha=0.7)
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Y (km)')
    ax2.set_zlabel('Z (km)')
    ax2.set_title(f'Smoothed Surface ({smoothing_info["method"]})')
    ax2.set_box_aspect([1,1,1])
    
    # Plot 3: Difference map
    ax3 = fig.add_subplot(233, projection='3d')
    y_diff = smoothed_vertices[:, 1] - original_vertices[:, 1]
    scatter = ax3.scatter(original_vertices[:, 0], original_vertices[:, 1], original_vertices[:, 2],
                         c=y_diff, cmap='RdBu_r', s=1, alpha=0.7)
    ax3.set_xlabel('X (km)')
    ax3.set_ylabel('Y (km)')
    ax3.set_zlabel('Z (km)')
    ax3.set_title('Smoothing Effect (Y-difference)')
    ax3.set_box_aspect([1,1,1])
    plt.colorbar(scatter, ax=ax3, shrink=0.6, label='Y change (km)')
    
    # Plot 4: Y-coordinate histograms
    ax4 = fig.add_subplot(234)
    ax4.hist(original_vertices[:, 1], bins=50, alpha=0.7, label='Original', density=True)
    ax4.hist(smoothed_vertices[:, 1], bins=50, alpha=0.7, label='Smoothed', density=True)
    ax4.set_xlabel('Y coordinate (km)')
    ax4.set_ylabel('Density')
    ax4.set_title('Y-coordinate Distribution')
    ax4.legend()
    
    # Plot 5: Smoothing difference histogram
    ax5 = fig.add_subplot(235)
    ax5.hist(y_diff, bins=50, alpha=0.7, edgecolor='black')
    ax5.set_xlabel('Y change (km)')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Distribution of Smoothing Changes')
    ax5.axvline(0, color='red', linestyle='--', alpha=0.7)
    
    # Plot 6: Summary statistics
    ax6 = fig.add_subplot(236)
    ax6.axis('off')
    
    # Create summary text
    summary_text = f"""Smoothing Summary
    
Method: {smoothing_info['method']}

Original Surface:
• Y std: {np.std(original_vertices[:, 1]):.6f} km
• Y range: {np.ptp(original_vertices[:, 1]):.6f} km

Smoothed Surface:
• Y std: {np.std(smoothed_vertices[:, 1]):.6f} km  
• Y range: {np.ptp(smoothed_vertices[:, 1]):.6f} km

Smoothing Effects:
• Mean Y change: {np.mean(y_diff):.6f} km
• RMS Y change: {np.sqrt(np.mean(y_diff**2)):.6f} km
• Max |Y change|: {np.max(np.abs(y_diff)):.6f} km"""

    if 'smoothing_factor' in smoothing_info:
        summary_text += f"\n• Smoothing factor: {smoothing_info['smoothing_factor']:.3f}"
    
    if smoothing_info['method'] == 'adaptive_smoothing':
        summary_text += f"\n• Target roughness: {smoothing_info['target_roughness']:.6f}"
        summary_text += f"\n• Final roughness: {smoothing_info['final_roughness']:.6f}"
        summary_text += f"\n• Iterations: {smoothing_info['iterations_used']}"
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    plot_filename = f"{output_prefix}_smoothing_effects.png"
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Smoothing effects plot saved to {plot_filename}")
    plt.close()

def calculate_vertical_plane_target(vertices, method='centroid_line'):
    """
    Calculate target vertical plane geometry for transition
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates  
    method : str
        Method for determining vertical plane:
        - 'centroid_line': Vertical line through surface centroid
        - 'best_fit_line': Best-fit vertical line through all points
        - 'strike_aligned': Vertical plane aligned with strike direction
        
    Returns:
    --------
    plane_info : dict
        Information about the target vertical plane
    """
    vertices = np.array(vertices)
    
    if method == 'centroid_line':
        # Use surface centroid projected vertically
        surface_vertices = vertices[vertices[:, 2] > np.percentile(vertices[:, 2], 75)]
        if len(surface_vertices) < 10:
            surface_vertices = vertices
        
        centroid_x = np.mean(surface_vertices[:, 0])
        centroid_y = np.mean(surface_vertices[:, 1])
        
        plane_info = {
            'method': method,
            'type': 'vertical_line',
            'x_position': centroid_x,
            'y_position': centroid_y,
            'description': f'Vertical line at x={centroid_x:.3f}, y={centroid_y:.3f}'
        }
        
    elif method == 'best_fit_line':
        # Find best-fit line through y-positions at each x
        unique_x = np.unique(np.round(vertices[:, 0], 3))
        if len(unique_x) > 3:
            x_centers = []
            y_means = []
            
            for x in unique_x:
                x_mask = np.abs(vertices[:, 0] - x) < 0.1
                if np.sum(x_mask) > 0:
                    x_centers.append(x)
                    y_means.append(np.mean(vertices[x_mask, 1]))
            
            if len(x_centers) > 2:
                # Fit line y = mx + b
                coeffs = np.polyfit(x_centers, y_means, 1)
                slope, intercept = coeffs
                
                plane_info = {
                    'method': method,
                    'type': 'vertical_plane',
                    'slope': slope,
                    'intercept': intercept,
                    'description': f'Vertical plane: y = {slope:.6f}*x + {intercept:.3f}'
                }
            else:
                # Fallback to centroid
                return calculate_vertical_plane_target(vertices, 'centroid_line')
        else:
            # Fallback to centroid
            return calculate_vertical_plane_target(vertices, 'centroid_line')
            
    elif method == 'strike_aligned':
        # Use strike direction to define vertical plane
        strike_angle = calculate_strike_angle_from_coords(vertices)
        
        # Calculate mean position
        mean_x = np.mean(vertices[:, 0])
        mean_y = np.mean(vertices[:, 1])
        
        # Strike-aligned vertical plane
        strike_rad = np.radians(strike_angle)
        normal_x = -np.sin(strike_rad)  # Normal to strike direction
        normal_y = np.cos(strike_rad)
        
        plane_info = {
            'method': method,
            'type': 'strike_aligned_plane',
            'strike_angle': strike_angle,
            'normal_x': normal_x,
            'normal_y': normal_y,
            'mean_x': mean_x,
            'mean_y': mean_y,
            'description': f'Strike-aligned vertical plane (strike={strike_angle:.1f}°)'
        }
        
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return plane_info

def calculate_y_position_on_vertical_plane(x_coords, z_coords, plane_info):
    """
    Calculate y-positions on the target vertical plane
    
    Parameters:
    -----------
    x_coords : array_like
        X coordinates
    z_coords : array_like  
        Z coordinates (not used for vertical planes, but kept for consistency)
    plane_info : dict
        Vertical plane information from calculate_vertical_plane_target
        
    Returns:
    --------
    y_coords : array_like
        Y coordinates on the vertical plane
    """
    x_coords = np.asarray(x_coords)
    
    if plane_info['type'] == 'vertical_line':
        # All points project to the same y-position
        y_coords = np.full_like(x_coords, plane_info['y_position'])
        
    elif plane_info['type'] == 'vertical_plane':
        # Linear relationship: y = mx + b
        y_coords = plane_info['slope'] * x_coords + plane_info['intercept']
        
    elif plane_info['type'] == 'strike_aligned_plane':
        # Project onto strike-aligned plane
        # Distance from mean position in strike-normal direction should be zero
        mean_x = plane_info['mean_x']
        mean_y = plane_info['mean_y']
        normal_x = plane_info['normal_x']
        normal_y = plane_info['normal_y']
        
        # Project points onto the plane
        # For vertical plane: find y such that (x-mean_x)*normal_x + (y-mean_y)*normal_y = 0
        # Solving for y: y = mean_y - (x-mean_x)*normal_x/normal_y
        if abs(normal_y) > 1e-10:
            y_coords = mean_y - (x_coords - mean_x) * normal_x / normal_y
        else:
            # Nearly horizontal normal, use constant y
            y_coords = np.full_like(x_coords, mean_y)
    
    else:
        raise ValueError(f"Unknown plane type: {plane_info['type']}")
    
    return y_coords

def apply_smooth_vertical_transition(vertices, grid_shape, transition_start_depth, 
                                   transition_distance=5.0, plane_method='best_fit_line',
                                   blend_function='cosine'):
    """
    Apply smooth transition from complex fault geometry to vertical plane
    
    Parameters:
    -----------
    vertices : array_like, shape (n, 3)
        Input vertices as [x, y, z] coordinates
    grid_shape : tuple
        Shape of the grid (nz, nx)
    transition_start_depth : float
        Depth (negative z) where transition begins (km)
    transition_distance : float
        Distance over which transition occurs (km)
    plane_method : str
        Method for calculating target vertical plane
    blend_function : str
        Blending function: 'linear', 'cosine', 'smooth_step', 'exponential'
        
    Returns:
    --------
    transitioned_vertices : array_like, shape (n, 3)
        Vertices with smooth vertical transition applied
    transition_info : dict
        Information about the transition
    """
    vertices = np.array(vertices)
    nz, nx = grid_shape
    
    # Calculate target vertical plane
    plane_info = calculate_vertical_plane_target(vertices, plane_method)
    
    # Transition parameters
    transition_end_depth = transition_start_depth - transition_distance
    
    print(f"Applying smooth vertical transition:")
    print(f"  Start depth: {transition_start_depth:.1f} km")
    print(f"  End depth: {transition_end_depth:.1f} km") 
    print(f"  Transition distance: {transition_distance:.1f} km")
    print(f"  Target plane: {plane_info['description']}")
    print(f"  Blend function: {blend_function}")
    
    try:
        # Sort vertices by z (descending) then x (ascending) to match grid structure
        sorted_indices = np.lexsort((vertices[:, 0], -vertices[:, 2]))
        sorted_vertices = vertices[sorted_indices]
        
        # Extract coordinates and reshape into grids
        x_grid = sorted_vertices[:, 0].reshape(nz, nx)
        y_grid = sorted_vertices[:, 1].reshape(nz, nx)
        z_grid = sorted_vertices[:, 2].reshape(nz, nx)
        
        # Calculate target y-positions on vertical plane
        y_target_grid = np.zeros_like(y_grid)
        for i in range(nz):
            for j in range(nx):
                y_target_grid[i, j] = calculate_y_position_on_vertical_plane(
                    x_grid[i, j], z_grid[i, j], plane_info)
        
        # Apply transition blending
        y_transitioned = y_grid.copy()
        
        for i in range(nz):
            for j in range(nx):
                z = z_grid[i, j]
                
                if z >= transition_start_depth:
                    # Above transition zone - keep original geometry
                    continue
                elif z <= transition_end_depth:
                    # Below transition zone - use vertical plane
                    y_transitioned[i, j] = y_target_grid[i, j]
                else:
                    # Within transition zone - blend between original and target
                    # Calculate blend factor (0 = original, 1 = target)
                    progress = (transition_start_depth - z) / transition_distance
                    
                    if blend_function == 'linear':
                        blend_factor = progress
                    elif blend_function == 'cosine':
                        blend_factor = 0.5 * (1 - np.cos(np.pi * progress))
                    elif blend_function == 'smooth_step':
                        blend_factor = progress * progress * (3 - 2 * progress)
                    elif blend_function == 'exponential':
                        # Exponential decay from original to target
                        blend_factor = 1 - np.exp(-3 * progress)
                    else:
                        blend_factor = progress  # Default to linear
                    
                    # Blend between original and target
                    y_original = y_grid[i, j]
                    y_target = y_target_grid[i, j]
                    y_transitioned[i, j] = y_original * (1 - blend_factor) + y_target * blend_factor
        
        # Reconstruct transitioned vertices
        transitioned_vertices = np.column_stack([
            x_grid.ravel(),
            y_transitioned.ravel(),
            z_grid.ravel()
        ])
        
        # Restore original order
        restore_indices = np.argsort(sorted_indices)
        transitioned_vertices = transitioned_vertices[restore_indices]
        
        # Calculate transition metrics
        y_change = y_transitioned - y_grid
        max_change = np.max(np.abs(y_change))
        rms_change = np.sqrt(np.mean(y_change**2))
        transition_mask = (z_grid >= transition_end_depth) & (z_grid <= transition_start_depth)
        transition_points = np.sum(transition_mask)
        
        transition_info = {
            'method': 'smooth_vertical_transition',
            'plane_info': plane_info,
            'transition_start_depth': transition_start_depth,
            'transition_end_depth': transition_end_depth,
            'transition_distance': transition_distance,
            'blend_function': blend_function,
            'max_y_change': max_change,
            'rms_y_change': rms_change,
            'transition_points': transition_points,
            'total_points': len(vertices),
            'transition_fraction': transition_points / len(vertices)
        }
        
        print(f"Transition applied successfully:")
        print(f"  Points in transition zone: {transition_points}")
        print(f"  Max Y change: {max_change:.6f} km")
        print(f"  RMS Y change: {rms_change:.6f} km")
        
    except Exception as e:
        print(f"Vertical transition failed: {str(e)}")
        transitioned_vertices = vertices.copy()
        transition_info = {'error': str(e)}
    
    return transitioned_vertices, transition_info

def plot_vertical_transition_effects(original_vertices, transitioned_vertices, 
                                   transition_info, output_prefix):
    """
    Create plots showing the effects of vertical transition
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(18, 12))
    
    # Plot 1: Original surface
    ax1 = fig.add_subplot(231, projection='3d')
    ax1.scatter(original_vertices[:, 0], original_vertices[:, 1], original_vertices[:, 2],
                c=original_vertices[:, 2], cmap='viridis', s=1, alpha=0.7)
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_zlabel('Z (km)')
    ax1.set_title('Original Surface')
    ax1.set_box_aspect([1,1,1])
    
    # Plot 2: Transitioned surface
    ax2 = fig.add_subplot(232, projection='3d')
    ax2.scatter(transitioned_vertices[:, 0], transitioned_vertices[:, 1], transitioned_vertices[:, 2],
                c=transitioned_vertices[:, 2], cmap='viridis', s=1, alpha=0.7)
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Y (km)')
    ax2.set_zlabel('Z (km)')
    ax2.set_title('Surface with Vertical Transition')
    ax2.set_box_aspect([1,1,1])
    
    # Plot 3: Transition zone and change
    ax3 = fig.add_subplot(233, projection='3d')
    y_diff = transitioned_vertices[:, 1] - original_vertices[:, 1]
    scatter = ax3.scatter(original_vertices[:, 0], original_vertices[:, 1], original_vertices[:, 2],
                         c=y_diff, cmap='RdBu_r', s=1, alpha=0.7)
    ax3.set_xlabel('X (km)')
    ax3.set_ylabel('Y (km)')
    ax3.set_zlabel('Z (km)')
    ax3.set_title('Transition Effect (Y-difference)')
    ax3.set_box_aspect([1,1,1])
    plt.colorbar(scatter, ax=ax3, shrink=0.6, label='Y change (km)')
    
    # Add transition zone boundaries
    if 'error' not in transition_info:
        start_depth = transition_info['transition_start_depth']
        end_depth = transition_info['transition_end_depth']
        
        # Add horizontal planes to show transition boundaries
        x_range = [np.min(original_vertices[:, 0]), np.max(original_vertices[:, 0])]
        y_range = [np.min(original_vertices[:, 1]), np.max(original_vertices[:, 1])]
        xx, yy = np.meshgrid(x_range, y_range)
        
        ax3.plot_surface(xx, yy, np.full_like(xx, start_depth), alpha=0.3, color='green')
        ax3.plot_surface(xx, yy, np.full_like(xx, end_depth), alpha=0.3, color='red')
    
    # Plot 4: Y-position vs depth profile
    ax4 = fig.add_subplot(234)
    
    # Sample profiles at different x-positions
    x_positions = np.percentile(original_vertices[:, 0], [25, 50, 75])
    
    for i, x_pos in enumerate(x_positions):
        # Find points near this x-position
        x_mask = np.abs(original_vertices[:, 0] - x_pos) < 5.0  # Within 5km
        if np.sum(x_mask) > 10:
            orig_subset = original_vertices[x_mask]
            trans_subset = transitioned_vertices[x_mask]
            
            # Sort by depth
            z_sort = np.argsort(orig_subset[:, 2])
            
            ax4.plot(orig_subset[z_sort, 1], orig_subset[z_sort, 2], 
                    label=f'Original (x≈{x_pos:.0f}km)', alpha=0.7, linestyle='-')
            ax4.plot(trans_subset[z_sort, 1], trans_subset[z_sort, 2], 
                    label=f'Transitioned (x≈{x_pos:.0f}km)', alpha=0.7, linestyle='--')
    
    ax4.set_xlabel('Y Position (km)')
    ax4.set_ylabel('Depth (km)')
    ax4.set_title('Y-Position vs Depth Profiles')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    # Add transition zone
    if 'error' not in transition_info:
        ax4.axhline(transition_info['transition_start_depth'], color='green', 
                   linestyle=':', alpha=0.7, label='Transition start')
        ax4.axhline(transition_info['transition_end_depth'], color='red', 
                   linestyle=':', alpha=0.7, label='Transition end')
    
    # Plot 5: Transition change histogram
    ax5 = fig.add_subplot(235)
    ax5.hist(y_diff, bins=50, alpha=0.7, edgecolor='black')
    ax5.set_xlabel('Y change (km)')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Distribution of Transition Changes')
    ax5.axvline(0, color='red', linestyle='--', alpha=0.7)
    
    # Plot 6: Summary information
    ax6 = fig.add_subplot(236)
    ax6.axis('off')
    
    if 'error' not in transition_info:
        plane_info = transition_info['plane_info']
        summary_text = f"""Vertical Transition Summary

Target Plane: {plane_info['description']}
Method: {plane_info['method']}

Transition Parameters:
• Start depth: {transition_info['transition_start_depth']:.1f} km
• End depth: {transition_info['transition_end_depth']:.1f} km  
• Distance: {transition_info['transition_distance']:.1f} km
• Blend function: {transition_info['blend_function']}

Effects:
• Points in transition: {transition_info['transition_points']}
• Transition fraction: {transition_info['transition_fraction']:.1%}
• Max Y change: {transition_info['max_y_change']:.6f} km
• RMS Y change: {transition_info['rms_y_change']:.6f} km

Original Surface:
• Y std: {np.std(original_vertices[:, 1]):.6f} km
• Y range: {np.ptp(original_vertices[:, 1]):.6f} km

Transitioned Surface:
• Y std: {np.std(transitioned_vertices[:, 1]):.6f} km
• Y range: {np.ptp(transitioned_vertices[:, 1]):.6f} km"""
    else:
        summary_text = f"Vertical Transition Failed:\n{transition_info['error']}"
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    plot_filename = f"{output_prefix}_vertical_transition.png"
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Vertical transition plot saved to {plot_filename}")
    plt.close()

def create_output_folder(stl_filename, processing_params):
    """
    Create organized output folder with descriptive name
    
    Parameters:
    -----------
    stl_filename : str
        Path to STL file
    processing_params : dict
        Dictionary of processing parameters
        
    Returns:
    --------
    folder_path : str
        Path to created output folder
    """
    # Extract base filename without extension
    base_name = os.path.splitext(os.path.basename(stl_filename))[0]
    
    # Build folder name with key features
    folder_components = [base_name]
    
    # Add resolution
    resolution = processing_params.get('resolution', 1.0)
    folder_components.append(f"res{resolution:g}km")
    
    # Add smoothing info
    smoothing_method = processing_params.get('smoothing_method', 'none')
    if smoothing_method != 'none':
        folder_components.append(f"smooth_{smoothing_method}")
        
        smoothing_params = processing_params.get('smoothing_params', {})
        if 'sigma' in smoothing_params:
            folder_components.append(f"sigma{smoothing_params['sigma']:g}")
        if 'target_roughness' in smoothing_params:
            folder_components.append(f"rough{smoothing_params['target_roughness']:g}")
    
    # Add vertical transition info
    if processing_params.get('enable_vertical_transition', False):
        transition_start = processing_params.get('transition_start_depth', None)
        transition_dist = processing_params.get('transition_distance', 5.0)
        plane_method = processing_params.get('plane_method', 'best_fit_line')
        
        folder_components.append("vertical_transition")
        if transition_start is not None:
            folder_components.append(f"start{abs(transition_start):g}km")
        else:
            folder_components.append("start_auto")
        folder_components.append(f"dist{transition_dist:g}km")
        folder_components.append(f"{plane_method}")
    
    # Create folder name
    folder_name = ".".join(folder_components)
    
    # Create the folder
    folder_path = os.path.join(os.getcwd(), folder_name)
    os.makedirs(folder_path, exist_ok=True)
    
    print(f"Created output folder: {folder_name}")
    return folder_path

def plot_smoothness_analysis(smoothness_data, vertices, output_prefix, title_suffix=""):
    """
    Create visualizations for smoothness analysis
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    vertices = np.array(vertices)
    distances = smoothness_data['distances_to_plane']
    plane_params = smoothness_data['planar_reference']
    
    fig = plt.figure(figsize=(15, 10))
    
    # Plot 1: 3D surface with deviation from plane
    ax1 = fig.add_subplot(221, projection='3d')
    scatter = ax1.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], 
                         c=distances, cmap='RdBu_r', s=1, alpha=0.7)
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_zlabel('Z (km)')
    ax1.set_title(f'Surface Deviation from Reference Plane{title_suffix}')
    ax1.set_box_aspect([1,1,1])
    plt.colorbar(scatter, ax=ax1, shrink=0.6, label='Distance to plane (km)')
    
    # Plot 2: Histogram of deviations
    ax2 = fig.add_subplot(222)
    ax2.hist(distances, bins=50, alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Distance to Reference Plane (km)')
    ax2.set_ylabel('Frequency')
    ax2.set_title(f'Distribution of Deviations{title_suffix}')
    ax2.axvline(0, color='red', linestyle='--', alpha=0.7, label='Reference Plane')
    ax2.legend()
    
    # Plot 3: Surface roughness distribution
    ax3 = fig.add_subplot(223)
    if 'roughness_distribution' in smoothness_data['surface_roughness_metrics']:
        roughness_dist = smoothness_data['surface_roughness_metrics']['roughness_distribution']
        ax3.hist(roughness_dist, bins=50, alpha=0.7, edgecolor='black', color='orange')
        ax3.set_xlabel('Local Roughness')
        ax3.set_ylabel('Frequency')
        ax3.set_title(f'Local Surface Roughness Distribution{title_suffix}')
    
    # Plot 4: Summary metrics
    ax4 = fig.add_subplot(224)
    ax4.axis('off')
    
    # Create summary text
    planar_metrics = smoothness_data['planar_deviation_metrics']
    roughness_metrics = smoothness_data['surface_roughness_metrics']
    
    summary_text = f"""Smoothness Analysis Summary{title_suffix}
    
Planar Deviation Metrics:
• RMS Deviation: {planar_metrics['rms_deviation']:.4f} km
• Max Deviation: {planar_metrics['max_deviation']:.4f} km
• Mean Abs Dev: {planar_metrics['mean_abs_deviation']:.4f} km
• Standard Dev: {planar_metrics['std_deviation']:.4f} km

Surface Roughness Metrics:
• RMS Roughness: {roughness_metrics['rms_roughness']:.4f} km
• Mean Local Roughness: {roughness_metrics['mean_local_roughness']:.4f} km
• Max Local Roughness: {roughness_metrics['max_local_roughness']:.4f} km

Overall Smoothness Score: {smoothness_data['overall_smoothness_score']:.4f}

Reference Plane: {plane_params['method']}
{plane_params['equation']}"""
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    plot_filename = f"{output_prefix}_smoothness_analysis.png"
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Smoothness analysis plot saved to {plot_filename}")
    plt.close()

def plot_smoothness_comparison(comparison_data, output_prefix):
    """
    Create comparison plots for smoothness preservation
    """
    import matplotlib.pyplot as plt
    
    orig_data = comparison_data['original_smoothness']
    proc_data = comparison_data['processed_smoothness']
    ratios = comparison_data['preservation_ratios']
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # Plot 1: Deviation comparison
    ax1 = axes[0, 0]
    metrics = ['rms_deviation', 'max_deviation', 'mean_abs_deviation', 'std_deviation']
    orig_values = [orig_data['planar_deviation_metrics'][m] for m in metrics]
    proc_values = [proc_data['planar_deviation_metrics'][m] for m in metrics]
    
    x = np.arange(len(metrics))
    width = 0.35
    ax1.bar(x - width/2, orig_values, width, label='Original', alpha=0.8)
    ax1.bar(x + width/2, proc_values, width, label='Processed', alpha=0.8)
    ax1.set_xlabel('Deviation Metrics')
    ax1.set_ylabel('Value (km)')
    ax1.set_title('Planar Deviation Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels([m.replace('_', '\n') for m in metrics], fontsize=8)
    ax1.legend()
    
    # Plot 2: Roughness comparison
    ax2 = axes[0, 1]
    rough_metrics = ['rms_roughness', 'mean_local_roughness', 'max_local_roughness']
    orig_rough = [orig_data['surface_roughness_metrics'][m] for m in rough_metrics]
    proc_rough = [proc_data['surface_roughness_metrics'][m] for m in rough_metrics]
    
    x = np.arange(len(rough_metrics))
    ax2.bar(x - width/2, orig_rough, width, label='Original', alpha=0.8)
    ax2.bar(x + width/2, proc_rough, width, label='Processed', alpha=0.8)
    ax2.set_xlabel('Roughness Metrics')
    ax2.set_ylabel('Value (km)')
    ax2.set_title('Surface Roughness Comparison')
    ax2.set_xticks(x)
    ax2.set_xticklabels([m.replace('_', '\n') for m in rough_metrics], fontsize=8)
    ax2.legend()
    
    # Plot 3: Preservation ratios
    ax3 = axes[1, 0]
    ratio_names = list(ratios.keys())
    ratio_values = list(ratios.values())
    colors = ['green' if r < 1.2 else 'orange' if r < 1.5 else 'red' for r in ratio_values]
    
    bars = ax3.bar(range(len(ratio_names)), ratio_values, color=colors, alpha=0.8)
    ax3.set_xlabel('Preservation Metrics')
    ax3.set_ylabel('Ratio (Processed/Original)')
    ax3.set_title('Smoothness Preservation Ratios')
    ax3.set_xticks(range(len(ratio_names)))
    ax3.set_xticklabels([n.replace('_', '\n') for n in ratio_names], fontsize=8)
    ax3.axhline(y=1.0, color='black', linestyle='--', alpha=0.7, label='Perfect preservation')
    ax3.legend()
    
    # Add value labels on bars
    for bar, value in zip(bars, ratio_values):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{value:.2f}', ha='center', va='bottom', fontsize=9)
    
    # Plot 4: Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary = comparison_data['summary']
    summary_text = f"""Smoothness Preservation Summary

Original Smoothness Score: {summary['original_score']:.4f}
Processed Smoothness Score: {summary['processed_score']:.4f}
Score Change: {summary['score_change']:.4f}
Relative Change: {summary['relative_change_percent']:.1f}%

Average Preservation Ratio: {comparison_data['average_preservation_ratio']:.3f}
Preservation Quality: {comparison_data['preservation_quality'].upper()}

Interpretation:
• Ratio < 1.2: Excellent preservation
• Ratio 1.2-1.5: Good preservation  
• Ratio 1.5-2.0: Moderate preservation
• Ratio > 2.0: Poor preservation

Lower smoothness scores indicate smoother surfaces."""
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    plot_filename = f"{output_prefix}_smoothness_comparison.png"
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Smoothness comparison plot saved to {plot_filename}")
    plt.close()

def process_stl_to_grid(stl_filename, 
                       resolution=1.0, 
                       x_range=(-380, 380), 
                       z_range=(0, -50),
                       max_vertices=50000,
                       smoothing_method='none',
                       smoothing_params=None,
                       enable_vertical_transition=False,
                       transition_start_depth=None,
                       transition_distance=5.0,
                       plane_method='best_fit_line',
                       blend_function='cosine'):
    """
    Process STL file and project it onto a regular grid
    
    Parameters:
    -----------
    stl_filename : str
        Path to the STL file
    resolution : float
        Grid resolution in km (default: 1.0)
    x_range : tuple
        (x_min, x_max) in km (default: (-350, 350))
    z_range : tuple  
        (z_max, z_min) in km (default: (0, -50))
    max_vertices : int
        Maximum vertices to read from STL (default: 50000)
    smoothing_method : str
        Smoothing method to apply: 'gaussian', 'adaptive', 'rbf', 'constrained', 'none'
    smoothing_params : dict, optional
        Parameters for the smoothing method
    enable_vertical_transition : bool
        Whether to apply smooth vertical transition (default: False)
    transition_start_depth : float, optional
        Depth where transition begins (km). If None, uses STL bottom - 5km
    transition_distance : float
        Distance over which transition occurs (default: 5.0 km)
    plane_method : str
        Method for target vertical plane: 'best_fit_line', 'centroid_line', 'strike_aligned'
    blend_function : str
        Blending function: 'linear', 'cosine', 'smooth_step', 'exponential'
        
    Returns:
    --------
    grid_data : pandas.DataFrame
        Columns: ['x', 'z', 'y'] where y is the fault position at y=0 plane
        Ordered by z (bottom first), then x (-x to +x)
    metadata : dict
        Processing metadata including original strike angle
    """
    
    print(f"Processing {stl_filename} to regular grid")
    print(f"Resolution: {resolution} km")
    print(f"X range: {x_range[0]} to {x_range[1]} km")
    print(f"Z range: {z_range[0]} to {z_range[1]} km")
    
    # Read vertices from STL
    vertices = read_ascii_stl_vertices(stl_filename, max_vertices)
    
    # Convert to km and get basic info
    vertices_km = vertices / 1000.0
    
    print(f"Original coordinate ranges:")
    print(f"X: {np.min(vertices_km[:, 0]):.1f} to {np.max(vertices_km[:, 0]):.1f} km")
    print(f"Y: {np.min(vertices_km[:, 1]):.1f} to {np.max(vertices_km[:, 1]):.1f} km")
    print(f"Z: {np.min(vertices_km[:, 2]):.1f} to {np.max(vertices_km[:, 2]):.1f} km")
    
    # Calculate smoothness metrics for original STL
    print("\n=== ORIGINAL STL SMOOTHNESS ANALYSIS ===")
    original_smoothness = calculate_smoothness_metrics(vertices_km, resolution=resolution)
    
    # Calculate original strike angle
    original_strike = calculate_strike_angle_from_coords(vertices_km)
    print(f"Original average strike angle: {original_strike:.2f}°")
    
    # Rotate to align strike with x-axis
    rotation_angle = -original_strike
    angle_rad = np.radians(rotation_angle)
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad), 0],
        [np.sin(angle_rad), np.cos(angle_rad), 0],
        [0, 0, 1]
    ])
    vertices_rotated = np.dot(vertices_km, rotation_matrix.T)
    print(f"Rotated mesh by {rotation_angle:.2f}° to align with x-axis")
    
    # Recenter mesh so top center is at (0,0,0)
    max_z = np.max(vertices_rotated[:, 2])
    top_vertices = vertices_rotated[vertices_rotated[:, 2] > max_z - 1.0]
    
    if len(top_vertices) > 0:
        center_x = np.mean(top_vertices[:, 0])
        center_y = np.mean(top_vertices[:, 1])
    else:
        center_x = np.mean(vertices_rotated[:, 0])
        center_y = np.mean(vertices_rotated[:, 1])
    
    vertices_centered = vertices_rotated.copy()
    vertices_centered[:, 0] -= center_x
    vertices_centered[:, 1] -= center_y
    vertices_centered[:, 2] -= max_z
    print("Recentered mesh so top center is at (0,0,0)")
    
    print(f"Processed coordinate ranges:")
    print(f"X: {np.min(vertices_centered[:, 0]):.1f} to {np.max(vertices_centered[:, 0]):.1f} km")
    print(f"Y: {np.min(vertices_centered[:, 1]):.1f} to {np.max(vertices_centered[:, 1]):.1f} km")
    print(f"Z: {np.min(vertices_centered[:, 2]):.1f} to {np.max(vertices_centered[:, 2]):.1f} km")
    
    # Determine the depth boundary for STL data vs vertical extension
    stl_min_z = np.min(vertices_centered[:, 2])
    extension_boundary = min(-20, stl_min_z)  # Use STL min depth or -20km, whichever is deeper
    
    print(f"STL data extends to {stl_min_z:.1f} km depth")
    print(f"Will extend vertically below {extension_boundary:.1f} km")
    
    # Create regular grid
    x_grid = np.arange(x_range[0], x_range[1] + resolution, resolution)
    z_grid = np.arange(z_range[0], z_range[1] - resolution, -resolution)  # Top to bottom
    
    print(f"Grid dimensions: {len(x_grid)} x {len(z_grid)} = {len(x_grid) * len(z_grid)} points")
    
    # Split grid into two parts: STL interpolation zone and vertical extension zone
    stl_zone_mask = z_grid >= extension_boundary
    extension_zone_mask = z_grid < extension_boundary
    
    stl_z_grid = z_grid[stl_zone_mask]
    extension_z_grid = z_grid[extension_zone_mask]
    
    print(f"STL interpolation zone: {len(stl_z_grid)} z-levels ({np.max(stl_z_grid):.1f} to {np.min(stl_z_grid):.1f} km)")
    print(f"Vertical extension zone: {len(extension_z_grid)} z-levels ({np.max(extension_z_grid):.1f} to {np.min(extension_z_grid):.1f} km)")
    
    # Process STL interpolation zone
    if len(stl_z_grid) > 0:
        print("Interpolating STL data onto regular grid...")
        
        X_stl, Z_stl = np.meshgrid(x_grid, stl_z_grid)
        stl_grid_points = np.column_stack([X_stl.ravel(), Z_stl.ravel()])
        
        # Use only vertices that are within reasonable bounds for interpolation
        valid_mask = ((vertices_centered[:, 0] >= x_range[0] - 50) & 
                      (vertices_centered[:, 0] <= x_range[1] + 50) &
                      (vertices_centered[:, 2] >= np.min(stl_z_grid) - 5) & 
                      (vertices_centered[:, 2] <= np.max(stl_z_grid) + 5))
        
        valid_vertices = vertices_centered[valid_mask]
        print(f"Using {len(valid_vertices)} vertices for STL interpolation")
        
        # Perform interpolation for STL zone
        y_stl_interpolated = griddata(
            points=valid_vertices[:, [0, 2]],  # x, z coordinates
            values=valid_vertices[:, 1],       # y values
            xi=stl_grid_points,                # regular grid points
            method='linear',                   # linear interpolation
            fill_value=np.nan                  # fill with NaN outside convex hull
        )
        
        # Handle NaN values by using nearest neighbor interpolation
        nan_mask = np.isnan(y_stl_interpolated)
        if np.any(nan_mask):
            print(f"Filling {np.sum(nan_mask)} NaN values with nearest neighbor interpolation")
            y_nearest = griddata(
                points=valid_vertices[:, [0, 2]],
                values=valid_vertices[:, 1],
                xi=stl_grid_points[nan_mask],
                method='nearest'
            )
            y_stl_interpolated[nan_mask] = y_nearest
    else:
        y_stl_interpolated = np.array([])
        stl_grid_points = np.empty((0, 2))
    
    # Process vertical extension zone
    if len(extension_z_grid) > 0:
        print("Creating vertical extension below STL data...")
        
        # Get y-values at the extension boundary (bottom of STL zone)
        boundary_y_values = {}
        
        if len(stl_z_grid) > 0:
            # Use the bottommost interpolated values from STL zone
            boundary_z = np.min(stl_z_grid)
            boundary_idx_start = (len(stl_z_grid) - 1) * len(x_grid)  # Last z-level in STL zone
            boundary_idx_end = boundary_idx_start + len(x_grid)
            
            for i, x in enumerate(x_grid):
                boundary_y_values[x] = y_stl_interpolated[boundary_idx_start + i]
        else:
            # Fallback: interpolate at extension boundary from original vertices
            print("No STL interpolation zone, using direct interpolation at boundary")
            boundary_points = np.column_stack([x_grid, np.full(len(x_grid), extension_boundary)])
            boundary_y_interp = griddata(
                points=valid_vertices[:, [0, 2]],
                values=valid_vertices[:, 1],
                xi=boundary_points,
                method='linear',
                fill_value=np.nan
            )
            # Fill NaN with nearest neighbor
            nan_mask = np.isnan(boundary_y_interp)
            if np.any(nan_mask):
                boundary_y_nearest = griddata(
                    points=valid_vertices[:, [0, 2]],
                    values=valid_vertices[:, 1],
                    xi=boundary_points[nan_mask],
                    method='nearest'
                )
                boundary_y_interp[nan_mask] = boundary_y_nearest
            
            for i, x in enumerate(x_grid):
                boundary_y_values[x] = boundary_y_interp[i]
        
        # Create vertical extension with smooth transition to optimal vertical plane
        # Find the optimal y-value for vertical plane that minimizes transition discontinuity
        X_ext, Z_ext = np.meshgrid(x_grid, extension_z_grid)
        extension_grid_points = np.column_stack([X_ext.ravel(), Z_ext.ravel()])
        
        y_extension = np.zeros(len(extension_grid_points))
        
        # Find optimal constant y-value or linear trend for vertical plane
        # Use boundary y-values to determine best-fit vertical plane
        boundary_x_values = list(boundary_y_values.keys())
        boundary_y_vals = list(boundary_y_values.values())
        
        if len(boundary_x_values) > 2:
            # Fit linear trend y = mx + b to boundary values
            coeffs = np.polyfit(boundary_x_values, boundary_y_vals, 1)
            slope, intercept = coeffs
            
            # Decide whether to use linear trend or mean value based on fit quality
            y_predicted = np.polyval(coeffs, boundary_x_values)
            r_squared = 1 - np.sum((boundary_y_vals - y_predicted)**2) / np.sum((boundary_y_vals - np.mean(boundary_y_vals))**2)
            
            if r_squared > 0.7:  # Good linear fit, use linear trend
                print(f"Using linear vertical plane: y = {slope:.6f}*x + {intercept:.3f} (R²={r_squared:.3f})")
                use_linear = True
                # Store plane parameters for lateral extension
                vertical_plane_slope = slope
                vertical_plane_intercept = intercept  
                vertical_plane_is_linear = True
            else:  # Poor linear fit, use mean value
                mean_y = np.mean(boundary_y_vals)
                print(f"Using constant vertical plane: y = {mean_y:.3f} (poor linear fit, R²={r_squared:.3f})")
                use_linear = False
                # Store plane parameters for lateral extension
                vertical_plane_slope = 0.0
                vertical_plane_intercept = mean_y
                vertical_plane_is_linear = False
        else:
            # Not enough points for linear fit, use mean
            mean_y = np.mean(boundary_y_vals) if boundary_y_vals else 0.0
            print(f"Using constant vertical plane: y = {mean_y:.3f} (insufficient points for linear fit)")
            use_linear = False
            # Store plane parameters for lateral extension
            vertical_plane_slope = 0.0
            vertical_plane_intercept = mean_y
            vertical_plane_is_linear = False
        
        # Calculate transition distance for smooth decay to vertical plane
        total_extension_depth = abs(extension_boundary - np.min(extension_z_grid))
        transition_depth = min(10.0, total_extension_depth * 0.5)  # Use up to 10km or half the extension
        transition_start = extension_boundary
        transition_end = extension_boundary - transition_depth
        
        for i, (x, z) in enumerate(extension_grid_points):
            # Calculate target y-value for vertical plane
            if use_linear:
                target_y = slope * x + intercept
            else:
                target_y = mean_y
            
            if x in boundary_y_values:
                boundary_y = boundary_y_values[x]
                
                if z >= transition_start:
                    # At boundary, use boundary value
                    y_extension[i] = boundary_y
                elif z <= transition_end:
                    # Below transition, use vertical plane value
                    y_extension[i] = target_y
                else:
                    # In transition zone, smoothly blend from boundary_y to target_y
                    progress = (transition_start - z) / transition_depth
                    # Use cosine transition for smoothness
                    blend_factor = 0.5 * (1 - np.cos(np.pi * progress))
                    y_extension[i] = boundary_y * (1 - blend_factor) + target_y * blend_factor
            else:
                # If x not in boundary_y_values, use vertical plane value
                y_extension[i] = target_y
        
        print(f"Created {len(y_extension)} points in vertical extension zone")
    else:
        y_extension = np.array([])
        extension_grid_points = np.empty((0, 2))
    
    # Initialize combined arrays first
    all_grid_points = np.vstack([stl_grid_points, extension_grid_points]) if len(extension_grid_points) > 0 else stl_grid_points
    all_y_values = np.concatenate([y_stl_interpolated, y_extension]) if len(y_extension) > 0 else y_stl_interpolated
    
    # Create comprehensive grid with proper zone-based extensions
    print("\\n=== CREATING COMPREHENSIVE GRID WITH ZONE-BASED EXTENSIONS ===")
    
    # Define zones clearly
    stl_x_min, stl_x_max = -350, 350  # STL interpolation zone
    stl_z_min, stl_z_max = 0, -20     # STL depth zone (corrected: z_max > z_min)  
    
    extension_x_inner = 365   # Inner boundary for pure lateral extension
    extension_z_inner = -35   # Inner boundary for pure vertical extension
    
    transition_width_x = 15   # X transition: 350 to 365 and -365 to -350  
    transition_width_z = 15   # Z transition: -20 to -35
    
    # Use the same vertical plane parameters for all extensions
    if 'vertical_plane_is_linear' in locals() and vertical_plane_is_linear:
        plane_func = lambda x: vertical_plane_slope * x + vertical_plane_intercept
        print(f"Using extension plane: y = {vertical_plane_slope:.6f}*x + {vertical_plane_intercept:.3f}")
    else:
        plane_func = lambda x: vertical_plane_intercept
        print(f"Using extension plane: y = {vertical_plane_intercept:.3f}")
    
    # Create comprehensive grid for all zones
    all_grid_points_comprehensive = []
    all_y_values_comprehensive = []
    
    # Process every point in the extended grid
    for i, z in enumerate(z_grid):
        for j, x in enumerate(x_grid):
            point_added = False
            
            # Zone 1: STL interpolation zone (x ∈ [-350,350], z > -20)
            if stl_x_min <= x <= stl_x_max and z > stl_z_max:
                # Use existing STL interpolation
                idx = i * len(x_grid) + j
                if idx < len(all_y_values):
                    all_grid_points_comprehensive.append([x, z])
                    all_y_values_comprehensive.append(all_y_values[idx])
                    point_added = True
            
            # Zone 2: Pure extension zone (x < -365 OR x > 365, AND z < -35)
            elif (x < -extension_x_inner or x > extension_x_inner) and z < extension_z_inner:
                all_grid_points_comprehensive.append([x, z])
                all_y_values_comprehensive.append(plane_func(x))
                point_added = True
            
            # Zone 3: All other areas - interpolate from STL boundary at z=-20 to plane at z=-35
            else:
                # Get the STL y-value at z=-20 for the current x position
                if len(y_stl_interpolated) > 0 and len(stl_z_grid) > 0:
                    bottom_z_idx = len(stl_z_grid) - 1  # Bottom of STL zone (z=-20)
                    
                    if stl_x_min <= x <= stl_x_max:
                        # Within STL x-range: get the actual STL y-value at this x position
                        stl_boundary_idx = bottom_z_idx * len(x_grid) + j
                        if stl_boundary_idx < len(y_stl_interpolated):
                            stl_y_at_boundary = y_stl_interpolated[stl_boundary_idx]
                        else:
                            stl_y_at_boundary = plane_func(x)
                    else:
                        # Outside STL x-range: use plane value
                        stl_y_at_boundary = plane_func(x)
                else:
                    stl_y_at_boundary = plane_func(x)
                
                # Get the plane y-value at z=-35
                plane_y_at_depth = plane_func(x)
                
                # Linear interpolation from STL boundary (z=-20) to plane (z=-35)
                z_progress = (stl_z_max - z) / (stl_z_max - extension_z_inner)  # 0 at z=-20, 1 at z=-35
                z_progress = np.clip(z_progress, 0, 1)
                
                y_value = stl_y_at_boundary * (1 - z_progress) + plane_y_at_depth * z_progress
                
                if y_value is not None:
                    all_grid_points_comprehensive.append([x, z])
                    all_y_values_comprehensive.append(y_value)
                    point_added = True
                else:
                    # Fallback: if no rule matched, use the plane function
                    all_grid_points_comprehensive.append([x, z])
                    all_y_values_comprehensive.append(plane_func(x))
                    point_added = True
    
    print(f"Created comprehensive grid with {len(all_y_values_comprehensive)} points")
    print(f"Grid zones: STL + transitions + extensions")
    
    # Update the grid arrays
    all_grid_points = np.array(all_grid_points_comprehensive)
    all_y_values = np.array(all_y_values_comprehensive)
    
    # Update x_grid to include the full range
    x_grid = np.arange(x_range[0], x_range[1] + resolution, resolution)
    
    # Create 2D grid structure for proper derivative calculation
    print("Creating 2D grid structure...")
    
    # Create the 2D grid: Y_grid[z_idx, x_idx]
    Y_grid = np.full((len(z_grid), len(x_grid)), np.nan)
    
    # Fill the 2D grid with interpolated/extended y values
    for i, z in enumerate(z_grid):  # z_grid goes from 0 to -50 (top to bottom)
        for j, x in enumerate(x_grid):  # x_grid goes from -350 to +350
            idx = i * len(x_grid) + j
            if idx < len(all_y_values):
                Y_grid[i, j] = all_y_values[idx]
    
    print(f"2D grid shape: {Y_grid.shape} ({len(z_grid)} z-levels × {len(x_grid)} x-points)")
    print(f"Grid coverage: {np.sum(~np.isnan(Y_grid))} / {Y_grid.size} points ({100*np.sum(~np.isnan(Y_grid))/Y_grid.size:.1f}%)")
    
    # Calculate derivatives using central differences on the 2D grid
    print("Calculating derivatives dy/dx and dy/dz on 2D grid...")
    
    # Initialize derivative grids
    dydx_grid = np.full_like(Y_grid, np.nan)
    dydz_grid = np.full_like(Y_grid, np.nan)
    
    # dy/dx: derivative along x direction (axis=1, columns)
    # Central difference: (y[i,j+1] - y[i,j-1]) / (2*dx)
    for i in range(len(z_grid)):
        for j in range(len(x_grid)):
            if not np.isnan(Y_grid[i, j]):
                if 1 <= j <= len(x_grid)-2:  # Central difference
                    if not np.isnan(Y_grid[i, j-1]) and not np.isnan(Y_grid[i, j+1]):
                        dydx_grid[i, j] = (Y_grid[i, j+1] - Y_grid[i, j-1]) / (2 * resolution)
                elif j == 0:  # Forward difference at left boundary
                    if not np.isnan(Y_grid[i, j+1]):
                        dydx_grid[i, j] = (Y_grid[i, j+1] - Y_grid[i, j]) / resolution
                elif j == len(x_grid)-1:  # Backward difference at right boundary
                    if not np.isnan(Y_grid[i, j-1]):
                        dydx_grid[i, j] = (Y_grid[i, j] - Y_grid[i, j-1]) / resolution
    
    # dy/dz: derivative along z direction (axis=0, rows)
    # Central difference: (y[i+1,j] - y[i-1,j]) / (2*dz)
    # Note: z decreases with increasing i, so dz = -resolution
    for i in range(len(z_grid)):
        for j in range(len(x_grid)):
            if not np.isnan(Y_grid[i, j]):
                if 1 <= i <= len(z_grid)-2:  # Central difference
                    if not np.isnan(Y_grid[i-1, j]) and not np.isnan(Y_grid[i+1, j]):
                        # Since z decreases with i: dz = z[i+1] - z[i-1] = -2*resolution
                        dydz_grid[i, j] = (Y_grid[i+1, j] - Y_grid[i-1, j]) / (-2 * resolution)
                elif i == 0:  # Forward difference at top boundary
                    if not np.isnan(Y_grid[i+1, j]):
                        # dz = z[i+1] - z[i] = -resolution
                        dydz_grid[i, j] = (Y_grid[i+1, j] - Y_grid[i, j]) / (-resolution)
                elif i == len(z_grid)-1:  # Backward difference at bottom boundary
                    if not np.isnan(Y_grid[i-1, j]):
                        # dz = z[i] - z[i-1] = -resolution
                        dydz_grid[i, j] = (Y_grid[i, j] - Y_grid[i-1, j]) / (-resolution)
    
    # Create final data structure: order bottom z first, then -x to +x
    print("Creating final ordered data structure...")
    grid_data = []
    
    # Process all z-levels from bottom to top (reverse order)
    for i in reversed(range(len(z_grid))):  # Bottom to top
        z = z_grid[i]
        for j in range(len(x_grid)):  # -x to +x
            x = x_grid[j]
            
            y_val = Y_grid[i, j]
            dydx_val = dydx_grid[i, j]
            dydz_val = dydz_grid[i, j]
            
            # Only include points where y is defined
            if not np.isnan(y_val):
                # Use 0.0 for undefined derivatives instead of NaN
                if np.isnan(dydx_val):
                    dydx_val = 0.0
                if np.isnan(dydz_val):
                    dydz_val = 0.0
                    
                grid_data.append([x, z, y_val, dydx_val, dydz_val])
    
    # Convert to DataFrame
    df = pd.DataFrame(grid_data, columns=['x', 'z', 'y', 'dy_dx', 'dy_dz'])
    
    print(f"Final grid data: {len(df)} points")
    print(f"Grid y range: {df['y'].min():.1f} to {df['y'].max():.1f} km")
    print(f"dy/dx range: {df['dy_dx'].min():.3f} to {df['dy_dx'].max():.3f}")
    print(f"dy/dz range: {df['dy_dz'].min():.3f} to {df['dy_dz'].max():.3f}")
    
    # Apply smoothing control if requested
    smoothing_info = None
    if smoothing_method != 'none':
        print(f"\n=== APPLYING SMOOTHING CONTROL ({smoothing_method}) ===")
        processed_vertices_raw = np.column_stack([df['x'], df['y'], df['z']])
        processed_grid_shape = (len(z_grid), len(x_grid))
        
        # Apply smoothing
        processed_vertices_smoothed, smoothing_info = apply_smoothing_control(
            processed_vertices_raw, smoothing_method, smoothing_params, processed_grid_shape)
        
        # Update grid data with smoothed vertices
        df['y'] = processed_vertices_smoothed[:, 1]
        
        # Recalculate derivatives for smoothed surface
        print("Recalculating derivatives for smoothed surface...")
        Y_grid_smoothed = processed_vertices_smoothed[:, 1].reshape(len(z_grid), len(x_grid))
        
        # Recalculate derivatives
        dydx_grid = np.full_like(Y_grid_smoothed, np.nan)
        dydz_grid = np.full_like(Y_grid_smoothed, np.nan)
        
        # dy/dx: derivative along x direction
        for i in range(len(z_grid)):
            for j in range(len(x_grid)):
                if not np.isnan(Y_grid_smoothed[i, j]):
                    if 1 <= j <= len(x_grid)-2:  # Central difference
                        if not np.isnan(Y_grid_smoothed[i, j-1]) and not np.isnan(Y_grid_smoothed[i, j+1]):
                            dydx_grid[i, j] = (Y_grid_smoothed[i, j+1] - Y_grid_smoothed[i, j-1]) / (2 * resolution)
                    elif j == 0:  # Forward difference
                        if not np.isnan(Y_grid_smoothed[i, j+1]):
                            dydx_grid[i, j] = (Y_grid_smoothed[i, j+1] - Y_grid_smoothed[i, j]) / resolution
                    elif j == len(x_grid)-1:  # Backward difference
                        if not np.isnan(Y_grid_smoothed[i, j-1]):
                            dydx_grid[i, j] = (Y_grid_smoothed[i, j] - Y_grid_smoothed[i, j-1]) / resolution
        
        # dy/dz: derivative along z direction
        for i in range(len(z_grid)):
            for j in range(len(x_grid)):
                if not np.isnan(Y_grid_smoothed[i, j]):
                    if 1 <= i <= len(z_grid)-2:  # Central difference
                        if not np.isnan(Y_grid_smoothed[i-1, j]) and not np.isnan(Y_grid_smoothed[i+1, j]):
                            dydz_grid[i, j] = (Y_grid_smoothed[i+1, j] - Y_grid_smoothed[i-1, j]) / (-2 * resolution)
                    elif i == 0:  # Forward difference
                        if not np.isnan(Y_grid_smoothed[i+1, j]):
                            dydz_grid[i, j] = (Y_grid_smoothed[i+1, j] - Y_grid_smoothed[i, j]) / (-resolution)
                    elif i == len(z_grid)-1:  # Backward difference
                        if not np.isnan(Y_grid_smoothed[i-1, j]):
                            dydz_grid[i, j] = (Y_grid_smoothed[i, j] - Y_grid_smoothed[i-1, j]) / (-resolution)
        
        # Update derivatives in DataFrame
        point_idx = 0
        for i in reversed(range(len(z_grid))):  # Bottom to top
            for j in range(len(x_grid)):  # -x to +x
                if point_idx < len(df):
                    dydx_val = dydx_grid[i, j]
                    dydz_val = dydz_grid[i, j]
                    
                    if np.isnan(dydx_val):
                        dydx_val = 0.0
                    if np.isnan(dydz_val):
                        dydz_val = 0.0
                    
                    df.iloc[point_idx, df.columns.get_loc('dy_dx')] = dydx_val
                    df.iloc[point_idx, df.columns.get_loc('dy_dz')] = dydz_val
                    point_idx += 1
        
        print(f"Applied {smoothing_method} smoothing successfully")
        processed_vertices = processed_vertices_smoothed
    else:
        processed_vertices = np.column_stack([df['x'], df['y'], df['z']])
        processed_grid_shape = (len(z_grid), len(x_grid))
    
    # Apply vertical transition if requested
    transition_info = None
    if enable_vertical_transition:
        print(f"\n=== APPLYING VERTICAL TRANSITION ===")
        
        # Determine transition start depth if not specified
        if transition_start_depth is None:
            stl_bottom = np.min(vertices_centered[:, 2])
            auto_transition_start = stl_bottom - 5.0  # 5km below STL bottom
            print(f"Auto-determined transition start depth: {auto_transition_start:.1f} km (STL bottom: {stl_bottom:.1f} km)")
        else:
            auto_transition_start = transition_start_depth
            print(f"Using specified transition start depth: {auto_transition_start:.1f} km")
        
        # Apply vertical transition
        processed_vertices_pre_transition = processed_vertices.copy()
        processed_vertices_transitioned, transition_info = apply_smooth_vertical_transition(
            processed_vertices, processed_grid_shape, auto_transition_start, 
            transition_distance, plane_method, blend_function)
        
        # Update grid data with transitioned vertices
        if 'error' not in transition_info:
            df['y'] = processed_vertices_transitioned[:, 1]
            
            # Recalculate derivatives for transitioned surface
            print("Recalculating derivatives for transitioned surface...")
            Y_grid_transitioned = processed_vertices_transitioned[:, 1].reshape(len(z_grid), len(x_grid))
            
            # Recalculate derivatives (similar to smoothing section)
            dydx_grid = np.full_like(Y_grid_transitioned, np.nan)
            dydz_grid = np.full_like(Y_grid_transitioned, np.nan)
            
            # dy/dx: derivative along x direction
            for i in range(len(z_grid)):
                for j in range(len(x_grid)):
                    if not np.isnan(Y_grid_transitioned[i, j]):
                        if 1 <= j <= len(x_grid)-2:  # Central difference
                            if not np.isnan(Y_grid_transitioned[i, j-1]) and not np.isnan(Y_grid_transitioned[i, j+1]):
                                dydx_grid[i, j] = (Y_grid_transitioned[i, j+1] - Y_grid_transitioned[i, j-1]) / (2 * resolution)
                        elif j == 0:  # Forward difference
                            if not np.isnan(Y_grid_transitioned[i, j+1]):
                                dydx_grid[i, j] = (Y_grid_transitioned[i, j+1] - Y_grid_transitioned[i, j]) / resolution
                        elif j == len(x_grid)-1:  # Backward difference
                            if not np.isnan(Y_grid_transitioned[i, j-1]):
                                dydx_grid[i, j] = (Y_grid_transitioned[i, j] - Y_grid_transitioned[i, j-1]) / resolution
            
            # dy/dz: derivative along z direction
            for i in range(len(z_grid)):
                for j in range(len(x_grid)):
                    if not np.isnan(Y_grid_transitioned[i, j]):
                        if 1 <= i <= len(z_grid)-2:  # Central difference
                            if not np.isnan(Y_grid_transitioned[i-1, j]) and not np.isnan(Y_grid_transitioned[i+1, j]):
                                dydz_grid[i, j] = (Y_grid_transitioned[i+1, j] - Y_grid_transitioned[i-1, j]) / (-2 * resolution)
                        elif i == 0:  # Forward difference
                            if not np.isnan(Y_grid_transitioned[i+1, j]):
                                dydz_grid[i, j] = (Y_grid_transitioned[i+1, j] - Y_grid_transitioned[i, j]) / (-resolution)
                        elif i == len(z_grid)-1:  # Backward difference
                            if not np.isnan(Y_grid_transitioned[i-1, j]):
                                dydz_grid[i, j] = (Y_grid_transitioned[i, j] - Y_grid_transitioned[i-1, j]) / (-resolution)
            
            # Update derivatives in DataFrame
            point_idx = 0
            for i in reversed(range(len(z_grid))):  # Bottom to top
                for j in range(len(x_grid)):  # -x to +x
                    if point_idx < len(df):
                        dydx_val = dydx_grid[i, j]
                        dydz_val = dydz_grid[i, j]
                        
                        if np.isnan(dydx_val):
                            dydx_val = 0.0
                        if np.isnan(dydz_val):
                            dydz_val = 0.0
                        
                        df.iloc[point_idx, df.columns.get_loc('dy_dx')] = dydx_val
                        df.iloc[point_idx, df.columns.get_loc('dy_dz')] = dydz_val
                        point_idx += 1
            
            processed_vertices = processed_vertices_transitioned
            print(f"Applied vertical transition successfully")
        else:
            print(f"Vertical transition failed: {transition_info['error']}")
            processed_vertices = processed_vertices_pre_transition
    
    # Calculate smoothness metrics for processed mesh
    print("\n=== PROCESSED MESH SMOOTHNESS ANALYSIS ===")
    processed_smoothness = calculate_smoothness_metrics(
        processed_vertices, grid_shape=processed_grid_shape, resolution=resolution)
    
    # Compare smoothness preservation
    print("\n=== SMOOTHNESS PRESERVATION ANALYSIS ===")
    smoothness_comparison = compare_smoothness_preservation(
        vertices_km, processed_vertices, 
        processed_grid_shape=processed_grid_shape, 
        resolution=resolution)
    
    # Create metadata
    metadata = {
        'stl_filename': stl_filename,
        'original_strike_angle': original_strike,
        'rotation_angle': rotation_angle,
        'resolution': resolution,
        'x_range': x_range,
        'z_range': z_range,
        'extension_boundary': extension_boundary,
        'grid_dimensions': (len(x_grid), len(z_grid)),
        'total_points': len(df),
        'y_range': (df['y'].min(), df['y'].max()),
        'processing_info': {
            'vertices_loaded': len(vertices),
            'vertices_used_for_interpolation': len(valid_vertices) if 'valid_vertices' in locals() else 0,
            'center_offsets': (center_x, center_y, max_z),
            'stl_zone_points': len(stl_z_grid) * len(x_grid) if len(stl_z_grid) > 0 else 0,
            'extension_zone_points': len(extension_z_grid) * len(x_grid) if len(extension_z_grid) > 0 else 0
        },
        # New smoothness analysis data
        'smoothness_analysis': {
            'original_smoothness': original_smoothness,
            'processed_smoothness': processed_smoothness,
            'smoothness_comparison': smoothness_comparison
        },
        # Smoothing control information
        'smoothing_control': {
            'method': smoothing_method,
            'parameters': smoothing_params,
            'smoothing_info': smoothing_info
        },
        # Vertical transition information
        'vertical_transition': {
            'enabled': enable_vertical_transition,
            'start_depth': transition_start_depth,
            'distance': transition_distance,
            'plane_method': plane_method,
            'blend_function': blend_function,
            'transition_info': transition_info
        }
    }
    
    # Return all data for plotting
    processing_data = {
        'vertices_original': vertices_km,
        'vertices_rotated': vertices_rotated, 
        'vertices_centered': vertices_centered,
        'grid_data': df,
        'metadata': metadata
    }
    
    return processing_data

def plot_processing_steps(processing_data, output_prefix):
    """
    Plot original nodes, rotated/centered nodes, and final mesh
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    vertices_original = processing_data['vertices_original']
    vertices_rotated = processing_data['vertices_rotated'] 
    vertices_centered = processing_data['vertices_centered']
    grid_data = processing_data['grid_data']
    metadata = processing_data['metadata']
    
    fig = plt.figure(figsize=(15, 5))
    
    # Plot 1: Original nodes
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.scatter(vertices_original[:, 0], vertices_original[:, 1], vertices_original[:, 2], 
                c=vertices_original[:, 2], cmap='viridis', s=0.5, alpha=0.6)
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_zlabel('Z (km)')
    ax1.set_title(f'Original STL Nodes\n{len(vertices_original)} points')
    ax1.set_box_aspect([1,1,1])
    
    # Plot 2: Rotated and centered nodes
    ax2 = fig.add_subplot(132, projection='3d')
    ax2.scatter(vertices_centered[:, 0], vertices_centered[:, 1], vertices_centered[:, 2], 
                c=vertices_centered[:, 2], cmap='viridis', s=0.5, alpha=0.6)
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Y (km)')
    ax2.set_zlabel('Z (km)')
    ax2.set_title(f'Rotated & Centered\nStrike: {metadata["original_strike_angle"]:.1f}° → 0°')
    ax2.set_box_aspect([1,1,1])
    
    # Plot 3: Final regular grid
    ax3 = fig.add_subplot(133, projection='3d')
    ax3.scatter(grid_data['x'], grid_data['y'], grid_data['z'], 
                c=grid_data['z'], cmap='viridis', s=0.5, alpha=0.6)
    ax3.set_xlabel('X (km)')
    ax3.set_ylabel('Y (km)')
    ax3.set_zlabel('Z (km)')
    ax3.set_title(f'Regular Grid\n{len(grid_data)} points, {metadata["resolution"]}km res')
    ax3.set_box_aspect([1,1,1])
    
    plt.tight_layout()
    
    plot_filename = f"{output_prefix}_processing.png"
    plt.savefig(plot_filename, dpi=200, bbox_inches='tight')
    print(f"Processing plots saved to {plot_filename}")
    plt.close()

def save_geometry_txt(grid_data, metadata, output_filename):
    """
    Save grid data to Geometry.txt format with headers
    Format:
    Line 1: nx nz 0
    Line 2: resolution min_x min_z
    Line 3+: y(m) dy/dx dy/dz
    """
    with open(output_filename, 'w') as f:
        # Extract grid parameters
        nx, nz = metadata['grid_dimensions']
        resolution = metadata['resolution']  # in km
        x_min = metadata['x_range'][0]       # in km
        z_min = metadata['z_range'][1]       # in km (most negative value)
        
        # Line 1: number of nodes along x, number of nodes along z, and a zero
        f.write(f"{nx} {nz} 0\n")
        
        # Line 2: grid resolution (m), minimum x (m), minimum z (m)
        f.write(f"{resolution*1000} {x_min*1000} {z_min*1000}\n")
        
        # Line 3+: Write data, convert km to meters
        for _, row in grid_data.iterrows():
            y_m = row['y'] * 1000  # Convert km to m
            dydx = row['dy_dx']    # Already dimensionless
            dydz = row['dy_dz']    # Already dimensionless
            
            f.write(f"{y_m:10.3f} {dydx:12.6f} {dydz:12.6f}\n")
    
    print(f"Geometry data saved to {output_filename} (header + y(m) dy/dx dy/dz format)")

def save_smoothness_report(processing_data, output_filename):
    """
    Save detailed smoothness analysis report
    """
    metadata = processing_data['metadata']
    smoothness_data = metadata['smoothness_analysis']
    
    with open(output_filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("FAULT SURFACE SMOOTHNESS ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"STL File: {metadata['stl_filename']}\n")
        f.write(f"Processing Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Grid Resolution: {metadata['resolution']} km\n")
        f.write(f"Grid Dimensions: {metadata['grid_dimensions'][0]} x {metadata['grid_dimensions'][1]}\n")
        f.write(f"Total Points: {metadata['total_points']}\n\n")
        
        # Original STL smoothness
        orig_smooth = smoothness_data['original_smoothness']
        f.write("ORIGINAL STL SMOOTHNESS METRICS:\n")
        f.write("-" * 40 + "\n")
        
        planar_orig = orig_smooth['planar_deviation_metrics']
        f.write(f"Planar Deviation Metrics:\n")
        f.write(f"  RMS Deviation: {planar_orig['rms_deviation']:.6f} km\n")
        f.write(f"  Max Deviation: {planar_orig['max_deviation']:.6f} km\n")
        f.write(f"  Mean Abs Deviation: {planar_orig['mean_abs_deviation']:.6f} km\n")
        f.write(f"  Standard Deviation: {planar_orig['std_deviation']:.6f} km\n")
        f.write(f"  Deviation Range: {planar_orig['deviation_range']:.6f} km\n\n")
        
        rough_orig = orig_smooth['surface_roughness_metrics']
        f.write(f"Surface Roughness Metrics:\n")
        f.write(f"  RMS Roughness: {rough_orig['rms_roughness']:.6f} km\n")
        f.write(f"  Mean Local Roughness: {rough_orig['mean_local_roughness']:.6f} km\n")
        f.write(f"  Max Local Roughness: {rough_orig['max_local_roughness']:.6f} km\n")
        f.write(f"  Neighborhood Size: {rough_orig['neighborhood_size']}\n\n")
        
        f.write(f"Overall Smoothness Score: {orig_smooth['overall_smoothness_score']:.6f}\n\n")
        
        # Processed mesh smoothness
        proc_smooth = smoothness_data['processed_smoothness']
        f.write("PROCESSED MESH SMOOTHNESS METRICS:\n")
        f.write("-" * 40 + "\n")
        
        planar_proc = proc_smooth['planar_deviation_metrics']
        f.write(f"Planar Deviation Metrics:\n")
        f.write(f"  RMS Deviation: {planar_proc['rms_deviation']:.6f} km\n")
        f.write(f"  Max Deviation: {planar_proc['max_deviation']:.6f} km\n")
        f.write(f"  Mean Abs Deviation: {planar_proc['mean_abs_deviation']:.6f} km\n")
        f.write(f"  Standard Deviation: {planar_proc['std_deviation']:.6f} km\n")
        f.write(f"  Deviation Range: {planar_proc['deviation_range']:.6f} km\n\n")
        
        rough_proc = proc_smooth['surface_roughness_metrics']
        f.write(f"Surface Roughness Metrics:\n")
        f.write(f"  RMS Roughness: {rough_proc['rms_roughness']:.6f} km\n")
        f.write(f"  Mean Local Roughness: {rough_proc['mean_local_roughness']:.6f} km\n")
        f.write(f"  Max Local Roughness: {rough_proc['max_local_roughness']:.6f} km\n")
        f.write(f"  Neighborhood Size: {rough_proc['neighborhood_size']}\n\n")
        
        f.write(f"Overall Smoothness Score: {proc_smooth['overall_smoothness_score']:.6f}\n\n")
        
        # Smoothness preservation analysis
        comparison = smoothness_data['smoothness_comparison']
        f.write("SMOOTHNESS PRESERVATION ANALYSIS:\n")
        f.write("-" * 40 + "\n")
        
        ratios = comparison['preservation_ratios']
        f.write(f"Preservation Ratios (Processed/Original):\n")
        f.write(f"  RMS Deviation Ratio: {ratios['rms_deviation_ratio']:.3f}\n")
        f.write(f"  Max Deviation Ratio: {ratios['max_deviation_ratio']:.3f}\n")
        f.write(f"  Roughness Ratio: {ratios['roughness_ratio']:.3f}\n")
        f.write(f"  Smoothness Score Ratio: {ratios['smoothness_score_ratio']:.3f}\n\n")
        
        f.write(f"Average Preservation Ratio: {comparison['average_preservation_ratio']:.3f}\n")
        f.write(f"Preservation Quality: {comparison['preservation_quality'].upper()}\n\n")
        
        summary = comparison['summary']
        f.write(f"Summary:\n")
        f.write(f"  Score Change: {summary['score_change']:.6f}\n")
        f.write(f"  Relative Change: {summary['relative_change_percent']:.2f}%\n\n")
        
        # Reference plane information
        plane_orig = orig_smooth['planar_reference']
        f.write("REFERENCE PLANE INFORMATION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Original STL Reference Plane:\n")
        f.write(f"  Method: {plane_orig['method']}\n")
        f.write(f"  Equation: {plane_orig['equation']}\n")
        f.write(f"  Normal Vector: [{plane_orig['normal'][0]:.6f}, {plane_orig['normal'][1]:.6f}, {plane_orig['normal'][2]:.6f}]\n")
        f.write(f"  Centroid: [{plane_orig['centroid'][0]:.3f}, {plane_orig['centroid'][1]:.3f}, {plane_orig['centroid'][2]:.3f}] km\n\n")
        
        # Interpretation guide
        f.write("INTERPRETATION GUIDE:\n")
        f.write("-" * 40 + "\n")
        f.write("Preservation Ratios:\n")
        f.write("  < 1.2: Excellent preservation\n")
        f.write("  1.2-1.5: Good preservation\n")
        f.write("  1.5-2.0: Moderate preservation\n")
        f.write("  > 2.0: Poor preservation\n\n")
        f.write("Smoothness Scores:\n")
        f.write("  Lower values indicate smoother surfaces\n")
        f.write("  Scores combine planar deviation + surface roughness + curvature\n\n")
        f.write("Surface Roughness:\n")
        f.write("  Measured using local neighborhood analysis\n")
        f.write("  RMS roughness provides overall surface texture measure\n\n")
    
    print(f"Detailed smoothness report saved to {output_filename}")

def save_grid_data(processing_data, output_folder_path):
    """
    Save all grid data and metadata to files in organized folder structure
    """
    grid_data = processing_data['grid_data'] 
    metadata = processing_data['metadata']
    
    # Generate file prefix from folder name
    folder_name = os.path.basename(output_folder_path)
    output_prefix = os.path.join(output_folder_path, folder_name)
    
    # Save CSV data
    csv_filename = f"{output_prefix}_grid.csv"
    grid_data.to_csv(csv_filename, index=False)
    print(f"Grid data saved to {csv_filename}")
    
    # Save Geometry.txt
    geometry_filename = f"{output_prefix}_Geometry.txt"
    save_geometry_txt(grid_data, metadata, geometry_filename)
    
    # Save numpy arrays and metadata
    npz_filename = f"{output_prefix}_grid.npz"
    np.savez(npz_filename,
             x=grid_data['x'].values,
             z=grid_data['z'].values, 
             y=grid_data['y'].values,
             dy_dx=grid_data['dy_dx'].values,
             dy_dz=grid_data['dy_dz'].values,
             metadata=metadata)
    print(f"Grid arrays and metadata saved to {npz_filename}")
    
    # Save detailed smoothness report
    smoothness_report_filename = f"{output_prefix}_smoothness_report.txt"
    save_smoothness_report(processing_data, smoothness_report_filename)
    
    # Create plots
    plot_processing_steps(processing_data, output_prefix)
    
    # Create smoothness analysis plots
    smoothness_data = metadata['smoothness_analysis']
    
    # Plot original STL smoothness
    plot_smoothness_analysis(
        smoothness_data['original_smoothness'], 
        processing_data['vertices_original'], 
        output_prefix, 
        title_suffix=" (Original STL)")
    
    # Plot processed mesh smoothness
    processed_vertices = np.column_stack([grid_data['x'], grid_data['y'], grid_data['z']])
    plot_smoothness_analysis(
        smoothness_data['processed_smoothness'], 
        processed_vertices, 
        f"{output_prefix}_processed", 
        title_suffix=" (Processed Mesh)")
    
    # Plot smoothness comparison
    plot_smoothness_comparison(smoothness_data['smoothness_comparison'], output_prefix)
    
    # Plot smoothing effects if smoothing was applied
    smoothing_control = metadata.get('smoothing_control', {})
    if smoothing_control.get('method', 'none') != 'none' and smoothing_control.get('smoothing_info'):
        # For smoothing visualization, we need to create a "before smoothing" version
        # Since we don't store the pre-smoothing processed vertices, we'll skip this visualization
        # or create a simplified version
        print("Note: Smoothing effects visualization skipped (would require storing intermediate results)")
        # plot_smoothing_effects could be called here if we stored pre-smoothing vertices
    
    # Plot vertical transition effects if transition was applied
    transition_data = metadata.get('vertical_transition', {})
    if transition_data.get('enabled', False) and transition_data.get('transition_info'):
        if 'error' not in transition_data['transition_info']:
            print("Note: Vertical transition effects visualization skipped (would require storing intermediate results)")
            # To properly show transition effects, we would need to store the pre-transition processed vertices
            # plot_vertical_transition_effects could be called here if we stored pre-transition vertices
    
    # Create summary report
    summary_filename = f"{output_prefix}_summary.txt"
    create_processing_summary(processing_data, summary_filename)
    
    # Copy Geometry.txt as bFault_Rough_Geometry.txt for convenience
    copy_geometry_as_bfault(processing_data, output_folder_path)
    
    # Print summary
    print(f"\nProcessing Summary:")
    print(f"Original STL: {metadata['stl_filename']}")
    print(f"Strike angle: {metadata['original_strike_angle']:.2f}°")
    print(f"Grid resolution: {metadata['resolution']} km")
    print(f"Grid size: {metadata['grid_dimensions'][0]} x {metadata['grid_dimensions'][1]}")
    print(f"STL zone points: {metadata['processing_info']['stl_zone_points']}")
    print(f"Extension zone points: {metadata['processing_info']['extension_zone_points']}")
    print(f"Total data points: {metadata['total_points']}")
    print(f"Y range: {metadata['y_range'][0]:.1f} to {metadata['y_range'][1]:.1f} km")
    
    # Print smoothness summary
    comparison = smoothness_data['smoothness_comparison']
    print(f"\nSmoothness Analysis:")
    print(f"Original smoothness score: {smoothness_data['original_smoothness']['overall_smoothness_score']:.4f}")
    print(f"Processed smoothness score: {smoothness_data['processed_smoothness']['overall_smoothness_score']:.4f}")
    print(f"Average preservation ratio: {comparison['average_preservation_ratio']:.3f}")
    print(f"Preservation quality: {comparison['preservation_quality'].upper()}")
    
    # Print smoothing summary if applied
    if smoothing_control.get('method', 'none') != 'none':
        print(f"\nSmoothing Applied:")
        print(f"Method: {smoothing_control['method']}")
        smoothing_info = smoothing_control.get('smoothing_info', {})
        if 'error' not in smoothing_info:
            if 'smoothing_factor' in smoothing_info:
                print(f"Smoothing factor: {smoothing_info['smoothing_factor']:.3f}")
    
    # Print transition summary if applied  
    if transition_data.get('enabled', False):
        print(f"\nVertical Transition Applied:")
        transition_info = transition_data.get('transition_info', {})
        if 'error' not in transition_info:
            print(f"Plane method: {transition_data['plane_method']}")
            print(f"Start depth: {transition_info['transition_start_depth']:.1f} km")
            print(f"Distance: {transition_info['transition_distance']:.1f} km")
            print(f"Points affected: {transition_info['transition_points']} ({transition_info['transition_fraction']:.1%})")
            print(f"Max Y change: {transition_info['max_y_change']:.6f} km")
        else:
            print(f"ERROR: {transition_info['error']}")

def create_processing_summary(processing_data, summary_filename):
    """
    Create a comprehensive processing summary report
    """
    metadata = processing_data['metadata']
    
    with open(summary_filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("STL TO FAULT GRID PROCESSING SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Processing Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"STL File: {metadata['stl_filename']}\n")
        f.write(f"Output Folder: {os.path.dirname(summary_filename)}\n\n")
        
        # Basic processing info
        f.write("BASIC PROCESSING PARAMETERS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Grid Resolution: {metadata['resolution']} km\n")
        f.write(f"Grid Dimensions: {metadata['grid_dimensions'][0]} x {metadata['grid_dimensions'][1]}\n")
        f.write(f"X Range: {metadata['x_range'][0]} to {metadata['x_range'][1]} km\n")
        f.write(f"Z Range: {metadata['z_range'][0]} to {metadata['z_range'][1]} km\n")
        f.write(f"Total Points: {metadata['total_points']}\n")
        f.write(f"Original Strike Angle: {metadata['original_strike_angle']:.2f}°\n\n")
        
        # Smoothing info
        smoothing_control = metadata.get('smoothing_control', {})
        f.write("SMOOTHING CONTROL:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Method: {smoothing_control.get('method', 'none')}\n")
        if smoothing_control.get('method', 'none') != 'none':
            params = smoothing_control.get('parameters', {})
            for key, value in params.items():
                f.write(f"{key}: {value}\n")
            
            smoothing_info = smoothing_control.get('smoothing_info', {})
            if 'error' not in smoothing_info and 'smoothing_factor' in smoothing_info:
                f.write(f"Smoothing factor achieved: {smoothing_info['smoothing_factor']:.3f}\n")
        f.write("\n")
        
        # Vertical transition info
        transition_data = metadata.get('vertical_transition', {})
        f.write("VERTICAL TRANSITION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Enabled: {transition_data.get('enabled', False)}\n")
        if transition_data.get('enabled', False):
            f.write(f"Plane Method: {transition_data.get('plane_method', 'N/A')}\n")
            f.write(f"Blend Function: {transition_data.get('blend_function', 'N/A')}\n")
            
            transition_info = transition_data.get('transition_info', {})
            if 'error' not in transition_info:
                f.write(f"Start Depth: {transition_info['transition_start_depth']:.1f} km\n")
                f.write(f"End Depth: {transition_info['transition_end_depth']:.1f} km\n")
                f.write(f"Distance: {transition_info['transition_distance']:.1f} km\n")
                f.write(f"Points in Transition: {transition_info['transition_points']}\n")
                f.write(f"Transition Fraction: {transition_info['transition_fraction']:.2%}\n")
                f.write(f"Max Y Change: {transition_info['max_y_change']:.6f} km\n")
                f.write(f"RMS Y Change: {transition_info['rms_y_change']:.6f} km\n")
                
                plane_info = transition_info['plane_info']
                f.write(f"Target Plane: {plane_info['description']}\n")
            else:
                f.write(f"ERROR: {transition_info['error']}\n")
        f.write("\n")
        
        # Smoothness analysis
        smoothness_data = metadata['smoothness_analysis']
        f.write("SMOOTHNESS ANALYSIS:\n")
        f.write("-" * 40 + "\n")
        
        orig_smooth = smoothness_data['original_smoothness']
        proc_smooth = smoothness_data['processed_smoothness']
        comparison = smoothness_data['smoothness_comparison']
        
        f.write(f"Original Smoothness Score: {orig_smooth['overall_smoothness_score']:.6f}\n")
        f.write(f"Processed Smoothness Score: {proc_smooth['overall_smoothness_score']:.6f}\n")
        f.write(f"Average Preservation Ratio: {comparison['average_preservation_ratio']:.3f}\n")
        f.write(f"Preservation Quality: {comparison['preservation_quality'].upper()}\n\n")
        
        # Output files
        f.write("OUTPUT FILES:\n")
        f.write("-" * 40 + "\n")
        folder_name = os.path.basename(os.path.dirname(summary_filename))
        f.write(f"• {folder_name}_Geometry.txt (main output for EQQuasi)\n")
        f.write(f"• {folder_name}_grid.csv (tabular data)\n")
        f.write(f"• {folder_name}_grid.npz (numpy arrays)\n")
        f.write(f"• {folder_name}_smoothness_report.txt (detailed smoothness analysis)\n")
        f.write(f"• {folder_name}_processing.png (processing visualization)\n")
        f.write(f"• {folder_name}_smoothness_analysis.png (original STL smoothness)\n")
        f.write(f"• {folder_name}_processed_smoothness_analysis.png (processed mesh smoothness)\n")
        f.write(f"• {folder_name}_smoothness_comparison.png (preservation comparison)\n")
        
        # Note: smoothing and transition effect plots are currently disabled due to 
        # vertex array size mismatches - would need intermediate result storage
        # if smoothing_control.get('method', 'none') != 'none':
        #     f.write(f"• {folder_name}_smoothing_effects.png (smoothing effects)\n")
        # 
        # if transition_data.get('enabled', False) and 'error' not in transition_data.get('transition_info', {}):
        #     f.write(f"• {folder_name}_vertical_transition.png (transition effects)\n")
        
        f.write(f"• {folder_name}_summary.txt (this file)\n")
    
    print(f"Processing summary saved to {summary_filename}")

def copy_geometry_as_bfault(processing_data, output_folder_path):
    """
    Copy the main Geometry.txt file as bFault_Rough_Geometry.txt for convenience
    """
    import shutil
    
    folder_name = os.path.basename(output_folder_path)
    geometry_filename = os.path.join(output_folder_path, f"{folder_name}_Geometry.txt")
    bfault_filename = os.path.join(output_folder_path, "bFault_Rough_Geometry.txt")
    
    if os.path.exists(geometry_filename):
        shutil.copy2(geometry_filename, bfault_filename)
        print(f"Copied Geometry.txt as bFault_Rough_Geometry.txt")
    
    return bfault_filename

if __name__ == "__main__":
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Enhanced STL to fault grid converter with smoothness analysis and control',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic processing with no smoothing
  python stl_to_eqquasi_fault_grid_mod.py fault.stl --resolution 1.0
  
  # Apply Gaussian smoothing
  python stl_to_eqquasi_fault_grid_mod.py fault.stl --smoothing gaussian --sigma 1.5
  
  # Apply adaptive smoothing to target roughness
  python stl_to_eqquasi_fault_grid_mod.py fault.stl --smoothing adaptive --target-roughness 0.05
  
  # Apply RBF smoothing
  python stl_to_eqquasi_fault_grid_mod.py fault.stl --smoothing rbf --sigma 2.0
  
  # Apply vertical transition to planar geometry
  python stl_to_eqquasi_fault_grid_mod.py fault.stl --vertical-transition --transition-distance 3.0
  
  # Combined smoothing and vertical transition
  python stl_to_eqquasi_fault_grid_mod.py fault.stl --smoothing gaussian --sigma 1.2 --vertical-transition --plane-method strike_aligned

Smoothing Methods:
  - none: No smoothing (default)
  - gaussian: Gaussian kernel smoothing (requires --sigma)
  - adaptive: Iterative smoothing to target roughness (requires --target-roughness)
  - rbf: Radial basis function smoothing (requires --sigma)
  - constrained: Optimization-based smoothing with constraints

Vertical Transition Methods:
  - best_fit_line: Fit line through fault geometry (default)
  - centroid_line: Vertical line through surface centroid
  - strike_aligned: Vertical plane aligned with strike direction
        """)
    
    parser.add_argument('stl_filename', help='Path to the STL file')
    parser.add_argument('--resolution', type=float, default=1.0, 
                       help='Grid resolution in km (default: 1.0)')
    parser.add_argument('--smoothing', choices=['none', 'gaussian', 'adaptive', 'rbf', 'constrained'],
                       default='none', help='Smoothing method to apply (default: none)')
    
    # Smoothing parameters
    parser.add_argument('--sigma', type=float, default=1.0,
                       help='Smoothing parameter sigma for gaussian/rbf methods (default: 1.0)')
    parser.add_argument('--target-roughness', type=float, default=0.1,
                       help='Target RMS roughness for adaptive smoothing (default: 0.1)')
    parser.add_argument('--max-iterations', type=int, default=10,
                       help='Maximum iterations for adaptive smoothing (default: 10)')
    parser.add_argument('--preserve-boundaries', action='store_true', default=True,
                       help='Preserve boundary values during smoothing (default: True)')
    
    # Vertical transition parameters
    parser.add_argument('--vertical-transition', action='store_true',
                       help='Enable smooth vertical transition to planar geometry')
    parser.add_argument('--transition-start-depth', type=float, default=None,
                       help='Depth where transition begins (km). If not specified, auto-determined from STL')
    parser.add_argument('--transition-distance', type=float, default=5.0,
                       help='Distance over which transition occurs (default: 5.0 km)')
    parser.add_argument('--plane-method', choices=['best_fit_line', 'centroid_line', 'strike_aligned'],
                       default='best_fit_line', help='Method for target vertical plane (default: best_fit_line)')
    parser.add_argument('--blend-function', choices=['linear', 'cosine', 'smooth_step', 'exponential'],
                       default='cosine', help='Blending function for transition (default: cosine)')
    
    args = parser.parse_args()
    
    # Prepare smoothing parameters
    smoothing_params = {}
    if args.smoothing in ['gaussian', 'rbf']:
        smoothing_params['sigma'] = args.sigma
        if args.smoothing == 'gaussian':
            smoothing_params['preserve_boundaries'] = args.preserve_boundaries
    elif args.smoothing == 'adaptive':
        smoothing_params['target_roughness'] = args.target_roughness
        smoothing_params['max_iterations'] = args.max_iterations
    
    # Configuration
    config = {
        'resolution': args.resolution,
        'x_range': (-380, 380),     # km
        'z_range': (0, -50),        # km (surface to depth)
        'max_vertices': 50000,      # limit for performance
        'smoothing_method': args.smoothing,
        'smoothing_params': smoothing_params if smoothing_params else None,
        'enable_vertical_transition': args.vertical_transition,
        'transition_start_depth': args.transition_start_depth,
        'transition_distance': args.transition_distance,
        'plane_method': args.plane_method,
        'blend_function': args.blend_function
    }
    
    # Create output folder with descriptive name
    output_folder_path = create_output_folder(args.stl_filename, config)
    
    print(f"Processing {args.stl_filename} with {args.resolution} km resolution...")
    print("Enhanced version with smoothness analysis and control enabled.")
    if args.smoothing != 'none':
        print(f"Smoothing method: {args.smoothing}")
        print(f"Smoothing parameters: {smoothing_params}")
    if args.vertical_transition:
        print(f"Vertical transition enabled:")
        print(f"  Plane method: {args.plane_method}")
        print(f"  Transition distance: {args.transition_distance} km")
        print(f"  Blend function: {args.blend_function}")
    
    # Process STL to grid
    processing_data = process_stl_to_grid(args.stl_filename, **config)
    
    # Save results to organized folder
    save_grid_data(processing_data, output_folder_path)
    
    # Get folder name for output summary
    folder_name = os.path.basename(output_folder_path)
    
    print(f"\nProcessing complete!")
    print(f"All output files saved to folder: {folder_name}/")
    print(f"Key output files:")
    print(f"  - {folder_name}_Geometry.txt (main output for EQQuasi)")
    print(f"  - {folder_name}_grid.csv (tabular data)")
    print(f"  - {folder_name}_smoothness_report.txt (detailed analysis)")
    print(f"  - {folder_name}_summary.txt (processing summary)")
    print(f"  - Various visualization plots (.png files)")
    
    # Print final results summary
    metadata = processing_data['metadata']
    
    if args.smoothing != 'none':
        smoothing_info = metadata['smoothing_control']['smoothing_info']
        print(f"\nSmoothing Results:")
        if 'error' in smoothing_info:
            print(f"  ERROR: {smoothing_info['error']}")
        else:
            print(f"  Method: {smoothing_info['method']}")
            if 'smoothing_factor' in smoothing_info:
                print(f"  Smoothing factor: {smoothing_info['smoothing_factor']:.3f}")
                print(f"  Original roughness: {smoothing_info['original_roughness']:.6f}")
                print(f"  Smoothed roughness: {smoothing_info['smoothed_roughness']:.6f}")
            if smoothing_info['method'] == 'adaptive_smoothing':
                print(f"  Target achieved: {smoothing_info['target_achieved']}")
                print(f"  Iterations used: {smoothing_info['iterations_used']}")
                print(f"  Final roughness: {smoothing_info['final_roughness']:.6f}")
    
    if args.vertical_transition:
        transition_info = metadata['vertical_transition']['transition_info']
        print(f"\nVertical Transition Results:")
        if 'error' in transition_info:
            print(f"  ERROR: {transition_info['error']}")
        else:
            print(f"  Plane: {transition_info['plane_info']['description']}")
            print(f"  Start depth: {transition_info['transition_start_depth']:.1f} km")
            print(f"  End depth: {transition_info['transition_end_depth']:.1f} km")
            print(f"  Points affected: {transition_info['transition_points']} ({transition_info['transition_fraction']:.1%})")
            print(f"  Max Y change: {transition_info['max_y_change']:.6f} km")
            print(f"  RMS Y change: {transition_info['rms_y_change']:.6f} km")