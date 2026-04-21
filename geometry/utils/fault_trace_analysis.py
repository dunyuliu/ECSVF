#!/usr/bin/env python3
"""
Script to analyze fault trace from STL file, find top trace at depth=0,
recenter and rotate it to align with strike direction, and transform station coordinates.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

class FaultTraceAnalyzer:
    def __init__(self, stl_file):
        self.stl_file = stl_file
        self.mesh = None
        self.top_trace = None
        self.strike_angle = None
        self.center = None
        self.rotation_matrix = None
        
    def load_stl(self):
        """Load the STL file and extract vertices only"""
        print(f"Loading STL file: {self.stl_file}")
        
        # Use numpy to read STL file more efficiently  
        with open(self.stl_file, 'rb') as f:
            # Skip header (80 bytes)
            f.read(80)
            # Read number of triangles (4 bytes)
            num_triangles = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            print(f"Reading {num_triangles} triangles...")
            
            # Read all triangle data at once
            # Each triangle: normal(3*4) + vertex1(3*4) + vertex2(3*4) + vertex3(3*4) + attribute(2) = 50 bytes
            triangle_data = f.read(num_triangles * 50)
            
        # Parse triangle data efficiently
        vertices = []
        offset = 0
        for i in range(num_triangles):
            # Skip normal vector (12 bytes)
            offset += 12
            
            # Read 3 vertices (36 bytes)
            for j in range(3):
                vertex_bytes = triangle_data[offset:offset+12]
                if len(vertex_bytes) == 12:
                    x, y, z = np.frombuffer(vertex_bytes, dtype=np.float32)
                    vertices.append([x, y, z])
                offset += 12
            
            # Skip attribute byte count (2 bytes) 
            offset += 2
            
            # Progress indicator
            if (i + 1) % 50000 == 0:
                print(f"  Processed {i+1}/{num_triangles} triangles...")
        
        self.vertices = np.array(vertices)
        print(f"Loaded {len(self.vertices)} vertices")
        
    def find_top_trace(self):
        """Find the top trace (surface trace where depth = 0, i.e., maximum z coordinates)"""
        print("Finding top trace...")
        
        # Find unique vertices using numpy (much faster)
        print("Removing duplicate vertices...")
        # Round vertices to remove very close duplicates (1mm precision)
        rounded_vertices = np.round(self.vertices, 3)  
        unique_vertices, unique_indices = np.unique(rounded_vertices, axis=0, return_index=True)
        print(f"Found {len(unique_vertices)} unique vertices (from {len(self.vertices)} total)")
        
        # Find the maximum z-coordinate (surface)
        max_z = np.max(unique_vertices[:, 2])
        z_tolerance = 100  # meters tolerance for "surface" points
        
        print(f"Maximum Z coordinate: {max_z:.2f}")
        print(f"Using Z tolerance: {z_tolerance:.2f} meters")
        
        # Get vertices near the surface (top trace)
        surface_mask = unique_vertices[:, 2] >= (max_z - z_tolerance)
        self.top_trace = unique_vertices[surface_mask]
        
        print(f"Found {len(self.top_trace)} surface points")
        print(f"Z range of surface points: {np.min(self.top_trace[:, 2]):.2f} to {np.max(self.top_trace[:, 2]):.2f}")
        
    def calculate_strike_direction(self):
        """Calculate the average strike direction of the trace"""
        print("Calculating strike direction...")
        
        # Sort points by x-coordinate to get a rough along-strike ordering
        sorted_indices = np.argsort(self.top_trace[:, 0])
        sorted_trace = self.top_trace[sorted_indices]
        
        # Calculate vector differences between consecutive points
        vectors = np.diff(sorted_trace[:, :2], axis=0)  # Only x,y components
        
        # Calculate angles for each segment
        angles = np.arctan2(vectors[:, 1], vectors[:, 0])
        
        # Use circular statistics to find average angle
        # Convert to unit vectors and average
        unit_vectors = np.column_stack([np.cos(angles), np.sin(angles)])
        mean_vector = np.mean(unit_vectors, axis=0)
        
        # Get the average strike angle
        self.strike_angle = np.arctan2(mean_vector[1], mean_vector[0])
        
        print(f"Average strike angle: {np.degrees(self.strike_angle):.2f} degrees")
        
    def recenter_and_rotate(self):
        """Recenter the trace to origin and rotate to align strike with x-axis"""
        print("Recentering and rotating trace...")
        
        # Calculate center (centroid) of the trace
        self.center = np.mean(self.top_trace, axis=0)
        print(f"Trace center: x={self.center[0]:.2f}, y={self.center[1]:.2f}, z={self.center[2]:.2f}")
        
        # Create 2D rotation matrix to align strike with x-axis
        cos_angle = np.cos(-self.strike_angle)  # Negative to rotate strike to x-axis
        sin_angle = np.sin(-self.strike_angle)
        
        self.rotation_matrix = np.array([
            [cos_angle, -sin_angle],
            [sin_angle, cos_angle]
        ])
        
        # Apply transformation to trace
        centered_trace = self.top_trace - self.center
        self.rotated_trace = centered_trace.copy()
        self.rotated_trace[:, :2] = (self.rotation_matrix @ centered_trace[:, :2].T).T
        
        print("Transformation completed")
        
    def plot_traces(self):
        """Plot original and transformed traces"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot original trace
        ax1.scatter(self.top_trace[:, 0], self.top_trace[:, 1], c=self.top_trace[:, 2], 
                   cmap='viridis', s=20, alpha=0.7)
        ax1.set_xlabel('X (m)')
        ax1.set_ylabel('Y (m)')
        ax1.set_title('Original Top Trace')
        ax1.axis('equal')
        ax1.grid(True, alpha=0.3)
        
        # Add center point
        ax1.plot(self.center[0], self.center[1], 'r*', markersize=10, label='Center')
        
        # Add strike direction arrow
        arrow_length = np.ptp(self.top_trace[:, 0]) * 0.2
        dx = arrow_length * np.cos(self.strike_angle)
        dy = arrow_length * np.sin(self.strike_angle)
        ax1.arrow(self.center[0], self.center[1], dx, dy, 
                 head_width=arrow_length*0.1, head_length=arrow_length*0.1, 
                 fc='red', ec='red', label='Strike Direction')
        ax1.legend()
        
        # Plot rotated trace
        ax2.scatter(self.rotated_trace[:, 0], self.rotated_trace[:, 1], 
                   c=self.rotated_trace[:, 2], cmap='viridis', s=20, alpha=0.7)
        ax2.set_xlabel('X (m) - Strike Direction')
        ax2.set_ylabel('Y (m) - Dip Direction')
        ax2.set_title('Recentered and Rotated Top Trace')
        ax2.axis('equal')
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5, label='Strike Line')
        ax2.axvline(x=0, color='r', linestyle='--', alpha=0.5, label='Center')
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig('fault_trace_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def transform_coordinates(self, coordinates):
        """
        Transform a list of coordinates using the same transformation as the trace
        
        Parameters:
        coordinates: array-like of shape (n, 2) or (n, 3) with [x, y] or [x, y, z]
        
        Returns:
        transformed_coords: array with same shape as input
        """
        coordinates = np.array(coordinates)
        if coordinates.ndim == 1:
            coordinates = coordinates.reshape(1, -1)
            
        # Handle both 2D and 3D coordinates
        if coordinates.shape[1] == 2:
            # Add zero z-coordinate
            coords_3d = np.column_stack([coordinates, np.zeros(len(coordinates))])
        else:
            coords_3d = coordinates.copy()
            
        # Apply same transformation as trace
        centered_coords = coords_3d - self.center
        transformed_coords = centered_coords.copy()
        transformed_coords[:, :2] = (self.rotation_matrix @ centered_coords[:, :2].T).T
        
        # Return same format as input
        if coordinates.shape[1] == 2:
            return transformed_coords[:, :2]
        else:
            return transformed_coords
    
    def run_analysis(self):
        """Run complete analysis"""
        self.load_stl()
        self.find_top_trace()
        self.calculate_strike_direction()
        self.recenter_and_rotate()
        self.plot_traces()
        
        print("\nAnalysis complete!")
        print(f"Transformation parameters:")
        print(f"  Center: ({self.center[0]:.2f}, {self.center[1]:.2f}, {self.center[2]:.2f})")
        print(f"  Strike angle: {np.degrees(self.strike_angle):.2f} degrees")
        print(f"  Rotation matrix: \n{self.rotation_matrix}")

def main():
    # Initialize analyzer
    analyzer = FaultTraceAnalyzer('no_dip_change_1km.stl')
    
    # Run the analysis
    analyzer.run_analysis()
    
    # Station coordinates from the provided image
    # Based on the table, extracting the coordinates (assuming these are x, y coordinates)
    station_coords = [
        [1382372.944, 5199926.971],  # Station 1
        [1347927.726, 5017777.326],  # Station 2
        [1536168.087, 5311181.708]   # Station 3
    ]
    
    print(f"\nOriginal station coordinates:")
    for i, coord in enumerate(station_coords, 1):
        print(f"  Station {i}: ({coord[0]:.3f}, {coord[1]:.3f})")
    
    # Transform station coordinates
    transformed_stations = analyzer.transform_coordinates(station_coords)
    
    print(f"\nTransformed station coordinates:")
    for i, coord in enumerate(transformed_stations, 1):
        print(f"  Station {i}: ({coord[0]:.3f}, {coord[1]:.3f})")
    
    # Plot stations on both coordinate systems
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Original coordinates
    ax1.scatter(analyzer.top_trace[:, 0], analyzer.top_trace[:, 1], 
               c='lightblue', s=10, alpha=0.5, label='Fault trace')
    for i, coord in enumerate(station_coords, 1):
        ax1.plot(coord[0], coord[1], 'ro', markersize=8)
        ax1.annotate(f'S{i}', (coord[0], coord[1]), xytext=(5, 5), 
                    textcoords='offset points', fontsize=10, fontweight='bold')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_title('Original Coordinates')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Transformed coordinates
    ax2.scatter(analyzer.rotated_trace[:, 0], analyzer.rotated_trace[:, 1], 
               c='lightblue', s=10, alpha=0.5, label='Fault trace')
    for i, coord in enumerate(transformed_stations, 1):
        ax2.plot(coord[0], coord[1], 'ro', markersize=8)
        ax2.annotate(f'S{i}', (coord[0], coord[1]), xytext=(5, 5), 
                    textcoords='offset points', fontsize=10, fontweight='bold')
    ax2.set_xlabel('X (m) - Strike Direction')
    ax2.set_ylabel('Y (m) - Dip Direction')
    ax2.set_title('Transformed Coordinates')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax2.axvline(x=0, color='r', linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('station_coordinates_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return analyzer, transformed_stations

if __name__ == "__main__":
    analyzer, transformed_stations = main()