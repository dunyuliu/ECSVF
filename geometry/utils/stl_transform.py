import numpy as np
from stl import mesh
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

class STLTransformer:
    def __init__(self, stl_file):
        self.stl_file = stl_file
        self.mesh = mesh.Mesh.from_file(stl_file)
        self.vertices = self._extract_unique_vertices()
        self.top_trace = None
        self.rotation_matrix = None
        self.translation_vector = None
        
    def _extract_unique_vertices(self):
        """Extract unique vertices from the STL mesh"""
        all_vertices = self.mesh.vectors.reshape(-1, 3)
        unique_vertices = np.unique(all_vertices, axis=0)
        return unique_vertices
    
    def find_top_surface_trace(self, z_threshold_percentile=95):
        """Find vertices on the top surface trace"""
        z_values = self.vertices[:, 2]
        z_threshold = np.percentile(z_values, z_threshold_percentile)
        
        # Get vertices near the top surface
        top_indices = np.where(z_values >= z_threshold)[0]
        top_vertices = self.vertices[top_indices]
        
        # Find the perimeter/trace by selecting vertices with fewer neighbors
        # This is a simplified approach - you might need to adjust based on your specific geometry
        trace_vertices = []
        
        for vertex in top_vertices:
            # Count neighbors within a small radius
            distances = cdist([vertex], top_vertices)[0]
            neighbor_count = np.sum(distances < np.std(distances) * 0.5)
            
            # Vertices on the trace typically have fewer neighbors
            if neighbor_count < len(top_vertices) * 0.3:
                trace_vertices.append(vertex)
        
        if not trace_vertices:
            # Fallback: use all top vertices
            trace_vertices = top_vertices
            
        self.top_trace = np.array(trace_vertices)
        return self.top_trace
    
    def calculate_average_strike(self):
        """Calculate the average strike direction of the top trace"""
        if self.top_trace is None:
            self.find_top_surface_trace()
        
        # Use PCA to find the principal direction
        centered_trace = self.top_trace - np.mean(self.top_trace, axis=0)
        covariance_matrix = np.cov(centered_trace.T)
        eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
        
        # The first principal component is the strike direction
        principal_direction = eigenvectors[:, np.argmax(eigenvalues)]
        
        # Calculate strike angle (azimuth from north)
        strike_angle = np.arctan2(principal_direction[1], principal_direction[0])
        
        return strike_angle, principal_direction
    
    def create_transformation_matrix(self):
        """Create transformation matrix to align with x-axis and center at origin"""
        strike_angle, _ = self.calculate_average_strike()
        
        # Rotation matrix to align strike with x-axis
        # We need to rotate by -strike_angle to align with x-axis
        rotation_angle = -strike_angle
        cos_theta = np.cos(rotation_angle)
        sin_theta = np.sin(rotation_angle)
        
        # 3D rotation matrix around z-axis
        self.rotation_matrix = np.array([
            [cos_theta, -sin_theta, 0],
            [sin_theta, cos_theta, 0],
            [0, 0, 1]
        ])
        
        # Translation to center at origin (only x,y - keep z unchanged)
        centroid = np.mean(self.vertices, axis=0)
        self.translation_vector = -centroid
        self.translation_vector[2] = 0  # Don't translate z-coordinate
        
        return self.rotation_matrix, self.translation_vector
    
    def transform_coordinates(self, coordinates):
        """Transform given coordinates using the same transformation"""
        if self.rotation_matrix is None or self.translation_vector is None:
            self.create_transformation_matrix()
        
        coordinates = np.array(coordinates)
        
        # Handle both single point and multiple points
        if coordinates.ndim == 1:
            coordinates = coordinates.reshape(1, -1)
        
        # Apply translation first, then rotation
        translated = coordinates + self.translation_vector
        transformed = np.dot(translated, self.rotation_matrix.T)
        
        return transformed
    
    def transform_all_vertices(self):
        """Transform all vertices for visualization"""
        return self.transform_coordinates(self.vertices)
    
    def visualize_transformation(self):
        """Visualize original and transformed geometry"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Original geometry
        ax1.scatter(self.vertices[:, 0], self.vertices[:, 1], 
                   c=self.vertices[:, 2], cmap='viridis', s=1)
        if self.top_trace is not None:
            ax1.scatter(self.top_trace[:, 0], self.top_trace[:, 1], 
                       c='red', s=10, alpha=0.7, label='Top trace')
        ax1.set_xlabel('X (original)')
        ax1.set_ylabel('Y (original)')
        ax1.set_title('Original Geometry')
        ax1.axis('equal')
        ax1.legend()
        
        # Transformed geometry
        transformed_vertices = self.transform_all_vertices()
        ax2.scatter(transformed_vertices[:, 0], transformed_vertices[:, 1], 
                   c=transformed_vertices[:, 2], cmap='viridis', s=1)
        
        if self.top_trace is not None:
            transformed_trace = self.transform_coordinates(self.top_trace)
            ax2.scatter(transformed_trace[:, 0], transformed_trace[:, 1], 
                       c='red', s=10, alpha=0.7, label='Top trace')
        
        ax2.set_xlabel('X (transformed)')
        ax2.set_ylabel('Y (transformed)')
        ax2.set_title('Transformed Geometry (aligned & centered)')
        ax2.axis('equal')
        ax2.legend()
        
        plt.tight_layout()
        plt.show()

# Example usage
if __name__ == "__main__":
    # Initialize transformer
    transformer = STLTransformer("no_dip_change_1km.stl")
    
    # Process the geometry
    print("Finding top surface trace...")
    top_trace = transformer.find_top_surface_trace()
    print(f"Found {len(top_trace)} vertices on top trace")
    
    print("Calculating strike direction...")
    strike_angle, strike_direction = transformer.calculate_average_strike()
    print(f"Average strike angle: {np.degrees(strike_angle):.2f} degrees")
    print(f"Strike direction vector: {strike_direction}")
    
    print("Creating transformation matrix...")
    rotation_matrix, translation_vector = transformer.create_transformation_matrix()
    print(f"Translation vector: {translation_vector}")
    print(f"Rotation matrix:\n{rotation_matrix}")
    
    # Example coordinate transformation
    # Replace these with your actual coordinates
    original_coordinates = [
        [100, 200, 50],
        [150, 250, 75],
        [200, 180, 60]
    ]
    
    transformed_coordinates = transformer.transform_coordinates(original_coordinates)
    
    print("\nCoordinate transformation results:")
    for i, (orig, trans) in enumerate(zip(original_coordinates, transformed_coordinates)):
        print(f"Point {i+1}: {orig} -> {trans}")
    
    # Visualize the transformation
    print("\nGenerating visualization...")
    transformer.visualize_transformation()
    
    # Function to transform new coordinates
    def transform_new_coordinates(coords_list):
        """Helper function to transform a list of coordinates"""
        return transformer.transform_coordinates(coords_list)
    
    print("\nTo transform new coordinates, use:")
    print("transformed_coords = transformer.transform_coordinates(your_coordinates)")