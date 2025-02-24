import numpy as np
from numba import jit
import matplotlib.pyplot as plt

def koch_generator(start:np.array, end:np.array)->np.array:
    start = np.array(start)
    end = np.array(end)
    
    # Compute segment vectors
    dx, dy = (end - start) / 4  # Step size

    # Define the 5 key points in order
    p1 = start
    p2 = start + np.array([dx, dy])  # First segment
    p3 = p2 + np.array([-dy, dx])    # Raise second part
    p4 = p3 + np.array([dx, dy])     # Lower third part
    p5 = p4 + 2*np.array([dy, -dx])    # Lower third part
    p6 = p5 + np.array([dx, dy])
    p7 = p6 + np.array([-dy, dx])
    p8 = end

    return [p1, p2, p3, p4, p5, p6, p7, p8]

def generate_fractal(level):
    """Generate the quadratic Koch fractal for a given level `ℓ`."""
    L = 1  # Side length of square
    square = [
        ([0, 0], [L, 0]),   # Bottom
        ([L, 0], [L, L]),   # Right
        ([L, L], [0, L]),   # Top
        ([0, L], [0, 0])    # Left
    ]
    
    new_edges = []
    for edge in square:
        new_edges.extend(koch_generator(*edge))  # Apply generator to each side

    return new_edges

# Generate fractal for l = 1
corners = generate_fractal(level=1)
# Convert list of points to x and y coordinates for plotting
x, y = zip(*corners)

# Plot the fractal
plt.figure(figsize=(6, 6))
plt.plot(x, y, 'k-', lw=1)
plt.axis("equal")
plt.title("Quadratic Koch Fractal (ℓ = 1)")
plt.show()
