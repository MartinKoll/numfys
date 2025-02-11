import matplotlib.pyplot as plt
import numpy as np

def generate_honeycomb_lattice(L):
    N = 2 * L * L  # Total number of nodes in the lattice
    bonds = set()  # Use a set to avoid duplicate bonds
    positions = {}  # Store node positions for visualization
    
    for y in range(L):
        for x in range(L):
            i = 2 * (y * L + x)  # First node in the hexagonal pair
            j = i + 1  # Second node in the hexagonal pair
            
            # Compute positions
            positions[i] = (x + (y % 2) * 0.5, y * np.sqrt(3) / 2)
            positions[j] = (x + (y % 2) * 0.5 + 0.5, y * np.sqrt(3) / 2)
            
            # Connect hexagon pair
            bonds.add((min(i, j), max(i, j)))
            
            # Connect to right neighbor (periodic in x)
            right = 2 * ((y * L + (x + 1) % L))
            bonds.add((min(j, right), max(j, right)))
            
            # Connect to bottom-right neighbor (periodic in y)
            bottom_right = 2 * (((y + 1) % L) * L + x) + 1
            bonds.add((min(j, bottom_right), max(j, bottom_right)))
    
    return N, len(bonds), sorted(bonds), positions

def plot_honeycomb_lattice(L):
    _, _, bonds, positions = generate_honeycomb_lattice(L)
    
    plt.figure(figsize=(8, 6))
    
    # Draw bonds
    for bond in bonds:
        x_values = [positions[bond[0]][0], positions[bond[1]][0]]
        y_values = [positions[bond[0]][1], positions[bond[1]][1]]
        plt.plot(x_values, y_values, 'k-', lw=1)
    
    # Draw nodes
    for node, (x, y) in positions.items():
        plt.scatter(x, y, color='blue', s=30)
    
    plt.axis('equal')
    plt.axis('off')
    plt.show()

def write_lattice_to_file(L, filename="honeycomb_lattice.txt"):
    N, num_bonds, bonds, _ = generate_honeycomb_lattice(L)
    
    with open(filename, "w") as f:
        f.write(f"{N} {num_bonds}\n")  # First line: number of nodes and number of bonds
        for bond in bonds:
            f.write(f"{bond[0]} {bond[1]}\n")
    
    print(f"Honeycomb lattice of size {L}x{L} written to {filename}")

# Example usage
L = 4  # Change this for different sizes
write_lattice_to_file(L)
plot_honeycomb_lattice(L)
