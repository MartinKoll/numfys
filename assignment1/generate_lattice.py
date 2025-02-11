import numpy as np
import matplotlib.pyplot as plt
from numba import jit, njit

@njit
def generate_square_lattice(L : int) -> list:
    N = L * L  # sites
    bonds = []

    for i in range(N):
        # Connect to the site to the right
        if (i + 1) % L == 0:
            bonds.append((i, i - L + 1))
        else:
            bonds.append((i, i + 1))
        if (i + L) >= N:
            bonds.append((i, i - N + L))
        else:
            bonds.append((i, i + L))
    return bonds, N, len(bonds)

def generate_triangular_lattice(L)  -> list:
    N = L * L  # Total number of nodes in the lattice
    bonds = []
    
    for i in range(N):
        x, y = i % L, i // L  # Get coordinates in the LxL grid
        
        right = ((x + 1) % L) + y * L  # Right neighbor (periodic boundary in x)
        down = x + ((y + 1) % L) * L  # Bottom neighbor (periodic boundary in y)
        diag = ((x + 1) % L) + ((y + 1) % L) * L  # Diagonal neighbor (triangular structure)
        
        bonds.append((i, right))
        bonds.append((i, down))
        bonds.append((i, diag))
    
    return bonds, N, len(bonds)


def generate_honeycomb_lattice(L:int) -> list:
    N = 2 * L * L  # Total number of nodes in the lattice
    bonds = set()  # Use a set to avoid duplicate bonds
    
    for y in range(L):
        for x in range(L):
            i = 2 * (y * L + x)  # First node in the hexagonal pair
            j = i + 1  # Second node in the hexagonal pair
            
            # Connect hexagon pair
            bonds.add((min(i, j), max(i, j)))
            
            # Connect to right neighbor (periodic in x)
            right = 2 * ((y * L + (x + 1) % L))
            bonds.add((min(j, right), max(j, right)))
            
            # Connect to bottom-right neighbor (periodic in y)
            bottom_right = 2 * (((y + 1) % L) * L + x) + 1
            bonds.add((min(j, bottom_right), max(j, bottom_right)))
    
    return list(bonds), N, len(bonds) 

def write_bonds_to_file(bonds):
    with open(f"./assignment1/bonds/honeycomb_lattice/bonds{bonds[1]}.txt", 'w') as f:
        f.write(f"{bonds[1]}\n")
        f.write(f"{bonds[2]}\n")
        for bond in bonds[0]:
            f.write(f"{bond[0]} {bond[1]}\n")
        f.close()
#test square lattice
#L = 1000
#bonds = generate_square_lattice(L)
#write_bonds_to_file(bonds)
# L(n) = 2n^2 square lattice

#test triangular lattice
#L = 1000
#bonds = generate_triangular_lattice(L)
#write_bonds_to_file(bonds)
# L(n) = 3n^2 triangular lattice

#test honeycomb lattice
L = 6
bonds = generate_honeycomb_lattice(L)
write_bonds_to_file(bonds)
#L(n) = ?

print("done")