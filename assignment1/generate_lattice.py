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

def honeycomb_lattice(L):
    N = L * L
    bonds = set()

    for i in range(N):
        if (i // L) % 2 == 0:
            if i % 2 == 0:
                if (i + 1) % L != 0:
                    bonds.add((min(i, i + 1), max(i, i + 1)))
                else:
                    bonds.add((min(i, i - L + 1), max(i, i - L + 1)))
        else:
            if i % 2:
                if (i + 1) % L != 0:
                    bonds.add((min(i, i + 1), max(i, i + 1)))
                else:
                    bonds.add((min(i, i - L + 1), max(i, i - L + 1)))
        
        if (i + L) < N:
            bonds.add((min(i, i + L), max(i, i + L)))
        else:
            bonds.add((min(i, i - N + L), max(i, i - N + L)))

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
#L(n) = 1.5 n
L=[100,200,300,400,500,600,700,800,900,1000]
for l in L:
    bonds = honeycomb_lattice(l)
    write_bonds_to_file(bonds)

print("done")