import numpy as np
from numba import njit
from scipy.special import comb


def ten_randints(seed : int = None, low:int = 0, high:int = 100):
    np.random.seed(seed)
    randints = np.random.randint(low, high, size=10)
    return randints


@njit(cache=True)
def find_root(sites:np.array, index:int):
    root = index
    while sites[root] >= 0:
        root = sites[root]
    while index != root: #compress path
        parent = sites[index]
        sites[index] = root
        index = parent
    return root

@njit(cache=True)
def process_bonds(bonds:np.array, sites:np.array, N:int):
    average_s = N  # Initial sum of squares (each cluster size 1)
    largest_cluster_size = 1
    #largest_cluster_index = 0
    steps = []
    p_inf = []
    p_inf2 = []
    avg_cluster_size = []
    total_bonds = bonds.shape[0]
    
    for i in range(total_bonds):
        bond = bonds[i]
        root1 = find_root(sites, bond[0])
        root2 = find_root(sites, bond[1])
        
        if root1 != root2:
            s1 = -sites[root1]
            s2 = -sites[root2]
            average_s -= (s1**2 + s2**2)
            combined_size = s1 + s2
            average_s += combined_size**2
            
            # Union by size
            if s1 >= s2:
                sites[root1] += sites[root2]
                sites[root2] = root1
                if combined_size > largest_cluster_size:
                    largest_cluster_size = combined_size
            else:
                sites[root2] += sites[root1]
                sites[root1] = root2
                if combined_size > largest_cluster_size:
                    largest_cluster_size = combined_size

        steps.append(i/total_bonds)
        p_inf.append(largest_cluster_size/N)
        p_inf2.append((largest_cluster_size/N)**2)
        if largest_cluster_size == N:
            avg = 0.0
        else:
            if N*(1 - p_inf[-1]) <= 0.000000000001 or N*p_inf[-1] <= 0.000000000001:
                avg = 0.0
            else:
                avg = (average_s - (N * p_inf[-1])**2) / (N * (1 - p_inf[-1]))
        avg_cluster_size.append(avg)
    return np.array(steps), np.array(p_inf), np.array(p_inf2), np.array(avg_cluster_size)

def load_results_percolation(filename):
    """Load results from an .npz file."""
    data = np.load(filename)
    return data["steps"], data["p_inf"], data["p_inf2"], data["suspect"], data["avg_cluster_size"]

def load_results_convolution(filename):
    """Load results from an .npz file."""
    data = np.load(filename)
    return data["q_values"], data["Q_p_inf"], data["Q_suspectability"], data["Q_avg_cluster_size"], data["Q_p_inf2"]

def suspectability(p_inf, p_inf2, N):
    suspectability = np.zeros_like(p_inf)
    for i in range(len(p_inf)):
        p_inf_temp = np.mean(p_inf[:i+1])
        p_inf2_temp = np.mean(p_inf2[:i+1])
        diff = p_inf2_temp - p_inf_temp**2
        suspectability[i] = N*np.sqrt(diff) if diff > 0 else 0.0
    return suspectability

