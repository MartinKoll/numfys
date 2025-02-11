import numpy as np
import matplotlib.pyplot as plt
from numba import jit, njit
from matplotlib.colors import ListedColormap
from tqdm import tqdm
from bonds_class import bonds
from functions import find_root, process_bonds, suspectability
import os

output_dir = "./assignment1/data/triangular/"
os.makedirs(output_dir, exist_ok=True)

#["bonds10000.txt", "bonds40000.txt", "bonds90000.txt", "bonds160000.txt", "bonds250000.txt", "bonds360000.txt", "bonds490000.txt", "bonds640000.txt", "bonds810000.txt", "bonds1000000.txt"]
txt_list = ["bonds10000.txt", "bonds40000.txt", "bonds90000.txt", "bonds160000.txt", "bonds250000.txt", "bonds360000.txt", "bonds490000.txt", "bonds640000.txt", "bonds810000.txt", "bonds1000000.txt"]

iterations = 1000

fig, ax = plt.subplots(2, 2, figsize=(12, 10))


for txt in txt_list:
    bonds1 = bonds(txt, type="triangular")
    bonbonds = bonds1.bonds
    sites = bonds1.sites
    N = bonds1.N

    steps_arr = np.zeros(bonbonds.shape[0])
    p_inf_arr = np.zeros(bonbonds.shape[0])
    p_inf2_arr = np.zeros(bonbonds.shape[0])
    avg_cluster_size_arr = np.zeros(bonbonds.shape[0])

    for _ in tqdm(range(iterations)):
        steps, p_inf, p_inf2, avg_cluster_size = process_bonds(bonbonds, sites.copy(), N)
        steps_arr += steps
        p_inf_arr += p_inf
        p_inf2_arr += p_inf2
        avg_cluster_size_arr += avg_cluster_size

    steps_avg = steps_arr / iterations
    p_inf_avg = p_inf_arr / iterations
    p_inf2_avg = p_inf2_arr / iterations
    suspect = suspectability(p_inf_avg, p_inf2_avg, N)
    avg_cluster_size_avg = avg_cluster_size_arr / iterations

    # Save results in a single .npz file
    np.savez(os.path.join(output_dir, f"{txt}_results.npz"),
             steps=steps_avg, p_inf=p_inf_avg, p_inf2=p_inf2_avg,
             suspect=suspect, avg_cluster_size=avg_cluster_size_avg)
    print(f"Results saved in {output_dir}{txt}_results.npz")
    print("completed")

#     ax[0][0].plot(steps_avg, p_inf_avg, label=f"{N}")
#     ax[0][1].plot(steps_avg, avg_cluster_size_avg, label=f"{N}")
#     ax[1][0].plot(steps_avg, suspect, label=f"{N}")
#     ax[1][1].plot(steps_avg, p_inf2_avg, label=f"{N}")


# ax[0][0].set_title("Percolation Probability")
# ax[0][0].legend()

# ax[0][1].set_title("Average Cluster Size")
# ax[0][1].legend()

# ax[1][0].set_title("Susceptibility")
# ax[1][0].legend()

# ax[1][1].set_title("Percolation Probability (P_inf2)")
# ax[1][1].legend()

# plt.savefig("./assignment1/pics/p_avg.png", dpi=500)
# plt.show()




