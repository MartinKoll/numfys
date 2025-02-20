import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
pre = "square"
prefix = f"./assignment1/convolution_results/{pre}/"
txt_list = ["bonds10000.txt_convolution.npz", "bonds40000.txt_convolution.npz", "bonds90000.txt_convolution.npz", "bonds160000.txt_convolution.npz", "bonds250000.txt_convolution.npz", "bonds360000.txt_convolution.npz", "bonds490000.txt_convolution.npz", "bonds640000.txt_convolution.npz", "bonds810000.txt_convolution.npz", "bonds1000000.txt_convolution.npz"]

def load_results_convolution(file):
    '''Load results from convolution'''
    results = np.load(file)
    q_values = results["q_values"]
    Q_p_inf = results["Q_p_inf"]
    Q_suspectability = results["Q_suspectability"]
    Q_avg_cluster_size = results["Q_avg_cluster_size"]
    Q_p_inf2 = results["Q_p_inf2"]
    return q_values, Q_p_inf, Q_suspectability, Q_avg_cluster_size, Q_p_inf2

def plot_results(p, P1, s_avg, susceptibility, L, ax):
    mask = (p >= 0.4) & (p <= 0.60) #square
    ax.plot(p[mask], P1[mask], label=f'L={L}')
    ax.set_xlabel('p')
    ax.set_ylabel('P1')
    ax.set_title('Giant Component')
    ax.legend()

def plot_avg_cluster_size(p, s_avg, L, ax):
    mask = (p >= 0.4) & (p <= 0.60) #square
    ax.plot(p[mask], s_avg[mask], label=f'L={L}')
    ax.set_xlabel('p')
    ax.set_ylabel('<s>')
    ax.set_title('Average Cluster Size')
    ax.legend()

def plot_susceptibility(p, susceptibility, L, ax):
    mask = (p >= 0.4)  #square
    ax.plot(p[mask], susceptibility[mask], label=f'L={L}')
    ax.set_xlabel('p')
    ax.set_ylabel('Susceptibility')
    ax.set_title('Susceptibility')
    ax.legend()

#fig1, ax1 = plt.subplots(figsize=(10, 6))
#fig2, ax2 = plt.subplots(figsize=(10, 6))
# fig3, ax3 = plt.subplots(figsize=(10, 6))

# for txt in txt_list:
#     q_values, Q_p_inf, Q_suspectability, Q_avg_cluster_size, Q_p_inf2 = load_results_convolution(prefix+txt) #load results
#     L = int(txt.split(".")[0][5:])  # Extract L from filename
#     #plot_results(q_values, Q_p_inf, Q_avg_cluster_size, Q_suspectability, L, ax1) # Plot results
#     #plot_avg_cluster_size(q_values, Q_avg_cluster_size, L, ax2) # Plot average cluster size
#     plot_susceptibility(q_values, Q_suspectability, L, ax3) # Plot susceptibility
#     print(f"{txt} done")

#plt.savefig("./assignment1/pics/square/suspectability.png", dpi=300)

#plt.tight_layout()
#plt.show()

# Plot p_inf for hexagonal, triangular, and square lattice
lattices = ["hexagonal", "triangular", "square"]
fig, ax = plt.subplots(1, 3, figsize=(12, 6))

for i, lattice in enumerate(lattices):
    prefix = f"./assignment1/convolution_results/{lattice}/"
    for txt in txt_list:
        q_values, Q_p_inf, Q_suspectability, Q_avg_cluster_size, Q_p_inf2 = load_results_convolution(prefix+txt) #load results
        L = int(txt.split(".")[0][5:])  # Extract L from filename
        ax[i].plot(q_values, Q_p_inf, label=f'{lattice.capitalize()} L={L}')
        print(f"{lattice} {txt} done")


plt.savefig("./assignment1/pics/comparison/p_inf_comparison.png", dpi=300)

#plt.tight_layout()
plt.show()