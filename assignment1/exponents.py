import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

prefix = "./assignment1/convolution_results/square/"
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

def plot_results(p, P1, s_avg, susceptibility, L):
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 3, 1)
    plt.plot(p, P1, label=f'L={L}')
    plt.xlabel('p')
    plt.ylabel('P1')
    plt.title('Giant Component')

    plt.subplot(1, 3, 2)
    plt.plot(p, s_avg, label=f'L={L}')
    plt.xlabel('p')
    plt.ylabel('<s>')
    plt.title('Average Cluster Size')

    plt.subplot(1, 3, 3)
    plt.plot(p, susceptibility, label=f'L={L}')
    plt.xlabel('p')
    plt.ylabel('Susceptibility')
    plt.title('Susceptibility')

for txt in txt_list:
    q_values, Q_p_inf, Q_suspectability, Q_avg_cluster_size, Q_p_inf2 = load_results_convolution(prefix+txt) #load results
    L = int(txt.split("_")[0][5:])  # Extract L from filename
    plot_results(q_values, Q_p_inf, Q_avg_cluster_size, Q_suspectability, L) # Plot results

plt.show()