import numpy as np
from scipy.special import gammaln
from functions import load_results_percolation
import matplotlib.pyplot as plt
from tqdm import tqdm
lattice = "hexagonal"
#square, triangular finished
txt_list = ["bonds1000000.txt"]
# "bonds10000.txt", "bonds40000.txt", "bonds90000.txt", "bonds160000.txt", "bonds250000.txt", "bonds360000.txt", "bonds490000.txt", "bonds640000.txt", "bonds810000.txt",  
def log_binomial_coefficients(M):
    """Precompute log binomial coefficients."""
    log_coeffs = gammaln(M + 1) - (gammaln(np.arange(M) + 1) + gammaln(M - np.arange(M) + 1))
    return log_coeffs


def compute_convolution_log(M, Q_n, num_q=10000):
    '''Convolution of Q_n with binomial distribution'''
    print("start convolution")
    q_values = np.linspace(0.001, 0.999, num_q)  # Probability range
    log_binom_coeffs = log_binomial_coefficients(M)  # Precompute log binomial coefficients
    M_array = np.arange(M)
    Q_q = np.zeros(num_q)
    log_q_values = np.log(q_values)
    log_1_minus_q_values = np.log(1 - q_values)
    print("start loop")
    for i in tqdm(range(num_q)):
        log_B_Mn = log_binom_coeffs + log_q_values[i] * M_array + log_1_minus_q_values[i] * (M - M_array)
        B_Mn = np.exp(log_B_Mn)  # Convert back from log space
        Q_q[i] = np.sum(B_Mn * Q_n)
    
    return q_values, Q_q



for txt in txt_list:
    prefix = f"./assignment1/data/{lattice}/"
    steps, p_inf, p_inf2, suspectability, avg_cluster_size = load_results_percolation(prefix+txt+"_results.npz")

    q_values, Q_p_inf = compute_convolution_log(len(p_inf), p_inf)
    q_values, Q_suspectability = compute_convolution_log(len(suspectability), suspectability)
    q_values, Q_avg_cluster_size = compute_convolution_log(len(avg_cluster_size), avg_cluster_size)
    q_values, Q_p_inf2 = compute_convolution_log(len(p_inf2), p_inf2)

    np.savez(f"./assignment1/convolution_results/{lattice}/{txt}_convolution.npz", 
             q_values=q_values, 
             Q_p_inf=Q_p_inf, 
             Q_suspectability=Q_suspectability, 
             Q_avg_cluster_size=Q_avg_cluster_size, 
             Q_p_inf2=Q_p_inf2)
    print(f"{txt} done and saved")

