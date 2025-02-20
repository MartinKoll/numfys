import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm
from scipy.stats import linregress

prefix = "./assignment1/convolution_results/triangular/"
txt_list = ["bonds10000.txt_convolution.npz", "bonds40000.txt_convolution.npz", "bonds90000.txt_convolution.npz", 
            "bonds160000.txt_convolution.npz", "bonds250000.txt_convolution.npz", "bonds360000.txt_convolution.npz", 
            "bonds490000.txt_convolution.npz", "bonds640000.txt_convolution.npz", "bonds810000.txt_convolution.npz", "bonds1000000.txt_convolution.npz"]

grid_size = [txt.split(".")[0][5:] for txt in txt_list]
epsilon = np.array([np.sqrt(int(L)) for L in grid_size])

def load_results_convolution(file):
    results = np.load(file)
    return results["q_values"], results["Q_p_inf"], results["Q_suspectability"], results["Q_avg_cluster_size"], results["Q_p_inf2"]

def find_pc(q_values, Q_p_inf):
    log_epsilon = np.log(epsilon)
    log_giant = np.log(Q_p_inf)
    log_giant_T = log_giant.T
    
    r2_vals, q_vals = [], []
    r2_max, pc, beta_over_nu = 0, 0, 0
    
    for i, q in enumerate(q_values):
        if Q_p_inf[-1][i] < 0.11:
            continue
        lin_fit = linregress(log_epsilon, log_giant_T[i])
        r2 = lin_fit.rvalue**2
        if r2 > r2_max:
            r2_max, pc, beta_over_nu = r2, q, lin_fit.slope
        r2_vals.append(r2)
        q_vals.append(q)
    
    plt.figure()
    plt.plot(q_vals, r2_vals, marker='o')
    plt.axvline(pc, color='r', linestyle='--', label=f'pc = {pc:.4f}')
    plt.xlabel('q values')
    plt.ylabel('r^2 values')
    plt.title('Finding pc')
    plt.legend()
    plt.show()
    
    return pc, -beta_over_nu

def max_s_function(Q_avg, q_values):
    max_s = [max(avg) for avg in Q_avg]
    max_q = [q_values[list(avg).index(max(avg))] for avg in Q_avg]
    return max_s, max_q

def finite_size_scaling(q_values, Q_avg):
    max_s, max_q = max_s_function(Q_avg, q_values)
    log_epsilon = np.log(epsilon)
    log_max_s = np.log(max_s)
    lin_fit = linregress(log_epsilon, log_max_s)
    
    plt.figure()
    plt.plot(log_epsilon, log_max_s, 'o', label='Data')
    plt.plot(log_epsilon, lin_fit.intercept + lin_fit.slope * log_epsilon, '--', label=f'Fit: slope={lin_fit.slope:.4f}')
    plt.xlabel('log(epsilon)')
    plt.ylabel('log(max_s)')
    plt.title('Finite Size Scaling')
    plt.legend()
    plt.show()
    
    return lin_fit.slope

def extract_nu(q_values, max_q, p_c):
    log_epsilon = np.log(epsilon)
    log_q_diff = np.log(np.abs(max_q - p_c) + 1e-10)  # Prevent log(0)
    lin_fit = linregress(log_epsilon, log_q_diff)
    
    plt.figure()
    plt.plot(log_epsilon, log_q_diff, 'o', label='Data')
    plt.plot(log_epsilon, lin_fit.intercept + lin_fit.slope * log_epsilon, '--', label=f'Fit: slope={lin_fit.slope:.4f}')
    plt.xlabel('log(epsilon)')
    plt.ylabel('log(|q_max - pc|)')
    plt.title('Extracting Nu')
    plt.legend()
    plt.show()
    
    return -1 / lin_fit.slope

Q_pinf, Q_sus, Q_avg, Q_pinf2 = [], [], [], []
for txt in txt_list:
    q_values, Q_p_inf, Q_suspectability, Q_avg_cluster_size, Q_p_inf2 = load_results_convolution(prefix+txt)
    Q_pinf.append(Q_p_inf)
    Q_sus.append(Q_suspectability)
    Q_avg.append(Q_avg_cluster_size)
    Q_pinf2.append(Q_p_inf2)

Q_pinf, Q_sus, Q_avg, Q_pinf2 = np.array(Q_pinf), np.array(Q_sus), np.array(Q_avg), np.array(Q_pinf2)

max_s, max_q = max_s_function(Q_avg, q_values)
gamma_over_nu = finite_size_scaling(q_values, Q_avg)
pc, beta_over_nu = find_pc(q_values, Q_pinf)
nu = extract_nu(q_values, max_q, pc)

print(f"Percolation threshold (pc): {pc}")
print(f"Critical exponent β: {beta_over_nu * nu}")
print(f"Critical exponent γ: {gamma_over_nu * nu}")
print(f"Critical exponent ν: {nu}")


# Hex
# Percolation threshold (pc): 0.6524597459745974
# Critical exponent β: 0.1465213428482498
# Critical exponent γ: 2.2428157962018647
# Critical exponent ν: 1.2834426330060353

# Square
# Percolation threshold (pc): 0.5003493349334933
# Critical exponent β: 0.12666768086590285
# Critical exponent γ: 2.4936796168078765
# Critical exponent ν: 1.420942821283726

# Triangular
# Percolation threshold (pc): 0.347041204120412
# Critical exponent β: 0.1581596970108306
# Critical exponent γ: 2.1990002993909576
# Critical exponent ν: 1.2583101530012957