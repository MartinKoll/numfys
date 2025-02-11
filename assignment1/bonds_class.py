import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from functions import find_root, process_bonds
from matplotlib.colors import ListedColormap

class bonds:
    def __init__(self, bonds_txt, type):
        #text reading
        if type == "square":
            prefix = "./assignment1/bonds/square_lattice/"
        elif type == "triangular":
            prefix = "./assignment1/bonds/triangular_lattice/"
        elif type == "honeycomb":
            prefix = "./assignment1/bonds/honeycomb_lattice/"
        self.bonds_txt = bonds_txt
        self.bonds = np.loadtxt(f"{prefix}{bonds_txt}", skiprows=2, dtype=int)
        self.N, self.number_of_bonds = np.loadtxt(f"{prefix}{bonds_txt}", max_rows=2, dtype=int)
        self.L = int(np.sqrt(self.N))
        #root creation
        self.sites = np.full(self.N, -1)
        self.cluster = np.zeros((self.L, self.L), dtype=int)
        
        self.largest_cluster_index = 0
        self.largest_cluster_size = 1

        self.number_of_activated_bonds = 0
        self.probability = []
        self.steps = []
        self.avg_cluster_size = []

        self.shuffle_bonds()

    def shuffle_bonds(self, seed=None):
        np.random.seed(seed)  # Set seed for reproducibility
        np.random.shuffle(self.bonds)  # Shuffle the entire list in place


    def activate_bonds(self, number_of_activations):
        if not isinstance(number_of_activations, int) or number_of_activations <= 0:
            number_of_activations = self.number_of_bonds
        for i in tqdm(range(number_of_activations)):
            bond = self.bonds[i]
            self.add_bonds(bond)
            self.number_of_activated_bonds += 1
            if i&((self.number_of_bonds)//200)==0:
                self.steps.append(i/self.number_of_bonds)
                self.probability.append(self.largest_cluster_size/self.N)
            # if i%((self.number_of_bonds)//20) == 0:
            #     self.visualize(f"{i}")

    def find_root(self, index):
        if self.sites[index] < 0:
            return index
        self.sites[index] = self.find_root(self.sites[index])  # Path compression
        return self.sites[index]
    
    def add_bonds(self, bond):
        root1 = find_root(self.sites, bond[0])
        root2 = find_root(self.sites, bond[1])
        if root1 == root2:
            return
        if root1 <= root2:
            self.sites[root1] += self.sites[root2]
            self.sites[root2] = root1
            #self.sites[self.sites == root2] = root1
            if -self.sites[root1] > self.largest_cluster_size:
                self.largest_cluster_size = -self.sites[root1]
                self.largest_cluster_index = root1
        else:
            self.sites[root2] += self.sites[root1]
            self.sites[root1] = root2
            #self.sites[self.sites == root1] = root2
            if -self.sites[root2] > self.largest_cluster_size:
                self.largest_cluster_size = -self.sites[root2]
                self.largest_cluster_index = root2

    def visualize(self, title=""):
        #cluster_nodes = np.where(self.sites == largest_cluster_index, 1, 0).reshape(self.L, self.L)
        cluster_nodes = np.zeros((self.L, self.L))
        for i in range(self.sites.shape[0]):
            site = find_root(i)
            if site == self.largest_cluster_index:
                cluster_nodes[i//self.L, i%self.L] = 1

        plt.imshow(cluster_nodes, cmap=ListedColormap(["white", "blue"]), vmin=0, vmax=1)
        plt.title(f"Largest Cluster at {title} iterations")
        plt.savefig(f"./assignment1/pics/frame_{title}.png")
        #plt.show()

    def plot_p(self):
        fig, ax = plt.subplots()
        ax.plot(self.steps, self.probability)
        ax.set_xlabel("Number of activated bonds")
        ax.set_ylabel("P")
        ax.set_title("Percolation Probability")
        #plt.savefig("./assignment1/pics/p.png")
        plt.show()

    def save_p(self):
        with open("./assignment1/pics/p.txt", "w") as f:
            for i in range(len(self.steps)):
                f.write(f"{self.steps[i]} {self.probability[i]}\n")      


#bonds1 = bonds("bonds1000.txt")
#bonds1.shuffle_bonds()
#bonds1.activate_bonds(0)
#bonds1.visualize("end")
#bonds1.plot_p()