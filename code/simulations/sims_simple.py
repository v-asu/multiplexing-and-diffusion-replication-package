################# Libraries ###############################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sio
import argparse
from pathlib import Path

################# Arguments ###############################

## we need to pass the village number as an argument
parser = argparse.ArgumentParser(description='Run simulations for gossip villages')
parser.add_argument('-vill', '--villageid', type=int)
args = parser.parse_args()

## village number
vlg = args.villageid

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
SIMS_DIR = PACKAGE_ROOT / "data" / "raw" / "sims"
SIMS_DIR.mkdir(parents=True, exist_ok=True)

################## Functions ###############################

## we need a measure for multiplexing
def mpex_intensity(G):
    """
    Calculate the intensity of multiplexing in a graph G
    """
 
    intersection = np.sum(G==2, axis=1) ## number of links in both
    union = np.sum(G>0, axis=1) ## total number of unique links

    return np.nanmean(intersection / union)

## out degree
def average_out_degree(adj_matrix):
    out_degrees = np.sum(adj_matrix, axis=1)
    return np.mean(out_degrees)

## sort networks based on out degree
def sort_networks(adj_matrices):
    avg_degrees = [(i, average_out_degree(adj_matrix)) for i, adj_matrix in enumerate(adj_matrices)]
    sorted_indices = [x[0] for x in sorted(avg_degrees, key=lambda y: y[1], reverse=True)]
    return [adj_matrices[i] for i in sorted_indices]

## remove links from G_h to match the out degree of G_l
def link_pruning(G_h, G_l):
    """
    This will prune links from G_h
    """
    G_hc = G_h.copy() ## make a copy
    avg_degree_high = average_out_degree(G_hc)
    avg_degree_low = average_out_degree(G_l)

    while avg_degree_high > avg_degree_low:
        # Find a random edge to remove
        row, col = np.nonzero(G_hc)
        index = np.random.choice(len(row))
        edge_to_remove = (row[index], col[index])
        
        # Remove the edge
        G_hc[edge_to_remove] = 0
        
        avg_degree_high = average_out_degree(G_hc)

    return G_hc

## output two multiplexed networks
def multiplexed_networks(adj_matrices):
    """
    adj_matrices: adjacency matrices of the three networks
    1. Sort the networks in descending order of average out-degree
    2. Prune links in the second network to match the average out-degree of the third network
    3. Return the two multiplexed networks
    """

    G_sort = sort_networks(adj_matrices[1:]) ## sort the networks
    G_h_pruned = link_pruning(G_sort[0], G_sort[1]) ## prune the links

    return adj_matrices[0]+G_h_pruned, adj_matrices[0]+G_sort[1]


## Single run of simulation
def single_iteration(n, G, q, tau, delta, infected):
    """
    n = number of nodes
    G = adjacency matrix
    q = probability of transmission
    tau = threshold
    delta = probability of recovery
    infected = vector of infected nodes at time pd t-1
    """

    # Step 1: Message propagation
    dosage = np.zeros(n)
    susceptible_indices = np.where(infected == 0.0)[0]

    ## can vectorize this further if needed
    for i in susceptible_indices:
        dosage[i] += np.sum(np.random.rand(np.sum(G[i] * infected).astype(int)) < q)

    # Step 2: Recovery for nodes infected up to previous step
    infected = (1 - (np.random.rand(n) < delta)*1)*infected

    # Step 3: Infection for new nodes
    infected = np.where(dosage>= tau, 1.0, infected)

    return infected

## diffusion simulation (full run)
def diffusion_simulation(G, q, tau, delta, seed_set, max_iter = 1000):
    """
    G = adjacency matrix
    q = probability of transmission
    tau = threshold
    delta = probability of recovery
    seed_set = initial set of infected nodes
    max_iter = maximum number of iterations
    """

    n = G.shape[0]
    infected = np.zeros(n)

    ## start the infection from a random seed set
    infected[seed_set] = 1.0

    ## run the simulation
    for step in range(max_iter):

        # get share of infected nodes in the system
        share_infected_old = np.mean(infected)

        # run the single iteration
        infected = single_iteration(n, G, q, tau, delta, infected)

        # get share of infected nodes
        share_infected_new = np.mean(infected)

        # check for convergence
        if np.abs(share_infected_new - share_infected_old) < 1e-8:
            break

    ## run the simulation for 100 times and take the average (for removing more noise)
    infected_list = np.zeros(100)
    for jj in range(100):
        infected = single_iteration(n, G, q, tau, delta, infected)
        infected_list[jj] = np.mean(infected)

    return np.mean(infected_list)

## Run multiple times for a single village for a given q and delta 
def diffusion_simulation_mpex(adj_list, q, tau, delta, max_iter = 1000, num_simulations = 10):

    ## generate empty matrix to store results (q,delta,mpex1,res1,mpex2,res2)
    res_out = np.zeros((num_simulations, 6))

    for i in range(num_simulations):

        M1, M2 = multiplexed_networks(adj_list) ## do pruning and return multiplexed networks

        n_g = M1.shape[0]
        seed_set = np.random.choice(np.arange(n_g), size = int(np.floor(np.sqrt(n_g))), replace=False)

        result1 = diffusion_simulation(M1, q, tau, delta, seed_set, max_iter)
        result2 = diffusion_simulation(M2, q, tau, delta, seed_set, max_iter)

        mpex1 = mpex_intensity(M1)
        mpex2 = mpex_intensity(M2)

        res_out[i] = np.array([q, delta, mpex1, result1, mpex2, result2])

    return res_out

## Run for a list of q and delta values
def run_multiple_simulations(adj_list, q_values, tau, delta_values, max_iter=1000, num_simulations = 10):
    num_q = len(q_values)
    num_delta = len(delta_values)
    total_simulations = num_q * num_delta

    # Initialize result matrix
    result_matrix = np.zeros((total_simulations * num_simulations, 6))

    # Iterate over each combination of q and delta
    for i, q in enumerate(q_values):
        for j, delta in enumerate(delta_values):
            index = i * num_delta + j
            # Run diffusion_simulation_mpex for each combination
            results = diffusion_simulation_mpex(adj_list, q, tau, delta, max_iter, num_simulations)
            result_matrix[(index * num_simulations):(index + 1) * num_simulations, :] = results

    return result_matrix


################# Data ###############################

## RCT village adjmats (directed)
adj_mat_rfe = sio.loadmat(PACKAGE_ROOT / "data" / "raw" / "rct_villages" / "rct_network_adjacency_layers.mat")
adj_mat_rfe = adj_mat_rfe["Z"] # 1 X 71 matrix (we have 70 villages)

## Directed Networks
G_rfe = {}
G_rfe["keroricecome"] = [adj_mat_rfe[0,i][0,1] for i in np.r_[0:26, 27:71]]
G_rfe["visitcome"] = [adj_mat_rfe[0,i][0,2] for i in np.r_[0:26, 27:71]]
G_rfe["keroricego"] = [adj_mat_rfe[0,i][0,0] for i in np.r_[0:26, 27:71]]
G_rfe["visitgo"] = [adj_mat_rfe[0,i][0,3] for i in np.r_[0:26, 27:71]]
G_rfe["decision"] = [adj_mat_rfe[0,i][0,4] for i in np.r_[0:26, 27:71]]
G_rfe["advise"] = [adj_mat_rfe[0,i][0,5] for i in np.r_[0:26, 27:71]]

focus_layers = ["keroricego", "visitgo", "advise"] ## we need three layers for the simulation

G_list = [[G_rfe[focus_layers[0]][i].copy(), G_rfe[focus_layers[1]][i].copy(),
            G_rfe[focus_layers[2]][i].copy()] for i in range(70)]

################### Parameters ###############################

q_list = [0.1, 0.2, 0.3, 0.4, 0.5]
delta_list = [0.1, 0.2, 0.3, 0.4, 0.5]
tau = 1.0
m_iter = 1000
n_sims = 500
 
out_simple = run_multiple_simulations(G_list[vlg], q_list, tau, delta_list, m_iter, n_sims)

## save the results
np.save(SIMS_DIR / f"sim_simple_contagion_village_{vlg}.npy", out_simple)
