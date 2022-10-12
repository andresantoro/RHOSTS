import sys
import random
import networkx as nx
import os
import numpy as np


def logistic(r, x):
    return (1 - r * (x**2))
    # alternative  chaotic map
    # return(4*x*(1-x))


# Reshuffle with probability p the links of the graph G yet ensuring that the graph
# remains connected (i.e. the degree distribution remains the same)
def double_switch_connected(G, p):
    N = len(G.nodes())
    nswap = int(p * N)
    counter_swaps = 0
    while counter_swaps < nswap:
        edges_to_switch = random.sample(list(G.edges()), 2)
        nodes_distinct = set(np.ravel(edges_to_switch))
        total_edges = G.edges()
        flag = 0
        if len(nodes_distinct) == 4:
            new_edge1 = (edges_to_switch[0][0], edges_to_switch[1][0])
            new_edge2 = (edges_to_switch[0][1], edges_to_switch[1][1])
            if new_edge1 not in total_edges and new_edge2 not in total_edges:
                G.add_edges_from([new_edge1, new_edge2])
                G.remove_edges_from(edges_to_switch)
                # Undo the swapping if the graph is not connected anymore
                if len(list(nx.connected_components(G))[0]) != N:
                    G.add_edges_from(edges_to_switch)
                    G.remove_edges_from([new_edge1, new_edge2])
                    flag = 1
                else:
                    counter_swaps += 1
                    flag = 2

            new_edge1 = (edges_to_switch[0][0], edges_to_switch[1][1])
            new_edge2 = (edges_to_switch[0][1], edges_to_switch[1][0])
            if flag <= 1:
                if new_edge1 not in total_edges and new_edge2 not in total_edges:
                    G.add_edges_from([new_edge1, new_edge2])
                    G.remove_edges_from(edges_to_switch)
                    # Undo the swapping if the graph is not connected anymore
                    if len(list(nx.connected_components(G))[0]) != N:
                        G.add_edges_from(edges_to_switch)
                        G.remove_edges_from([new_edge1, new_edge2])
                    else:
                        counter_swaps += 1
    return(G)


def generate_graph(N, m, p):
    G = nx.watts_strogatz_graph(n=N, k=m, p=p)
    # G=double_switch_connected(G,p=p)
    return(G)


# Generate couple map lattice according to this equation:
# x_i^t= (1-\eps)f[x_i^{t-1}] + \eps/(# neigbours) \sum_{j in \neighbours} f[x_j]^{t-1}
def generate_couple_map_fromG(T, N, r, epsilon, transient_time, G):
    series = {}

    # Filing the dictionary with N initial random values
    for index_series in range(0, N):
        s = random.random()
        series[index_series] = [s]

    # Generate the coupled maps on the graph G for a length of size T
    # (we discard the first transient_time elements to remove the transient)
    for t in range(1, T + transient_time + 1):
        for node_i in list(G.nodes()):
            prev_state = (1 - epsilon) * logistic(r, series[node_i][t - 1])

            list_neighbours_j = list(G.neighbors(node_i))
            N_j = len(list_neighbours_j)

            # Loop over the neighbours of node i
            neigh_state = 0
            for node_j in list_neighbours_j:
                neigh_state += logistic(r, series[node_j][t - 1])
            neigh_state *= ((1 / N_j) * epsilon)

            new_point = prev_state + neigh_state
            series[node_i].append(new_point)
    return(series)


def print_map(couple_map_dictionary, T, N, counter, epsilon_current, transient_time):
    for index_series in range(0, N):
        print("0", end=' ')
    print(epsilon_current)

    for t in range(transient_time + 1, T + transient_time + 1):
        print(t - transient_time + T * counter, end=' ')
        for index_series in range(0, N):
            print(couple_map_dictionary[index_series][t], end=' ')
        print()


if len(sys.argv) <= 5:
    print("Error! {0} <Ntimesteps> <Nnodes> <r> <WS prob. value> <epsilon_values> ".format(
        sys.argv[0]))
    exit(1)

#####################################################
T = int(sys.argv[1])  # Length of the time points
N = int(sys.argv[2])  # Number of nodes in the graph
r = float(sys.argv[3])  # Order parameter for the logistic map
p = float(sys.argv[4])
m = 2  # Number of neighbours in the WS network
epsilon_list = []
transient_time = 100000


for i in range(5, len(sys.argv)):
    epsilon_current = float(sys.argv[i])  # Epsilon value
    epsilon_list.append(epsilon_current)
sys.stderr.write("Epsilon values: %s\n" % str(epsilon_list))

# Generate a WS grapp
G = generate_graph(N, m, p)
# random.seed(10)

# Generating a multivariate time series by concatenating different regimes (if different value of eps are inserted)
# Each regime is separated by a row full of zeroes
for index_cycle, epsilon_current in enumerate(epsilon_list):
    coupled_map = generate_couple_map_fromG(
        T, N, r, epsilon_current, transient_time, G)
    print_map(coupled_map, T, N, index_cycle, epsilon_current, transient_time)
