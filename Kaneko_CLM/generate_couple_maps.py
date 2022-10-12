import sys
import random
import os
import numpy as np


def logistic(r, x):
    return (1 - r * (x**2))
    # alternative  chaotic map
    # return(4*x*(1-x))


# See book of Kaneko for parameters of eps and/or r to have different regimes
# or https://en.wikipedia.org/wiki/Coupled_map_lattice for notable regimes


# Generate couple map lattice according to this equation: x_i^t= (1-\eps)f[x_i^{t-1}] + \eps/order \sum_{j in \neighbours} f[x_j]^{t-1}
def generate_couple_map(T, N, epsilon, transient_time, r, order=2):
    series = {}

    # Filing the dictionary with N initial random values
    for index_series in range(0, N):
        s = random.random()
        series[index_series] = [s]

    # Generate the coupled maps for a length of size T (yet, we discard the first transient_time elements to remove the transient)
    for i in range(1, T + transient_time + 1):
        for index_series in range(0, N):
            order_k_term = compute_neighbours(
                N, series, epsilon, index_series, i - 1, r, order)
            new_point = (1 - epsilon) * logistic(r,
                                                 series[index_series][i - 1]) + order_k_term
            series[index_series].append(new_point)
    return(series)


def compute_neighbours(N, series, epsilon, index_series, i, r, order=2):
    eps_overN = epsilon * (1 / order)
    term_left_right = int(order / 2)

    term = 0
    # Sum over the left neighbors with periodic boundary conditions
    for s in range(1, term_left_right + 1):
        term += logistic(r, series[(index_series - s) % N][i])

    # Sum over the right neighbors with periodic boundary conditions
    for s in range(1, term_left_right + 1):
        term += logistic(r, series[(index_series + s) % N][i])

    # if order is odd, then take the neighbors in an asymmetric way, int(order/2) on the left, int(order/2)+1 on the right
    if order % 2 == 1:
        s = term_left_right + 1
        term += logistic(r, series[(index_series + s) % N][i])

    return(term * eps_overN)


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
    print("Error! {0} <T> <nodes> <order neighbors> <r> <epsilon_values> (e.g. 50 10 2 1.75 0.6) ".format(
        sys.argv[0]))
    exit(1)

#####################################################
T = int(sys.argv[1])  # Number of time points
N = int(sys.argv[2])  # Number of nodes in the ring network
# Dependency from the number of neighbours in the ring lattice for the diffusive dynamics
k_order = int(sys.argv[3])
r = float(sys.argv[4])  # Order parameter for the logistic map
epsilon_list = []
transient_time = 100000
# random.seed(10)

for i in range(5, len(sys.argv)):
    epsilon_current = float(sys.argv[i])  # Epsilon value
    epsilon_list.append(epsilon_current)
sys.stderr.write("Epsilon values: %s\n" % str(epsilon_list))


# Generating a multivariate time series by concatenating different regimes (if different value of eps are inserted)
# Each regime is separated by a row full of zeroes
for index_cycle, epsilon_current in enumerate(epsilon_list):
    coupled_map = generate_couple_map(
        T, N, epsilon_current, transient_time, r, k_order)
    print_map(coupled_map, T, N, index_cycle, epsilon_current, transient_time)
