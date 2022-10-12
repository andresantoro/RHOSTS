# Importing the libraries
import numpy as np
import sys
import scipy.io as sio  # For reading the matlab .mat format
from scipy.stats import zscore, entropy
from scipy.special import binom as binomial
import collections
import pickle as pk
import itertools
import persim
import cechmate as cm

# Libraries for the scaffold (piping the filtration file to a jython code)
import os
from collections import OrderedDict
import subprocess
from subprocess import Popen, PIPE, STDOUT


def parse_input(input):
    t_init = t_end = 0
    ncores = 1
    null_model_flag = False
    flag_edgeweight = False
    flag_edgeweight_fn = None
    n_input = len(input)
    # This is the filename containing the multivariate time series
    path_file = input[1]
    javaplex_path = False
    scaffold_path = False

    for s in range(n_input):
        if sys.argv[s] == '-t' or input[s] == '-T':
            t_init = int(input[s + 1])
            t_end = int(input[s + 2])
        if sys.argv[s] == '-p' or input[s] == '-P':
            ncores = int(input[s + 1])
        if sys.argv[s] == '-n' or input[s] == '-N':
            # ->  z-score of edges and triplets is computed from the shuffled original data
            null_model_flag = True
        if sys.argv[s] == '-s' or input[s] == '-S':
            # ->  save the weighted network for each time point on the hd5 file
            flag_edgeweight = True
            flag_edgeweight_fn = input[s + 1]
        if sys.argv[s] == '-j':
            # -> Load the javaplex path and the scaffold path
            javaplex_path = input[s + 1]
            scaffold_path = input[s + 2]

    t_total = [t for t in range(t_init, t_end)]

    return(path_file, t_init, t_end, t_total, ncores, null_model_flag,
           flag_edgeweight, flag_edgeweight_fn, javaplex_path, scaffold_path)


def load_data(path_single_file):
    extension_file = path_single_file.split('.')[-1]
    if extension_file == 'mat':
        data = load_data_mat(path_single_file)
    elif extension_file == 'txt_kaneko':
        data = load_data_synthetic_kaneko(path_single_file)
    elif extension_file == 'txt':
        data = load_normaltxt(path_single_file)
    # print(np.shape(data))
    return(data)


# Load brain data in .mat format (rows are ROI, columns are the time instants)
def load_data_mat(path_single_file):
    file_to_open = path_single_file
    data = sio.loadmat(file_to_open)
    key_data = list(data.keys())[-1]
    data = data[key_data]
    return(data)

# Load the synthetic data generated from the Kaneko maps (coupled map lattices)


def load_data_synthetic_kaneko(path_single_file):
    file_to_open = path_single_file
    data = np.loadtxt(file_to_open)
    data_cleaned = []
    eps_list = []
    eps_list.append([0, data[0][-1]])
    counter = 0
    for i in data[1:]:
        # If the data are in the shape of t 0 0 0... 0 eps then save the value of eps
        if i[0] == 0 and i[1] == 0:
            N = counter
            counter = 0
            c = len(eps_list)
            eps_list.append([c * N, i[-1]])
        else:
            counter += 1
            data_cleaned.append(i)
    M = len(data[0]) - 1
    return(np.transpose(np.array(data_cleaned)[:, 1:]))

# Load synthetic data (format: columns represents independent time series )


def load_normaltxt(path_single_file):
    file_to_open = path_single_file
    data = np.loadtxt(file_to_open)
    return(np.transpose(data))


class simplicial_complex_mvts():
    def __init__(self, multivariate_time_series, null_model_flag, folder_javaplex, scaffold_outdir):
        nR, T = np.shape(multivariate_time_series)

        # Variables
        self.raw_data = multivariate_time_series
        self.num_ROI = nR
        self.T = T

        # Edges
        self.ets_indexes = {}
        self.ets_zscore = []
        self.ets_max = None

        # Triplets
        self.triplets_indexes = {}
        self.triplets_ts_zscore = []
        self.triplets_max = None

        # Variables for the filtration
        self.list_simplices = []
        self.list_violations = []
        self.percentage_of_triangles_discarded = 0
        self.percentage_CC_triangles_positive = 0
        self.percentage_CC_triangles_negative = 0

        # Variables for the scaffold
        self.javaplex_path = folder_javaplex
        self.scaffold_outdir = scaffold_outdir

        # If null model is on, do an independent reshuffling of the original time series
        if null_model_flag == True:
            self.shuffle_original_data()

        # Computing the z-score of the initial data and replace the variable self.raw_data
        self.compute_zscore_data()

        # Initialising the variables by computing the edges and triplets
        self.compute_edges_triplets()

    def shuffle_original_data(self):
        # Shuffling the original
        data = np.array([list(np.random.permutation(row))
                         for row in self.raw_data])
        # Save it
        self.raw_data = data

    def compute_zscore_data(self):
        # Computing the z-score of the data
        self.raw_data = zscore(self.raw_data, axis=1)

    # Initial setup: computation of the edges and triplets

    def compute_edges_triplets(self):
        #-------------------------EDGES-----------------------------
        # Number of edges
        N_edges = int(binomial(self.num_ROI, 2))
        # Indices for the products
        u, v = np.triu_indices(self.num_ROI, k=1, m=self.num_ROI)
        self.ets_zscore = np.zeros((N_edges, 2))
        self.ets_max = np.zeros((self.T))
        l_index_prev = 0
        l_index_next = self.num_ROI - 1
        gap = l_index_next - l_index_prev

        # To save memory, ets_zscore will save the mean and std of each independent time series
        # instead of all the z-scored edges
        # ets_max -> is a vector 1xT containing the maximum between all the z-scored edges
        # To decrease the RAM usage, the product is done in batches of size < N
        for i in range(self.num_ROI):
            c_prod = self.raw_data[u[l_index_prev:l_index_next]
                                   ] * self.raw_data[v[l_index_prev:l_index_next]]
            self.ets_zscore[l_index_prev:l_index_next] = np.array(
                [np.mean(c_prod, axis=1), np.std(c_prod, axis=1)]).T
            self.ets_max = np.max(
                np.vstack((self.ets_max, np.abs((c_prod - np.tile(self.ets_zscore[l_index_prev:l_index_next, 0], (self.T, 1)).T) /
                                                np.tile(self.ets_zscore[l_index_prev:l_index_next, 1], (self.T, 1)).T))), axis=0)
            l_index_prev = l_index_next
            l_index_next += (self.num_ROI - i - 2)

        # Save in a dictionary the indexes of the edges (i,j), i.e. key: index (from 0 to N*(N-1)/2), value: (i,j)
        self.ets_indexes = dict(zip(np.arange(N_edges), zip(u, v)))

        #------------------------TRIPLETS----------------------------

        # Number of triplets
        N_triplets = int(binomial(self.num_ROI, 3))
        # Indices for the products
        self.idx_list_triplets = list(
            itertools.combinations(range(self.num_ROI), r=3))
        indices = np.array(self.idx_list_triplets)
        u, v, w = indices[:, 0], indices[:, 1], indices[:, 2]
        # Filling the rows of tts with the product
        self.triplets_ts_zscore = np.zeros((N_triplets, 2))
        self.triplets_max = np.zeros((self.T))
        l_index_prev = 0
        l_index_next = int((self.num_ROI - 1) * (self.num_ROI - 2) / 2)

        # Same as above,
        # To save memory, triplets_ts_zscore will save the mean and std of each independent triplet time series
        # instead of all the z-scored triplets
        # triplets_max -> is a vector 1xT containing the maximum between all the z-scored triplets
        # To decrease the RAM usage, the product is done in batches of size < N*(N-1)
        for i in range(self.num_ROI):
            #
            # print(i,l_index_prev,l_index_next)
            c_prod = self.raw_data[u[l_index_prev:l_index_next]] * \
                self.raw_data[v[l_index_prev:l_index_next]] * \
                self.raw_data[w[l_index_prev:l_index_next]]
            self.triplets_ts_zscore[l_index_prev:l_index_next] = np.array(
                [np.mean(c_prod, axis=1), np.std(c_prod, axis=1)]).T
            self.triplets_max = np.max(
                np.vstack((self.triplets_max, np.abs((c_prod - np.tile(self.triplets_ts_zscore[l_index_prev:l_index_next, 0], (self.T, 1)).T) /
                                                     np.tile(self.triplets_ts_zscore[l_index_prev:l_index_next, 1], (self.T, 1)).T))), axis=0)
            l_index_prev = l_index_next
            l_index_next += int((self.num_ROI - i - 2) *
                                (self.num_ROI - i - 3) / 2)
            gap = l_index_next - l_index_prev
        # Saving the indices of all the triplets
        self.triplets_indexes = dict(zip(np.arange(N_triplets), indices))

    # Function that, for a specific time t, computes the maximum between edges and triplets
    # This is used to replace the infty term after computing the persistence diagram

    def find_max_weight(self, t):
        edges_abs_max = self.ets_max[t]
        triplets_abs_max = self.triplets_max[t]
        m = np.max([edges_abs_max, triplets_abs_max])
        return(m)

    # Function that remaps the weight of a k-order products using the pure coherence rule.

    def correction_for_coherence(self, current_list_sign, current_weight):
        # If the original signals are fully coherent, then the corresponding weight becomes positive, otherwise negative
        flag = 0
        coherence = coherence_function(current_list_sign)
        # If all the signs are concordant then set the weight sign as positive, otherwise negative
        if coherence == 1:
            weight_corrected = np.abs(current_weight)
        else:
            weight_corrected = -np.abs(current_weight)
        return(weight_corrected)

    # Function that creates the list of simplices (and provide also the list of violations)

    def create_simplicial_complex(self, t_current):

        # Creating the list of simplicial complex with all the edges and triangles
        list_simplices = []

        # Selecting the extremal weight between edges and triplets. It will be assigned to all the nodes (i.e. nodes enter at the same instant)
        m_weight = np.max([np.ceil(self.triplets_max[t_current]), np.ceil(
            self.ets_max[t_current])])
        # Adding all the nodes from the beginning with the same weights
        for i in range(self.num_ROI):
            list_simplices.append(([i], m_weight))

        # Adding the edges:
        # Also, modify the signs of the weights to correct the z-score so that: if the edge signal is fully coherent, then assign a positive sign, otherwise negative
        for i in self.ets_indexes:
            indexes_ij = self.ets_indexes[i]
            c_mean = self.ets_zscore[i][0]
            c_std = self.ets_zscore[i][1]
            weight_current = (self.raw_data[indexes_ij[0]][t_current]
                              * self.raw_data[indexes_ij[1]][t_current] - c_mean) / c_std
            list_of_signs = [self.raw_data[indexes_ij[0]]
                             [t_current], self.raw_data[indexes_ij[1]][t_current]]
            weight_current_corrected = self.correction_for_coherence(
                list_of_signs, weight_current)
            list_simplices.append((indexes_ij, weight_current_corrected))

        # Adding the triplets
        # Here I modify the signs of the weights, if it is fully coherent I assign a positive sign, otherwise negative
        for i in self.triplets_indexes:
            indexes_ijk = self.triplets_indexes[i]
            c_mean = self.triplets_ts_zscore[i][0]
            c_std = self.triplets_ts_zscore[i][1]
            weight_current = (self.raw_data[indexes_ijk[0]][t_current] * self.raw_data[indexes_ijk[1]][t_current] *
                              self.raw_data[indexes_ijk[2]][t_current] - c_mean) / c_std
            list_of_signs = [self.raw_data[indexes_ijk[0]][t_current],
                             self.raw_data[indexes_ijk[1]][t_current], self.raw_data[indexes_ijk[2]][t_current]]
            weight_current_corrected = self.correction_for_coherence(
                list_of_signs, weight_current)
            list_simplices.append((indexes_ijk, weight_current_corrected))

        list_simplices_for_filtration, list_violations, percentage_of_triangles_discarded, list_simplices_scaffold_all = self.fix_violations(
            list_simplices, t_current)
        return(list_simplices_for_filtration, list_violations, percentage_of_triangles_discarded, list_simplices_scaffold_all)

    # Function that removes all the violating triangles to create a proper filtration
    def fix_violations(self, list_simplices, t_current):
        # Sorting the simplices in a descending order according to weights
        sorted_simplices = sorted(
            list_simplices, key=lambda x: x[1], reverse=True)

        # Remove the violations
        list_violating_triangles = []
        list_simplices_for_filtration = []
        set_simplices = set()
        counter = 0
        triangles_count = 0
        violation_triangles = 0
        violation_triangles_negativeterms = 0

        CC_triangles_positive = 0
        total_CC_triangles = 0
        CC_triangles_negative = 0

        # List all the valid simplices that will be
        # used for the computation of the scaffold
        list_simplices_scaffold_all = OrderedDict()
        counter_simplices_all = 0

        # Loop over the sorted simplices, and flipping the sign of all the weights (so that the points in the persistence diagram are above the diagonal)
        for index, i in enumerate(sorted_simplices):
            simplices, weight = i

            # If the current simplex is an edge or a node, then I will immediately include it
            if len(simplices) <= 2:
                # If it's a node I add it in the list of all simplices
                if len(simplices) == 1:
                    list_simplices_scaffold_all[str(list(simplices))] = [str(
                        counter_simplices_all), str(-weight)]
                else:
                    # If the simplex is an edge, check whether the weight of the last inserted element
                    # is different from the current one, if not use the same order idx of appearance (counter_simplices_all)
                    if weight != sorted_simplices[index - 1][1]:
                        counter_simplices_all += 1
                        list_simplices_scaffold_all[str(list(simplices))] = [str(
                            counter_simplices_all), str(-weight)]
                    else:
                        list_simplices_scaffold_all[str(list(simplices))] = [str(
                            counter_simplices_all), str(-weight)]

                list_simplices_for_filtration.append((simplices, -weight))
                set_simplices.add(tuple(simplices))
                counter += 1
            else:
                # If the current simplex is a triplet, I check whether all the sub-simplices have been included.
                flag = 0
                for t in itertools.combinations(simplices, 2):
                    if t in set_simplices:
                        flag += 1

                # If all the sub-simplices already belong to the set, then I add it in the filtration
                if flag == 3:
                    set_simplices.add(tuple(simplices))
                    list_simplices_for_filtration.append((simplices, -weight))
                    counter += 1
                    if weight != sorted_simplices[index - 1][1]:
                        counter_simplices_all += 1
                        list_simplices_scaffold_all[str(list(simplices))] = [str(
                            counter_simplices_all), str(-weight)]
                    else:
                        list_simplices_scaffold_all[str(list(simplices))] = [str(
                            counter_simplices_all), str(-weight)]

                    # Count the number of positive triangles that are in the filtration
                    if weight >= 0:
                        triangles_count += 1
                else:
                    # Count the violations only for fully coherent state (--- or +++)
                    if weight >= 0:
                        violation_triangles += 1
                        list_violating_triangles.append(
                            (simplices, np.abs(weight), 3 - flag))

                    else:
                        violation_triangles_negativeterms += 1

        # Fraction of positive triangle discarderd (a.k.a. the hyper coherence)
        hyper_coherence = (1.0 * violation_triangles) / \
            (triangles_count + violation_triangles)
        list_simplices_scaffold_all = dict(list_simplices_scaffold_all)

        return(list_simplices_for_filtration, list_violating_triangles, hyper_coherence, list_simplices_scaffold_all)


# Function that check for the pure coherence rule
def coherence_function(vector):
    n = len(vector)
    temp = 0
    for el in vector:
        temp += np.sign(el)
    exponent = np.sign(n - np.abs(temp))
    res = (-1)**exponent
    return(res)


# Clean the persistence diagram, i.e. replace the inf term with a value given in input
# In our case the value corresponds to the maximum weight at a specific time t
def clean_persistence_diagram_cechmate(dgms, max_filtration, order=1):
    pdgm = []
    for i in dgms[order]:
        if i[1] == np.inf:
            pdgm.append([i[0], max_filtration])
        else:
            pdgm.append(i)
    return(np.array(pdgm))


# Compute the persistence diagram using cechmate
def compute_persistence_diagram_cechmate(list_simplices_all):
    dgms = cm.phat_diagrams(list_simplices_all, show_inf=True, verbose=False)
    return(dgms)


# Compute the edge weight projection starting from the list of violating triangles
def compute_edgeweight(list_violations, num_ROI):
    Nviolations = len(list_violations)
    # Computing the edge weights from the violating triangles
    # (it also keeps the count on the number triangles one edge belongs to)
    edge_weight = {}
    for element in list_violations:
        triplets, weight, _ = element
        for edgeID in itertools.combinations(triplets, 2):
            if edgeID in edge_weight:
                edge_weight[edgeID] = np.add(
                    edge_weight[edgeID], [weight, 1.0])
            else:
                edge_weight[edgeID] = [weight, 1.0]
    return(edge_weight)


def compute_scaffold(clique_dic_file, dimension, directory='./', tag_name_output='_0', javaplex_path='/home/andrea/javaplex/lib/',
                     save_generators=True, verbose=False, python_persistenthomologypath='persistent_homology_calculation.py'):
    '''
    Function that calls the jython code to compute the scaffold.
    It requires a valid filtration "clique_dic_file", which is piped as an input string for the jython code
    persistent_homology_calculation.py
    Notice that clique_dic_file is a dictionary, where keys represent the simplices
    in the string form '[i,j,l]' (list not tuple!).
    Values are in the form: [str(idx_arrival), str(weight)]
    Example of a valid dictionary with three nodes (with the same arrival idx and weight),
    three edges and a triangle:
    {'[0]':['0','1.5'],'[1]':['0','1.5'],'[2]':['0','1.5'],'[0,2]':['1','2.5'],'[0,1]':['2','3.2'],
     '[1, 2]':['3','3.5'], '[0,1,2]':['4','5.0']]}
    The parameter dimension represents the homology group '''

    # Check that the persistent homology python scripts exists in the current directory:
    # print(os.path.exists(python_persistenthomologypath))
    if os.path.exists(python_persistenthomologypath) == False:
        sys.stderr.write("File {0} is not present in the current directory. I cannot launch the scaffold code! Skipping...\n".format(
            python_persistenthomologypath))
    else:
        Clique_dictionary = str(clique_dic_file)
        args = ["jython", python_persistenthomologypath]
        for opt in [dimension, directory, tag_name_output, javaplex_path, save_generators]:
            args.extend([str(opt)])
        s = subprocess.Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)

        # Here I'm piping the dictionary as input string for the jython code
        grep_stdout = s.communicate(input=Clique_dictionary.encode('utf-8'))[0]
        if verbose == True:
            print(grep_stdout.decode())
