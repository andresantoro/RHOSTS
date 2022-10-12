from utils import *
from multiprocessing import Pool
import h5py
import os


def create_simplicial_framework_from_data(path_file, null_model_flag, folder_javaplex, scaffold_outdir):
    global ts_simplicial
    # Loading the data from file
    data = load_data(path_file)
    # Create the ets and the triplets_ts
    ts_simplicial = simplicial_complex_mvts(
        data, null_model_flag, folder_javaplex, scaffold_outdir)
    # return(ts_simplicial)


# This function allows to save on .hd5 file the list of violating triangles when projected at the level of edges.
# Moreover, it saves on the standard output several global quantities (line 32):
# Time; Hyper complexity indic.; Hyper complexity FC; Hyper complexity CT;
# Hyper complexity FD; Hyper coherence; Average edge violation
def handle_output(result):
    global flag_edgeweight_fn
    if flag_edgeweight_fn != None:
        f2 = h5py.File('{0}.hd5'.format(flag_edgeweight_fn), 'a')
        current_time = int(result[0])
        n_lines = len(list(result[-1].items()))
        c_values = np.array(
            list(result[-1].items())).flatten().reshape(n_lines, 4)
        m, n = np.shape(c_values)
        dset1 = f2.create_dataset(
            "{0}".format(current_time), (m, n), dtype='f', data=c_values)
        f2.close()
    print(" ".join([str(el) for el in result[:-1]]))


def launch_code_one_t(t):
    # Computing the simplicial filtration for the time t
    list_simplices_positive, list_violation_fully_coherence, hyper_coherence, list_filtration_scaffold = ts_simplicial.create_simplicial_complex(
        t)
    # Computing the persistence diagram using cechmate
    dgms1 = compute_persistence_diagram_cechmate(list_simplices_positive)
    # Maximum value that will be used to replace the inf term (important for the WS distance)
    max_filtration_weight = ts_simplicial.find_max_weight(t)
    # Replace the inf value of the persistence diagram with maximum weight
    dgms1_clean = clean_persistence_diagram_cechmate(
        dgms1, max_filtration_weight)
    # dict_file = 'PD_{0}.pck'.format(t)
    # pk.dump(dgms1_clean,open(dict_file,'wb'))

    # If flag is activated, compute the scaffold and save the list of generators on file
    if ts_simplicial.javaplex_path != False:
        compute_scaffold(list_filtration_scaffold, dimension=1, directory=ts_simplicial.scaffold_outdir,
                         tag_name_output='_{0}'.format(t),
                         javaplex_path=ts_simplicial.javaplex_path, save_generators=True, verbose=False)

    hyper_complexity = persim.sliced_wasserstein(dgms1_clean, np.array([]))

    # Since the signs of the persistence diagram are flipped,
    # then Fully Coherent contributes identify points with birth and death <=0
    dgms1_complexity_FC = dgms1_clean[(
        dgms1_clean[:, 0] < 0) & (dgms1_clean[:, 1] <= 0)]
    # Coherent Transition contributes identify points with birth < 0 and death > 0
    dgms1_complexity_CT = dgms1_clean[(
        dgms1_clean[:, 0] < 0) & (dgms1_clean[:, 1] > 0)]
    # Fully Decoherence contributes identify points with birth > 0 and death > 0
    dgms1_complexity_FD = dgms1_clean[(
        dgms1_clean[:, 0] > 0) & (dgms1_clean[:, 1] > 0)]

    # Computing the Wasserstein distances
    complexity_FC = persim.sliced_wasserstein(
        dgms1_complexity_FC, np.array([]))
    complexity_CT = persim.sliced_wasserstein(
        dgms1_complexity_CT, np.array([]))
    complexity_FD = persim.sliced_wasserstein(
        dgms1_complexity_FD, np.array([]))

    flag_violations_list = np.array(
        list_violation_fully_coherence, dtype="object")[:, 2]
    # Average edge violation
    avg_edge_violation = np.mean(flag_violations_list)

    n_ROI = ts_simplicial.num_ROI
    # From the magnitude of the list of violating triangles $\Delta_v$,
    # we compute the downward projection at the level of edges
    edge_weights = compute_edgeweight(list_violation_fully_coherence, n_ROI)

    # Report the results in a vector and print everything (except the downward projections)
    # on output
    results = [t, hyper_complexity, complexity_FC, complexity_CT,
               complexity_FD, hyper_coherence, avg_edge_violation, edge_weights]

    return(results)


############# MAIN CODE #############

if len(sys.argv) <= 1:
    print(
        "******************************************************************************\n"
        "**                                                                          **\n"
        "**              Computation of the higher-order indicators                  **\n"
        "**               starting from a multivariate time series                   **\n"
        "**                                                                          **\n"
        "**                                                                          **\n"
        "**  <filename_multivariate_series> file containing the multiv. time series  **\n"
        "**                         Format currently accepted:                       **\n"
        "**        .txt:  where columns represents the independent time series       **\n"
        "**        .mat:  where rows are ROI, and columns are the time instants      **\n"
        "**                                                                          **\n"
        "**                                                                          **\n"
        "**                     ----   Optional Variables  ----                      **\n"
        "**                                                                          **\n"
        "**    <-t t0 T> restricts the output of the higher-order indicators         **\n"
        "**                  only for the time interval [t0,T]                       **\n"
        "**                                                                          **\n"
        "**   <-p #core> represents the number of cores used for the computation of  **\n"
        "**                     the higher-order indicators                          **\n"
        "**                                                                          **\n"
        "**     <-n > computes the higher-order indicators for the null model        **\n"
        "**           constructed by independently reshuffling each signal           **\n"
        "**                                                                          **\n"
        "**   <-s <filename>> saves on filename.hdf5 the weighted network obtained   **\n"
        "**    when projecting the magnitude of the list of violations on a graph    **\n"
        "**                                                                          **\n"
        "**    <-j -path_javaplex -outdir> launch the jython code for computing the  **\n"
        "**    homological scaffold and save it in the in the folder 'outdir', it    **\n"
        "**    relies on javaplex and requires a lot of RAM for this computation     **\n"
        "**                                                                          **\n"
        "**      OUTPUT: by default the algorithm returns the following info:        **\n"
        "** Time; Hyper complexity indic.; Hyper complexity FC; Hyper complexity CT; **\n"
        "**        Hyper complexity FD; Hyper coherence; Average edge violation      **\n"
        "**                                                                          **\n"
        "******************************************************************************\n"
        "Usage: %s <filename_multivariate_series>   [-t t0 T] [-p #core] [-n] [-s <filename>] [-j <path_javaplex> <name_outdir>]\n\n" % sys.argv[0]);
    exit(1)


if __name__ == "__main__":

    # Parsing the input (still to do with the argparse library)
    [path_file, t_init, t_end, t_total, ncores, null_model_flag,
        flag_edgeweight, flag_edgeweight_fn, folder_javaplex, scaffold_outdir] = parse_input(sys.argv)

    # Empty existing file
    if flag_edgeweight_fn != None:
        f1 = h5py.File("{0}.hd5".format(flag_edgeweight_fn), "w")
        f1.close()

    # Creating the structure containing the edge and triplet signals within the Pool process
    # with this syntax, it should create problems in OS systems
    pool = Pool(processes=ncores, initializer=create_simplicial_framework_from_data,
                initargs=(path_file, null_model_flag, folder_javaplex, scaffold_outdir))

    if t_init == 0 and t_end == 0:  # By default, the script does the analysis on all the time points
        t_total = [t for t in range(t_init, ts_simplicial.T)]

    # Parallel computation
    for i in t_total:
        pool.apply_async(launch_code_one_t, (i, ), callback=handle_output)
    pool.close()
    pool.join()
