# Julia Code for computing the homological scaffold and list of violating triangles from a multivariate time series
1. After downloading julia 1.10.2 (tested) and installing it, the line 4 of the file `compile_package.sh` should be modified accordingly to include the julia path, e.g. julia_path="/home/$user/julia-1.10.2/bin/julia". The package can then be compiled using the command `bash compile_package.sh` (wait around 3-10m for the compilation)
2. To launch the code, you can use one of the existing bash scripts, e.g. `example_launch_julia_code_parallel.sh` that will generate an output having the same format as the original python code of RHOSTS (the values of the hypercomplexity will be slightly different from the python code, whilethe last two columns will be identical. See NB below for an explanation).
3. A file .hd5 will contain the information about the homological scaffold (_default is the frequency scaffold_). The format is a dictionary, where the key corresponds to time (N.B. julia starts from 1 rather than from 0), whereas the value encodes the weighted adjacency matrix of the scaffold.
4. A file .hd5 will contain the information about the list of violating triangles. The format is a dictionary, where the key corresponds to time (N.B. julia starts from 1 rather than from 0), whereas the values are encoded in a the list of size (N choose 3) violating triangles, where zeros correspond to non-violating triangles.

If you use `-h`, it will show all the different options for launching the code.



## N.B. The Hypercomplexity and the corresponding contributes are computed considering an in-house conversion of the sliced Wasserstein distance between the persistence diagrams. The results are slightly different from the python code in the RHOSTS code for the following reasons:
- The sliced Wasserstein distance function (defult option) is a Julia conversion of the python code https://github.com/scikit-tda/persim/blob/master/persim/sliced_wasserstein.py 

- To speed-up the computation, the python code relies on the Sliced Wasserstein distance rather than the proper Wasserstein distance (The sliced Wasserstein is an approximation of the "true" Wasserstain metric). The julia code has the option to rely on the true Wasserstein metric, which is similar in terms of computational time. Use the option `-w` to rely on the true Wasserstein metric.

- The Wasserstein in python provides the shortest distance with respect to the empty persistence diagram (i.e. diagonal), which for a selected point in the diagram the distance *D* corresponds the segment perpendicular to the diagonal. The implementation in Julia provides instead the hypotenuse of the triangle, i.e. the Julia code provides $\bar{D}= D / \sqrt{2}$.
