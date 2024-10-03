# Julia version of RHOSTS for computing the Homological Scaffold and Violating Triangles from Multivariate Time Series

### Prerequisites:
1. **Julia Installation**: Download and install Julia version 1.10.2 (tested). Update line 4 in the `compile_package.sh` script with your Julia installation path, for example: 
    ```sh
    julia_path="/home/$user/julia-1.10.2/bin/julia"
    ``` 
2. **Compile the Package**: Run the following command to compile the package (compilation may take 3-10 minutes):
    ```sh
    bash compile_package.sh
    ```

### Running the Code:
1. **Launch the Code**: Use one of the provided bash scripts, such as:
    ```sh
    example_launch_julia_code_parallel.sh
    ```
    This script generates output in the same format as the original Python-based RHOSTS code. Hypercomplexity values are computed using a custom Julia adaptation of the Sliced Wasserstein distance from the `persim` Python library (see Note below).

2. **Output Files**: 
   - **Homological Scaffold**: Stored in an `.hd5` file. The format is a dictionary where:
      - Key: Time (Note: Julia starts indexing at 1, unlike Python which starts at 0).
      - Value: Weighted adjacency matrix of the scaffold (default is the frequency scaffold).
   - **Violating Triangles**: Saved in a separate `.hd5` file. The format is a dictionary where:
      - Key: Time (starting from 1).
      - Value: List of violating triangles, where zeros indicate non-violating triangles.
   - **Higher-order Indicators**: Print on the standard output the higher-order indicators for each time point.

3. **Help Options**: Use the `-h` flag for available launch options:
    ```sh
    ./example_launch_julia_code_parallel.sh -h
    ```

### Notes on Wasserstein Distance:
- **Hypercomplexity Computation**: By default, the hypercomplexity and its corresponding contributions are computed using a Julia conversion of the Sliced Wasserstein distance between persistence diagrams. This is adapted from the Python library `persim`, specifically [sliced_wasserstein.py](https://github.com/scikit-tda/persim/blob/master/persim/sliced_wasserstein.py), as used in the Python version of RHOSTS 

- **Wasserstein Metric Options**: The Julia implementation is significantly faster than the Python equivalent. To enable the true Wasserstein metric, use the `-w` option. This ensures higher accuracy without the performance limitations of Python.

- **Distance Calculation Differences**: Notice that the Python version of Wasserstein distance reports the shortest distance from a point to the diagonal of the persistence diagram, while the Julia implementation computes the hypotenuse of the triangle formed by the point and the diagonal. This means that in Julia, the distance is scaled by a factor of $\sqrt{2}$, i.e., $\bar{D} = D / \sqrt{2}$.

