#!/bin/bash
### Launch the code for the times 0-5 and using 5 cores 
### and save the magnitude of the projected violating triangles \Delta_v at the level of edges
### on the file "edges_projection.hd5"


##The code need to be launched from the directory "High_order_TS_with_scaffold" in order to include the
##different libraries. It is computationally expensive
cd ../High_order_TS_with_scaffold/
codepath="simplicial_multivariate.py"
filename="./../Kaneko_CLM/trial_N50_T240_r175_eps012_008_03_0068_005.txt_kaneko"
javaplexpath="javaplex/"
python ${codepath} ${filename} -t 0 5 -p 5 -j ${javaplexpath} scaffold_ 
mv scaffold_gen/ ../Sample_results/
cd ../Example/


### Alternatively, you can run directly this code in the correponding folder:
## python simplicial_multivariate.py ../Kaneko_CLM/trial_N50_T240_r175_eps012_008_03_0068_005.txt_kaneko -t 0 5 -p 5 -j javaplex trial_

