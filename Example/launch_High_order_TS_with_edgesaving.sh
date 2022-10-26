#!/bin/sh
### Launch the code for the times 0-5 and using 5 cores 
### and save the magnitude of the projected violating triangles \Delta_v at the level of edges
### on the file "edges_projection.hd5"
codepath="../High_order_TS/simplicial_multivariate.py"
filename="./../Kaneko_CLM/trial_N50_T240_r175_eps012_008_03_0068_005.txt_kaneko"
outputfile_edges="../Sample_results/edges_projection"
outputfile="../Sample_results/results_T0_5.txt"
python ${codepath} ${filename} -t 0 5 -p 5 -s ${outputfile_edges}


