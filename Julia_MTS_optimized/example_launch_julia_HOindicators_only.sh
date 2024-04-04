#!/bin/sh

juliapath=`which julia`
compiled_code="./SimplicialTS/Simplicial.so"
nthreads=6
input_file="./SimplicialTS/data/TS_Schaefer60S_gsr_bp_z.mat"



###Launching the code 
#time julia -J ${compiled_code} --threads ${nthreads} simplicial_multivariate.jl ${input_file}

##Launching the code only for time 1-20 using 6 threads
time julia -J ${compiled_code} --threads ${nthreads} simplicial_multivariate.jl ${input_file} -t 1 20

