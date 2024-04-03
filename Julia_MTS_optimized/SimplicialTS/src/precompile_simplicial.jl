current_path_simplicial_TS=pwd()
push!(LOAD_PATH,"$current_path_simplicial_TS")
using SimplicialTS
using HDF5
flag_scaffold=1
flag_triangles=1
flag_sliced_wasserstein = 0
t=1
data=load_data(current_path_simplicial_TS*"/data/TS_Schaefer60S_gsr_bp_z.mat");
simplicial_TS=create_data_structure(data);
list_all_simplices=create_simplicial_complex(simplicial_TS,t);
FID= h5open("trial.hd5", "w")
FID_triangles=h5open("trial_triangles.hd5","w")
fix_violations_and_compute_complexity(list_all_simplices, t, simplicial_TS, FID, flag_scaffold, FID_triangles, flag_triangles, flag_sliced_wasserstein)
close(FID)
close(FID_triangles)
