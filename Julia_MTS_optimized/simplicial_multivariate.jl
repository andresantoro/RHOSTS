current_path_simplicial_TS = pwd()
push!(LOAD_PATH, "$current_path_simplicial_TS" * "/SimplicialTS/")
using SimplicialTS
using Base.Threads
using HDF5
using ArgParse

function parse_commandline(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "filename"
            help = "multivariate time series input file"
            required = true
        "-t", "--time"
            help = "time interval for the multivariate time series"
            nargs=2
            default = [1,1]
        "-s", "--scaffold_name"
            help = "scaffold output name "
            default = "None"
        "-f", "--flag_scaffold"
            help = "scaffold flag (1-> frequency, 2 -> persistence)"
            arg_type = Int
            default = 0
        "-o", "--violating_triangles_name"
            help = "violating triangles output name "
            default = "None"
        "-w", "--wasserstein"
            help = "Use the wasserstein distance (default Sliced Wasserstein)"
            action = :store_true
            #default =  0
        "-v", "--verbose"
            help = "Activate verbose"
            # default =  0
            action = :store_true
    end

    return parse_args(s)
end




#### MAIN ### 

input_data=parse_commandline(ARGS)
#println(input_data)
filename=input_data["filename"]

if typeof(input_data["time"][1]) != Int64
    t0=parse(Int64,input_data["time"][1])
    t1=parse(Int64,input_data["time"][2])
else
    t0=input_data["time"][1]
    t1=input_data["time"][2]
end
scaffold_name=input_data["scaffold_name"]
flag_scaffold=input_data["flag_scaffold"]
if scaffold_name != "None" && flag_scaffold ==0
    flag_scaffold = 1 ## frequency scaffold as default 
end
triangles_name=input_data["violating_triangles_name"]
flag_triangles=0
if triangles_name!="None"
    flag_triangles=1
end
flag_sliced_wasserstein=1*input_data["wasserstein"]
verbose = input_data["verbose"]


if verbose == 1
    println(stderr, "Filename inserted: ", filename)
    if t0 != t1
        println(stderr, "Interval time points: ", t0, "-", t1)
    end
    println(stderr, "Scaffold Flag: ", flag_scaffold, ";   Scaffold name: ", scaffold_name)
    println(stderr, "Triangles Flag: ", flag_triangles, ";   Triangles name: ", triangles_name)
end


### Data loading ###
data = load_data(filename)
if t0 == t1
    t0 = t0
    t1 = size(data)[2]
end

simplicial_TS = create_data_structure(data)
if flag_scaffold != 0
    FID = h5open(scaffold_name * ".hd5", "w")
else
    FID= h5open("None.hd5","w")
end
if flag_triangles != 0
    FID_triangles = h5open(triangles_name * ".hd5", "w")
else
    FID_triangles= h5open("None1.hd5","w")
end

alldata= Array{Union{Nothing, String}}(nothing, t1-t0+1)
@threads for t in t0:t1
    list_all_simplices = create_simplicial_complex(simplicial_TS, t)
    local output=fix_violations_and_compute_complexity(list_all_simplices, t, simplicial_TS, FID, flag_scaffold, FID_triangles, flag_triangles, flag_sliced_wasserstein)
    alldata[t]=join(output, ' ')
end

# ##Printing the HOindicators on the output
for HO_values in alldata
    println(HO_values)
end
# Closing the files
close(FID)
close(FID_triangles)
