module SimplicialTS

using MAT
using Ripserer
using LinearAlgebra
using Combinatorics
using LoopVectorization
using PersistenceDiagrams
using Statistics
using Base.Threads
using HDF5
using DelimitedFiles
using Distances


export simplicial_complex_mvts,z_score, find_maximum_two_vecs!,create_data_structure,coherence_function,correction_for_coherence,create_simplicial_complex, fix_violations_and_compute_complexity, compute_complexity, load_data,load_data_mat,load_data_synthetic_kaneko, load_normaltxt, find_max_weight, compute_scaffold_fast, sliced_wasserstein


struct simplicial_complex_mvts
    raw_data::Matrix{Float64}
    num_ROI::Int64
    T::Int64
    ets_indexes::Dict{Int64,Vector{Int64}}
    ets_zscore::Matrix{Float64}
    ets_max::Vector{Float64}
    triplets_indexes::Dict{Int64,Vector{Int64}}
    triplets_indexes_reverse::Dict{Tuple{Int64,Int64,Int64},Int64}
    triplets_zscore::Matrix{Float64}
    triplets_max::Vector{Float64}
end

# function cityblock(x,y)
    
#     return sum(abs.(x-y))
# end

function sliced_wasserstein(PD1,PD2, M::Int=50)
    """
    Implementation of Sliced Wasserstein distance as described in 
    Sliced Wasserstein Kernel for Persistence Diagrams by Mathieu Carriere, Marco Cuturi, Steve Oudot (https://arxiv.org/abs/1706.03358)
    It is a Julia version of the python code https://github.com/scikit-tda/persim/blob/master/persim/sliced_wasserstein.py

    Parameters
    -----------
    PD1: Vector{Tuple{Real,Real}}
        Persistence diagram
    PD2: Vector{Tuple{Real,Real}}
        Persistence diagram
    M: Int, default is 50
        Iterations to run approximation.

    Returns
    --------
    sw::Float64
        Sliced Wasserstein distance between PD1 and PD2
    """
    
    diag_theta = Float32.([cos(0.25 * π), sin(0.25 * π)])

    l_theta1 = [dot(diag_theta, x) for x in PD1]
    l_theta2 = [dot(diag_theta, x) for x in PD2]

    PD_delta1 = [[sqrt(x^2 / 2.0)] .* [1, 1] for x in l_theta1]
    PD_delta2 = [[sqrt(x^2 / 2.0)] .* [1, 1] for x in l_theta2]

    # Compute the Sliced Wasserstein distance
    sw = 0.0
    theta = 0.5
    step = 1.0 / M
    for i in 1:M
        l_theta =Float32.([cos(theta * π), sin(theta * π)])

        V1 = [[dot(l_theta, x) for x in PD1] ; [dot(l_theta, x) for x in PD_delta2]]
        V2 = [[dot(l_theta, x) for x in PD2] ; [dot(l_theta, x) for x in PD_delta1]]

        sw += step * (cityblock(sort(V1), sort(V2)))
        theta += step
    end

    return sw
end


function z_score(data::Matrix{Float64}, N::Int64, T::Int64)::Matrix{Float64}
    zscore_raw=zeros((N,T))
    l=1
    @inbounds for row in eachrow(data)
        m=mean(row)
        s=std(row,mean=m,corrected=false)
        zscore_raw[l,:]=(row.-m)./s
        l+=1
    end
    return(zscore_raw)

end

function find_maximum_two_vecs!(vec1::Vector{Float64},vec2::Vector{Float64})
    @inbounds for i in eachindex(vec1)
        vec1[i]=max(abs(vec1[i]),abs(vec2[i]))
        end
end


function create_data_structure(data::Matrix{Float64})::simplicial_complex_mvts
    N,T=size(data)
    data=z_score(data,N,T)
    
    ##Compute the edges
    N_edges = div(N*(N-1),2)
    ets_zscore = similar(data,(N_edges, 2))
    ets_max=similar(data,T)
    ets_indexes=Dict{Int64, Vector{Int64}}()
    current_val=similar(data, T)
    l=1
    @inbounds for i in 1:N
         for j in  i+1:N
                @views current_val = data[i,:] .* data[j,:]
                m = mean(current_val)
                s = std(current_val,mean=m,corrected=false,)
                ets_zscore[l,1] = m
                ets_zscore[l,2] = s
                find_maximum_two_vecs!(ets_max,(current_val .- m) ./ s)
                ets_indexes[l] = [i,j]
                l+=1
        end
    end
        
    
    ##Compute the triplets
    N_triplets =binomial(N,3)
    triplets_zscore=similar(data,(N_triplets,2))
    triplets_max=similar(data,T)
    triplets_indexes=Dict{Int64, Vector{Int64}}()
    triplets_indexes_reverse=Dict{Tuple{Int64,Int64,Int64},Int64}()
    l=1
    @inbounds for i in 1:N
         for j in i+1:N
            for k in j+1:N
            @views current_val = data[i,:] .* data[j,:] .* data[k,:]
            m = mean(current_val)
            s = std(current_val,mean=m,corrected=false)
            triplets_zscore[l,1] = m
            triplets_zscore[l,2] = s
            find_maximum_two_vecs!(triplets_max,(current_val .- m) ./ s)
            triplets_indexes[l]=[i,j,k]
            triplets_indexes_reverse[(i,j,k)]=l
            l+=1
            end
        end
    end
    data_simplex=simplicial_complex_mvts(data,N,T,ets_indexes,ets_zscore,ets_max,triplets_indexes,triplets_indexes_reverse,triplets_zscore,triplets_max)
    return(data_simplex)
end


function coherence_function(vector::Vector{Float64})::Int64
    n = length(vector)
    temp = 0
    for el in vector
        temp += sign(el)
    end
    exponent = sign(n - abs(temp))
    res = (-1)^exponent
    return res
end

function correction_for_coherence(current_list_sign::Vector{Float64}, current_weight::Float64)::Float64
    # If the original signals are fully coherent, then the corresponding weight becomes positive, otherwise negative
    coherence = coherence_function(@views current_list_sign)
    # If all the signs are concordant then set the weight sign as positive, otherwise negative
    if coherence == 1
        weight_corrected = abs(current_weight)
    else
        weight_corrected = -abs(current_weight)
    end
    return weight_corrected
end


function create_simplicial_complex(simplicial_TS::simplicial_complex_mvts,t_current::Int64)::Vector{Tuple{Vector{Int64}, Float64}}
    list_simplices = Vector{Tuple{Vector{Int64}, Float64}}(undef, simplicial_TS.num_ROI + length(simplicial_TS.ets_indexes)+length(simplicial_TS.triplets_indexes))
    # compute the extremal weight
    m_weight = max(ceil(simplicial_TS.triplets_max[t_current]), ceil(simplicial_TS.ets_max[t_current]))
    
    #add the nodes to the list of simplices
    @inbounds for i in 1:simplicial_TS.num_ROI
        list_simplices[i]= ([i], m_weight)
    end
    
    # Adding the edges to the list of simplices
    j=simplicial_TS.num_ROI+1
    @inbounds for (i, indexes_ij) in simplicial_TS.ets_indexes
        c_mean = simplicial_TS.ets_zscore[i,1]
        c_std  = simplicial_TS.ets_zscore[i,2]
        raw_data_i = simplicial_TS.raw_data[indexes_ij[1],t_current]
        raw_data_j = simplicial_TS.raw_data[indexes_ij[2],t_current]
        weight_current = (raw_data_i * raw_data_j - c_mean) / c_std
        list_of_signs_edges = [raw_data_i, raw_data_j]
        @views weight_current_corrected = correction_for_coherence(list_of_signs_edges, weight_current)
        list_simplices[j] = (indexes_ij, weight_current_corrected)
        j+=1
    end    

#     println("j -> ",j)
    # Adding the triplets
    # Here I modify the signs of the weights, if it is fully coherent I assign a positive sign, otherwise negative
    @inbounds for (i,indexes_ijk) in simplicial_TS.triplets_indexes
        c_mean = simplicial_TS.triplets_zscore[i,1]
        c_std  = simplicial_TS.triplets_zscore[i,2]
        raw_data_i = simplicial_TS.raw_data[indexes_ijk[1],t_current]
        raw_data_j = simplicial_TS.raw_data[indexes_ijk[2],t_current]
        raw_data_k = simplicial_TS.raw_data[indexes_ijk[3],t_current]
        weight_current = (raw_data_i * raw_data_j * raw_data_k - c_mean) / c_std
        @views weight_current_corrected = correction_for_coherence([raw_data_i, raw_data_j,raw_data_k], weight_current)
        list_simplices[j] =(indexes_ijk, weight_current_corrected)
        j+=1
    end
    return(list_simplices)
    
end

# Function that remove all the violating triangles to create a proper filtration
function fix_violations_and_compute_complexity(sorted_simplices::Vector{Tuple{Vector{Int64}, Float64}},
 t_current::Int64, simplicial_TS::simplicial_complex_mvts,
 fid::HDF5.File,flag_scaffold::Int64, fid_triangles::HDF5.File, flag_triangles::Int64, flag_sliced_wasserstein::Int64)#::Tuple{Vector{Pair{Tuple{Int64, Vararg{Int64}}, Float64}},Vector{Tuple{Vector{Int64}, Float64, Float64}},Float64,Vector{Pair{Tuple{Int64, Vararg{Int64}}, Float64}}}
    # Sorting the simplices in a descending order according to the weights
    sort!(sorted_simplices, rev=true, by=x->x[2])
    
    # Remove the violations
    list_violating_triangles = Vector{Tuple{Vector{Int64}, Float64, Float64}}()
    list_simplices_for_filtration = Vector{Pair{Tuple{Int64, Vararg{Int64}}, Float64}}()
    set_simplices = Set()
    counter = 0
    triangles_count = 0
    violation_triangles = 0
    violation_triangles_negativeterms = 0
    CC_triangles_positive = 0
    total_CC_triangles = 0
    CC_triangles_negative = 0
    list_simplices_all = Vector{Pair{Tuple{Int64, Vararg{Int64}}, Float64}}()
    counter_simplices_all = 0

     # Loop over the sorted simplices, and flipping the sign of all the weights (so that the points in the persistence diagram are above the diagonal)
    for (index, i) in enumerate(sorted_simplices)
        simplices, weight = i

        # If the current simplex is an edge or a node, then I will immediately include it
        if length(simplices) <= 2
#             if length(simplices) == 1
            push!(list_simplices_all, (tuple(simplices...) => -weight))
#             else
#                 if weight != sorted_simplices[index - 1][2]
#                     counter_simplices_all += 1
#                     push!(list_simplices_all, (tuple(simplices...) => -weight))
#                 else
#                 push!(list_simplices_all, (tuple(simplices...) => -weight))
#                 end
#             end

            push!(list_simplices_for_filtration, (tuple(simplices...) => -weight))
            push!(set_simplices, simplices)
            counter += 1
        else
            # If the current simplex is a triplet, I check whether all the sub-simplices have been included.
            flag = 0
            for t in combinations(simplices, 2)
                if t in set_simplices
                    flag += 1
                end
            end
            

            # If all the sub-simplices already belong to the set, then I add it in the filtration
            if flag == 3
                push!(set_simplices, simplices)
                push!(list_simplices_for_filtration, (tuple(simplices...) => -weight))
                counter += 1
                if weight != sorted_simplices[index - 1][2]
                    counter_simplices_all += 1
                    # list_simplices_all[string(list(simplices))] = [string(counter_simplices_all), string(-weight)]
                    # list_simplices_all[Tuple(simplices)] = -weight
                    push!(list_simplices_all, (tuple(simplices...) => -weight))
                else
                    # list_simplices_all[string(list(simplices))] = [string(counter_simplices_all), string(-weight)]
                    # list_simplices_all[Tuple(simplices)] = -weight
                    push!(list_simplices_all, (tuple(simplices...) => -weight))
                end

                # Count the number of positive triangles that are in the filtration
                if weight >= 0
                    triangles_count += 1
                end
            else
                # Count the violations only for fully coherent state (--- or +++)
                if weight >= 0
                    violation_triangles += 1
                    push!(list_violating_triangles, (simplices, abs(weight), 3 - flag))
                else
                    violation_triangles_negativeterms += 1
                end
            end
        end
    end
     # Fraction of positive triangle discarderd (a.k.a. the hyper coherence)
    hyper_coherence = (1.0 * violation_triangles) / (triangles_count + violation_triangles)
#     println("Hyper Coherence: ",hyper_coherence)
#     println("Length:", length(list_simplices_all))
    # return(list_simplices_for_filtration, list_violating_triangles, hyper_coherence, list_simplices_all)

    compute_complexity(simplicial_TS,list_simplices_for_filtration, list_violating_triangles,
     hyper_coherence, list_simplices_all, t_current, flag_scaffold, fid, fid_triangles, flag_triangles, flag_sliced_wasserstein)
end

function compute_complexity(simplicial_TS::simplicial_complex_mvts,  list_simplices_positive::Vector{Pair{Tuple{Int64, Vararg{Int64}}, Float64}},
 list_violation_fully_coherence::Vector{Tuple{Vector{Int64}, Float64, Float64}},
  hyper_coherence::Float64,  list_filtration_scaffold::Vector{Pair{Tuple{Int64, Vararg{Int64}}, Float64}},
   t_current::Int64, flag_scaffold::Int64, fid::HDF5.File, fid_triangles::HDF5.File, flag_triangles::Int64, flag_sliced_wasserstein::Int64)
	num_nodes=simplicial_TS.num_ROI
    ##From here it starts the computation of hypercomplexity and scaffold
    dgms1_clean_FC=Vector{Tuple{Real,Real}}()
    dgms1_clean_CT=Vector{Tuple{Real,Real}}()
    dgms1_clean_FD=Vector{Tuple{Real,Real}}()
    
    dgms1=ripserer(Custom(list_simplices_positive),reps=1, alg=:cohomology, verbose=false)
    max_filtration_weight = find_max_weight(simplicial_TS,t_current)
    dgms1_clean = [(death(i) == Inf ? (birth(i), max_filtration_weight) : (birth(i), death(i))) for i in dgms1[2]]
    for i in dgms1[2]
        (a,b)=(birth(i), death(i))
        if b==Inf
            b=max_filtration_weight
        end
        # println(a, " ", b)
    end    
    for i in dgms1[2]
        (a,b)=(birth(i), death(i))
        if b==Inf
            b=max_filtration_weight
        end
        if a<0 
            if b<=0
                push!(dgms1_clean_FC,(a,b))
            else
                push!(dgms1_clean_CT,(a,b))
            end
        else
            push!(dgms1_clean_FD,(a,b))
        end
    end

    if flag_sliced_wasserstein ==1
        trivial_diagram=PersistenceDiagram([]; )
        
        complexity_FC = Wasserstein()(PersistenceDiagram(dgms1_clean_FC), trivial_diagram)
        complexity_CT = Wasserstein()(PersistenceDiagram(dgms1_clean_CT), trivial_diagram)
        complexity_FD = Wasserstein()(PersistenceDiagram(dgms1_clean_FD), trivial_diagram)
        #hyper_complexity = complexity_FC+ complexity_CT+complexity_FD
        hyper_complexity = Wasserstein()(PersistenceDiagram(dgms1_clean), trivial_diagram)
    else
        trivial_diagram = [(0.0,0.0),(1.0,1.0)]
        complexity_FC = sliced_wasserstein(dgms1_clean_FC, trivial_diagram)
        complexity_CT = sliced_wasserstein(dgms1_clean_CT, trivial_diagram)
        complexity_FD = sliced_wasserstein(dgms1_clean_FD, trivial_diagram)
        hyper_complexity = sliced_wasserstein(dgms1_clean, trivial_diagram)
    end
    
    if flag_scaffold !=0
        dgms1=ripserer(Custom(list_filtration_scaffold),reps=1, alg=:homology, verbose=false)
        f_tensor_scaffold_b,p_tensor_scaffold_b=compute_scaffold_fast(dgms1,num_nodes)
        if flag_scaffold == 1 # This i 
            fid[string(t_current)]=f_tensor_scaffold_b
        elseif flag_scaffold == 2
            fid[string(t_current)]=p_tensor_scaffold_b
        end
    end

    avg_edge_violation=0
    counter=0
    vec_triangles=zeros(floor(Int,num_nodes*(num_nodes-1)*(num_nodes-2)/6))

    for l in list_violation_fully_coherence
        avg_edge_violation+=l[3]
        counter+=1
        #println(l[1],l[2])
        current_index= simplicial_TS.triplets_indexes_reverse[tuple(l[1]...)]
        vec_triangles[current_index]=l[2]
    end
    #println(vec_triangles)
    avg_edge_violation/=counter

    ##Printing the triangles on a file
    if flag_triangles != 0
        fid_triangles[string(t_current)]=vec_triangles
    end


    println(t_current, " ", hyper_complexity, " ",complexity_FC, " ",complexity_CT, " ",complexity_FD, " ",hyper_coherence, " ",avg_edge_violation)
end




function load_data(path_single_file::String)::Matrix{Float64}
    extension_file = split(path_single_file, '.')[end]
    if extension_file == "mat"
        data = load_data_mat(path_single_file)
    elseif extension_file == "txt"
        data = load_data_synthetic_kaneko(path_single_file)
    elseif extension_file == "tx"
        data = load_normaltxt(path_single_file)
    end
    # print(size(data))
    return data
end

# Load brain data in .mat format (rows are ROI, columns are the time instants)
function load_data_mat(path_single_file::String)::Matrix{Float64}
    file_to_open = path_single_file
    data = matread(file_to_open)
    key_data = collect(keys(data))[end]
    data = data[key_data]
    return data
end



# Load the synthetic data generated from the Kaneko maps (coupled map lattices)
function load_data_synthetic_kaneko(path_single_file::String)::Matrix{Float64}
    file_to_open = path_single_file
    data = readdlm(file_to_open,Float64)
    data_cleaned = Any[]
    eps_list = [[0,data[1, end]]]
    
    # println(size(data[2:end,:]))
    # println(data[2:end])
    counter = 0
    for i in eachrow(data[2:end,:])
        # println(i)
        # If the data are in the shape of t 0 0 0... 0 eps then save the value of eps
        if i[1] == 0 && i[2] == 0
            N = counter
            counter = 0
            c = length(eps_list)
            push!(eps_list, [c * N, i[end]])
        else
            counter += 1
            push!(data_cleaned, i)
        end
    end
    M = size(data, 2) - 1
    data_cleaned=convert(Array{Float64,2},reduce(hcat,data_cleaned))

    return data_cleaned[2:end,:]
end

# Load synthetic data (format: columns represents independent time series)
function load_normaltxt(path_single_file::String)::Matrix{Float64}
    file_to_open = path_single_file
    data = readdlm(file_to_open)
    return permutedims(data)
end


function find_max_weight(list_simplices_positive, t)
    edges_abs_max = list_simplices_positive.ets_max[t]
    triplets_abs_max = list_simplices_positive.triplets_max[t]
    m = max(edges_abs_max, triplets_abs_max)
    return(m)
end

function compute_scaffold_fast(diagrams,num_regions)
    p_tensor_scaffold_b = zeros(num_regions,num_regions)
    f_tensor_scaffold_b = zeros(num_regions,num_regions)
    t = 0.0001
    eps = 0.01

    list_persistences_nb = Any[]
    list_persistences_b = Any[]
    list_persistences_c = Any[]
#     list_births_nb = Any[]
#     list_births_b = Any[]
#     list_births_c = Any[]

#     list_deaths_nb = Any[]
#     list_deaths_b = Any[]
#     list_deaths_c = Any[]

#     list_lengths_nb = Any[]
#     list_lengths_b = Any[]
#     list_lengths_c = Any[]

    count = 0
    c_fail = 0
    d = diagrams
    f = d[2].filtration
    persistences = []
    births = []
    deaths = []
    lengths = []
    for l in eachindex(d[2])
        g = d[2][l]
        p = persistence(g)
#         append!(persistences, p)
#         append!(births, birth(g))
#         append!(deaths, death(g))
#         append!(lengths,length(representative(g)))
        for rep in representative(g)
            count += 1
            i , j = vertices(rep)
            @inbounds p_tensor_scaffold_b[i,j] += p
            @inbounds f_tensor_scaffold_b[i,j] += 1
            @inbounds p_tensor_scaffold_b[j,i] += p
            @inbounds f_tensor_scaffold_b[j,i] += 1
        end
    end
#     push!(list_persistences_b, persistences)
#     push!(list_births_b, births)
#     push!(list_deaths_b, deaths)
#     push!(list_lengths_b, lengths)
    #return (p_tensor_scaffold_b,f_tensor_scaffold_b)
    return (f_tensor_scaffold_b,p_tensor_scaffold_b)
end
end