#!/bin/sh

path_simplicialTS="./SimplicialTS/"
julia_path="/home/andrea/julia-1.8.5/bin/julia"

cd $path_simplicialTS
${julia_path} <<EOF 
using Pkg
Pkg.activate(".")
using PackageCompiler
println("Compiling the package SimplicialTS...")
# path_simplicial="$path_simplicialTS"
create_sysimage(["SimplicialTS"]; sysimage_path="Simplicial.so",precompile_execution_file="./src/precompile_simplicial.jl")
EOF


