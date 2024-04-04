#!/bin/sh

path_simplicialTS="./SimplicialTS/"
julia_path=`which julia`


###Install the Package Compiler for Julia if it is not present
${julia_path} <<EOF
import Pkg; 
try
    @eval import PackageCompiler
    println("PackageCompiler is installed")
catch
    println("PackageCompiler is not installed!\n Installing it now...")
   	Pkg.add("PackageCompiler") 
end
EOF


###Install the Package SimplicialTS
cd $path_simplicialTS
${julia_path} <<EOF 
using Pkg
Pkg.activate(".")
Pkg.update()
using PackageCompiler
println("Compiling the package SimplicialTS...");
create_sysimage(["SimplicialTS"]; sysimage_path="Simplicial.so",precompile_execution_file="./src/precompile_simplicial.jl");
EOF
