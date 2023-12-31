#!/usr/bin/env bash
while getopts t:v:n: flag
do
    case "${flag}" in
        t) tag=${OPTARG};;
        v) jlversion=${OPTARG};;
        n) nthreads=${OPTARG};;
    esac
done
juliaup default $jlversion
mydir=$(mktemp -d)
git checkout $tag --force
resfile="tag=$tag-nthreads=$nthreads-julia=$jlversion.ttfxresults"
export JULIA_NUM_THREADS=$nthreads
cat << EOF | julia
using Pkg;
Pkg.activate("$mydir");
pkg"develop .."
Pkg.precompile()
EOF
for i in {1..4}
do
cat << EOF | julia -t$nthreads | tail -n1 >> using_$resfile
using Pkg;
Pkg.activate("$mydir");
println(@elapsed @eval using QuantumClifford)
EOF
cat << EOF | julia -t$nthreads | tail -n1 >> task_$resfile
using Pkg;
Pkg.activate("$mydir");
using QuantumClifford
println(@elapsed @eval canonicalize!(apply!(apply!(random_stabilizer(3),random_clifford(3)),sCNOT(3,1))))
EOF
done
git checkout master --force