#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia -e 'using Pkg; Pkg.add("Documenter")';
    julia -e 'using SimpleClifford; cd(joinpath(dirname(pathof(SimpleClifford)), "..", "docs")); include("make.jl")';
fi
