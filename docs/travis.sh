#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia -e 'if !(v"1.3.0" > VERSION > v"1.1.0") exit() end; using Pkg; Pkg.add("Documenter"); Pkg.add("Plots"); using SimpleClifford; cd(joinpath(dirname(pathof(SimpleClifford)), "..", "docs")); include("make.jl")';
fi
