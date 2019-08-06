#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia -e 'println(VERSION); if VERSION.minor!=2 println("not the desired version"); exit() end; println("building docs"); using Pkg; Pkg.add("Documenter"); Pkg.add("Plots"); using SimpleClifford; cd(joinpath(dirname(pathof(SimpleClifford)), "..", "docs")); include("make.jl")';
fi
