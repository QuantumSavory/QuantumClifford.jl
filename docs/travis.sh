#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --project=docs/ -e 'println(VERSION); if VERSION.minor!=3 println("not the desired version"); exit() end; println("building docs"); using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd())); using SimpleClifford; cd(joinpath(dirname(pathof(SimpleClifford)), "..", "docs")); include("make.jl")';
fi
