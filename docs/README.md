To generate the documentation run from the project folder
```
julia --project=docs -tauto -e '
            using Pkg
            Pkg.develop(PackageSpec(path=pwd()))
            Pkg.instantiate()';
julia --project=docs -tauto -i docs/make.jl
```

You might want to modify these lines of `make.jl`
```
warnonly = [:missing_docs],
linkcheck = true,
```
in order to make sure the build succeeds locally, if you are debugging any issues.