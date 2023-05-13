To generate the documentation run from the project folder
```
julia --project=docs -tauto -e '
            using Pkg
            Pkg.develop(PackageSpec(path=pwd()))
            Pkg.instantiate()';
julia --project=docs -tauto -i docs/make.jl
```