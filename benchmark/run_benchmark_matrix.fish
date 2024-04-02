#!/usr/bin/env fish
for jlversion in nightly 1.10.0 1.9.0 1.8.0 1.7.0 1.6.0
    for tag in (git for-each-ref --sort=creatordate --format '%(tag)' refs/tags | tail -n+15)
        for nthreads in 1 4
            for retune in true
                ./run_ttfx_benchmarks.sh -v$jlversion -t$tag -n$nthreads
                ./run_benchmarks.sh -v$jlversion -t$tag -n$nthreads -r$retune
            end
        end
    end
end

