"""
Benchmark quantum_mallows from randoms.jl.
"""

using QuantumClifford: quantum_mallows_bigint, quantum_mallows_int, quantum_mallows_float
using BenchmarkTools
using DataFrames
using Gadfly; import Cairo, Fontconfig
using ProgressLogging
using Random: GLOBAL_RNG

## Run the benchmarks and save results

results = []
@progress for n in [2:2:30; 31:20:200; 200:50:500; 500:100:800][end:-1:1]
    println(n)
    bench = @benchmark quantum_mallows_bigint($GLOBAL_RNG, $n)
    push!(results, (;type=:bigint, trial=bench, n=n))
    n>=500 && continue
    bench = @benchmark quantum_mallows_float($GLOBAL_RNG, $n)
    push!(results, (;type=:float, trial=bench, n=n))
    n>=30 && continue
    bench = @benchmark quantum_mallows_int($GLOBAL_RNG, $n)
    push!(results, (;type=:int, trial=bench, n=n))
end

## Prepare the results in a data frame

DF = DataFrame(results)
for f in [:mean, :maximum, :minimum, :median]
    DF[!,f] .= (x->eval(f)(x).time).(DF[!,:trial])
end

## Collate all benchmarks

Gadfly.with_theme(style(highlight_width=0.1mm)) do
    p = Gadfly.plot(
        DF,
        color=:type,
        x=:n, y=:median,
        Scale.x_log10(),
        Scale.y_log10(),
        Guide.xlabel("Bits"),
        Guide.ylabel("Time (ns)"),
        Guide.colorkey(title="Int Type"),
    )
    img = PNG("bench_mallows.png", 20cm, 15cm, dpi=200)
    draw(img, p)
end