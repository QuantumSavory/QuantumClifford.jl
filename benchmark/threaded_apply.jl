"""
Benchmark threaded Clifford application.
"""

using QuantumClifford
using Nemo
using StableRNGs
using BenchmarkTools
using DataFrames
using Gadfly; import Cairo, Fontconfig
using Serialization
using CSV
using Glob

## Pick a name for the logs

NAME = "$(length(ARGS)>0 ? ARGS[1] : "")_$(Threads.nthreads())"

## Run the benchmarks and save results

results = []
for _N in [2048:-128:1024+128; 1024:-32:256; 256-32:-16:16; 14:-2:2]
    global N = _N              # TODO understand why global is necessary for eval to work
    global rng = StableRNG(42)
    for s in [
        :(random_destabilizer(rng, N)),
        :(one(Stabilizer, N))
       ]
       s_e = eval(s)
       for  gate in [
                :(random_clifford(rng,N),),
                :(tensor_pow(CNOT,N÷2),),
                :(CNOT,[N÷2,N]),
                :(Hadamard,[N÷2]),
                :(sCNOT(N÷2,N),),
                :(sHadamard(N÷2),)
            ]
            gate_e = eval(gate)
            println((Threads.nthreads(), N, gate, s))
            bench_apply = @benchmark apply!(_s,_g...) setup=(_s=deepcopy($s_e);_g=deepcopy($gate_e))
            push!(results, (;N, gate=string(gate), s=string(s), threads=Threads.nthreads(), trial=bench_apply))
        end
    end
end

serialize("threaded_apply_benchmarks_$(NAME).serialize", results)

## Prepare the results in a data frame and store to disk

DF = DataFrame(results)
for f in [:mean, :maximum, :minimum, :median]
    DF[!,f] .= (x->eval(f)(x).time).(DF[!,:trial])
end

CSV.write("threaded_apply_benchmarks_$(NAME).csv", DF)

## Collate all benchmarks

dfs = []
for f in glob("threaded_apply_benchmarks_*csv")
    df = CSV.read(f, DataFrame)
    df[!,:system] .= split(split(f, "_")[end],".")[1]
    push!(dfs,df)
end
DF = vcat(dfs...)
DF[!, :threadtime] .= DF[:,:minimum] .* DF[:,:threads]
Gadfly.with_theme(style(highlight_width=0.1mm)) do
    p = Gadfly.plot(
        DF,
        color=:threads,
        x=:N, y=:minimum,
        xgroup=:s,
        ygroup=:gate,
        Geom.subplot_grid(Geom.line,free_y_axis=true),
        Scale.y_log10(),
        Scale.x_log10(),
        Scale.color_discrete(),
        Guide.xlabel("Number of qubits\n(grouped by whether the stabilizer is random or diagonal)"),
        Guide.ylabel("Time (ns)\n(grouped by type of Clifford operation)"),
        Guide.colorkey(title="Number of \nthreads"),
    )
    img = PNG("bench_threaded_apply.png", 18cm, 6*8cm, dpi=200)
    draw(img, p)
end
Gadfly.with_theme(style(highlight_width=0.1mm)) do
    p = Gadfly.plot(
        DF,
        color=:threads,
        x=:N, y=:threadtime,
        xgroup=:s,
        ygroup=:gate,
        Geom.subplot_grid(Geom.line,free_y_axis=true),
        Scale.y_log10(),
        Scale.x_log10(),
        Scale.color_discrete(),
        Guide.xlabel("Number of qubits\n(grouped by whether the stabilizer is random or diagonal)"),
        Guide.ylabel("Time × Threads (ns) \n(grouped by type of Clifford operation)"),
        Guide.colorkey(title="Number of \nthreads"),
    )
    img = PNG("bench_threaded_apply_threadtime.png", 18cm, 6*8cm, dpi=200)
    draw(img, p)
end
