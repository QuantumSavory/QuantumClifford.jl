using QuantumClifford
using StableRNGs
using BenchmarkTools
using DataFrames
using Gadfly; import Cairo, Fontconfig
using Serialization
using CSV
using Glob

## Pick a name for the logs

NAME = "NAME"

## Run the benchmarks and save results

results = []
for n in 1:20, int in [UInt8, UInt16, UInt32, UInt64], order in [:column,:row]
    @show n, int, order
    N = 64*n-2
    rng = StableRNG(42)
    s64 = random_stabilizer(rng,N,N);
    phases = s64.phases;
    xzs64_colmajor = s64.xzs;
    xzs64_rowmajor = collect(s64.xzs')';
    p64 = random_pauli(rng,N);

    p = PauliOperator(p64.phase, N, collect(reinterpret(int,p64.xz)));
    xzs_rowmajor = collect(reinterpret(int, collect(xzs64_colmajor')))';
    xzs_colmajor = collect(xzs_rowmajor);
    s_col = Stabilizer(phases,N,xzs_colmajor);
    s_row = Stabilizer(phases,N,xzs_rowmajor);
    s = order == :column ? s_col : s_row
    bench_apply_pauli = @benchmark apply!(_s,_p) setup=(_s=deepcopy($s);_p=deepcopy($p))
    bench_canon = @benchmark canonicalize!(_s) setup=(_s=deepcopy($s))
    bench_canon_nophase = @benchmark canonicalize!(_s, phases=false) setup=(_s=deepcopy($s))

    push!(results, (N=N, int=int, order=order, bench=:apply_pauli, trial=bench_apply_pauli))
    push!(results, (N=N, int=int, order=order, bench=:canon, trial=bench_canon))
    push!(results, (N=N, int=int, order=order, bench=:canon_nophase, trial=bench_canon_nophase))
end

serialize("intsize_benchmarks_$(NAME).serialize", results)

## Prepare the results in a data frame and store to disk

DF = DataFrame(results)
for f in [:mean, :maximum, :minimum, :median]
    DF[!,f] .= (x->eval(f)(x).time).(DF[!,:trial])
end
DF[!,:bits] .= sizeof.(DF[!,:int])*8

CSV.write("intsize_benchmarks_$(NAME).csv", DF)

## Collate all benchmarks

dfs = []
for f in glob("intsize_benchmarks_*csv")
    df = CSV.read(f, DataFrame)
    df[!,:system] .= split(split(f, "_")[end],".")[1]
    push!(dfs,df)
end
DF = vcat(dfs...)
Gadfly.with_theme(style(highlight_width=0.1mm)) do
    p = Gadfly.plot(
        DF,
        color=:bits,
        x=:N, y=:median,
        shape=:system,
        xgroup=:order,
        ygroup=:bench,
        Geom.subplot_grid(Geom.point,free_y_axis=true),
        Scale.color_discrete_hue(),
        Scale.shape_discrete(),
        Scale.y_log2(),
        Scale.x_log2(),
        Guide.colorkey(title="Int\ntype"),
        Guide.shapekey(title="System"),
        Guide.xlabel("Number of Qubits\n(grouped by array order)"),
        Guide.ylabel("Time (ns)\n(grouped by benchmark)"),
    )
    img = PNG("bench.png", 20cm, 15cm, dpi=200)
    draw(img, p)
end
for df in dfs
    DF = vcat(dfs...)
    Gadfly.with_theme(style(highlight_width=0.1mm)) do
        p = Gadfly.plot(
            df,
            color=:bits,
            x=:N, y=:median,
            xgroup=:order,
            ygroup=:bench,
            Geom.subplot_grid(Geom.point,free_y_axis=true),
            Scale.color_discrete_hue(),
            Scale.y_log2(),
            Scale.x_log2(),
            Guide.colorkey(title="Int\ntype"),
            Guide.xlabel("Number of Qubits\n(grouped by array order)"),
            Guide.ylabel("Time (ns)\n(grouped by benchmark)"),
        )
        system = df[1,:system]
        img = PNG("bench_$(system).png", 20cm, 15cm, dpi=200)
        draw(img, p)
    end
end

## Run this with custom settings for a quick-and-dirty benchmark

N = 3*64-2   # Pick a number of qubits
int = UInt32 # Pick a data type

rng = StableRNG(42)
s64 = random_stabilizer(rng,N,N);
phases = s64.phases;
xzs64_colmajor = s64.xzs;
xzs64_rowmajor = collect(s64.xzs')';
p64 = random_pauli(rng,N);

p = PauliOperator(p64.phase, N, collect(reinterpret(int,p64.xz)));
xzs_rowmajor = collect(reinterpret(int, collect(xzs64_colmajor')))';
xzs_colmajor = collect(xzs_rowmajor);
s_col = Stabilizer(phases,N,xzs_colmajor);
s_row = Stabilizer(phases,N,xzs_rowmajor);

# Both column-major and row-major orders will be tested
@benchmark canonicalize!(_s, phases=false) setup=(_s=deepcopy($s_row))
@benchmark canonicalize!(_s, phases=false) setup=(_s=deepcopy($s_col))
@benchmark canonicalize!(_s) setup=(_s=deepcopy($s_row))
@benchmark canonicalize!(_s) setup=(_s=deepcopy($s_col))
@benchmark apply!(_s,_p) setup=(_s=deepcopy($s_row);_p=deepcopy($p))
@benchmark apply!(_s,_p) setup=(_s=deepcopy($s_col);_p=deepcopy($p))
