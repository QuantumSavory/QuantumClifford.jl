"""
Benchmark Clifford application
with respect to the choice data structure (CliffordOperator or CliffordColumnForm)
and array (column-major or row-major)
and bitpacking (in UInt8, or 16, or 32, or 64 bit types).
"""

using QuantumClifford
using StableRNGs
using BenchmarkTools
using DataFrames
using Gadfly; import Cairo, Fontconfig
using Serialization
using CSV
using Glob

## Pick a name for the logs

NAME = "ryzen"

## Run the benchmarks and save results

results = []
for n in [2,20], 
    int in [UInt8, UInt16, UInt32, UInt64], 
    s_order in [:column,:row], 
    c_order in [:column,:row], 
    clifftype in [CliffordOperator, CliffordColumnForm]

    @show n, int, s_order, c_order, clifftype
    N = 64*n-2
    rng = StableRNG(42)
    s64 = random_stabilizer(rng,N,N);
    phases = s64.phases;
    xzs64_colmajor = s64.xzs;
    xzs64_rowmajor = collect(s64.xzs')';

    xzs_rowmajor = collect(reinterpret(int, collect(xzs64_colmajor')))';
    xzs_colmajor = collect(xzs_rowmajor);
    s_col = Stabilizer(phases,N,xzs_colmajor);
    s_row = Stabilizer(phases,N,xzs_rowmajor);
    s = s_order == :column ? s_col : s_row

    c_raw = Stabilizer(zeros(UInt8, 2N), N, collect(reinterpret(int, collect(Destabilizer(random_stabilizer(rng,N,N)).tab.xzs'))'))
    if clifftype == CliffordOperator
        c = CliffordOperator(c_raw)
        if c_order == :row
            c = CliffordOperator(Stabilizer(c.tab.phases, N, collect(c.tab.xzs')'))
        end
    elseif clifftype == CliffordColumnForm
        c_medium = CliffordColumnForm(c_raw)
        c = CliffordColumnForm(c_medium.phases, N,
            collect(reinterpret(int, c_medium.xztox')'),
            collect(reinterpret(int, c_medium.xztoz')'))
        if c_order == :row
            c = CliffordColumnForm(c.phases, N, collect(c.xztox')', collect(c.xztoz')')
        end
    end
    bench_apply = @benchmark apply!(_s,_c; phases=false) setup=(_s=deepcopy($s);_c=deepcopy($c))
    push!(results, (;N, int, s_order, c_order, clifftype, bench=:apply_cliff, trial=bench_apply))
end

serialize("intsize_clifford_benchmarks_$(NAME).serialize", results)

## Prepare the results in a data frame and store to disk

DF = DataFrame(results)
for f in [:mean, :maximum, :minimum, :median]
    DF[!,f] .= (x->eval(f)(x).time).(DF[!,:trial])
end
DF[!,:bits] .= sizeof.(DF[!,:int])*8

CSV.write("intsize_clifford_benchmarks_$(NAME).csv", DF)

## Collate all benchmarks

dfs = []
for f in glob("intsize_clifford_benchmarks_*csv")
    df = CSV.read(f, DataFrame)
    df[!,:system] .= split(split(f, "_")[end],".")[1]
    push!(dfs,df)
end
DF = vcat(dfs...)
DF = DF[DF[!,:N].==1278,:]
Gadfly.with_theme(style(highlight_width=0.1mm)) do
    p = Gadfly.plot(
        DF,
        color=:clifftype,
        x=:bits, y=:median,
        xgroup=:s_order,
        ygroup=:c_order,
        Geom.subplot_grid(Geom.point,free_y_axis=false),
        Scale.y_log10(),
        Scale.x_discrete(levels=unique(DF[:,:bits])),
        Guide.xlabel("Int type\n(grouped by Stabilizer array order)"),
        Guide.ylabel("Time (ns)\n(grouped by Clifford array order)"),
        Guide.colorkey(title="Clifford\ndatastructure"),
    )
    img = PNG("bench.png", 20cm, 15cm, dpi=200)
    draw(img, p)
end