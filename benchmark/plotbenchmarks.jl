using PkgBenchmark: BenchmarkTools
using Gadfly; import Cairo, Fontconfig
using DataFrames
using PkgBenchmark
using BenchmarkTools

versions = [
    "0.2.0",
    "0.2.1",
    "0.2.2",
    "0.2.3",
    "0.2.4",
    "0.2.5",
    "0.2.6",
    "0.2.7",
    "0.2.8",
    "0.2.9",
    "0.2.10",
    "master"
]

flatres(br::BenchmarkResults) = flatres(br.benchmarkgroup)
flatres(bg::BenchmarkGroup) = vcat([[(key*"⋅"*k, f) for (k,f) in flatres(subgroup)] for (key,subgroup) in bg]...)
flatres(t::BenchmarkTools.Trial) = [("",t)]

df = DataFrame([Dict(:version=>v,:benchmark=>b,:result=>r) for v in versions for (b,r) in flatres(readresults("bench_v$(v)"))])
df[:,:time] = (x->median(x).time).(df[:,:result])
df[:,:memory] = (x->median(x).memory).(df[:,:result])
df[:,:allocs] = (x->median(x).allocs).(df[:,:result])

df_stacked = stack(df, [:time, :memory, :allocs],)

Gadfly.with_theme(style(highlight_width=0.1mm)) do
    p = Gadfly.plot(
        df,
        x=:version, y=:time,
        ygroup=:benchmark,
        Geom.subplot_grid(Geom.line,free_y_axis=true,free_x_axis=true),
        #Scale.y_log10(),
        Scale.x_discrete(levels=unique(versions)),
        Guide.xlabel("version"),
        Guide.ylabel("time"),
    )
    img = PNG("all_bench.png", 20cm, 5cm*length(unique(df[:,:benchmark])), dpi=200)
    draw(img, p)
end