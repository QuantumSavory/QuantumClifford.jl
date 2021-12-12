using PkgBenchmark: BenchmarkTools
using Gadfly; import Cairo, Fontconfig
using DataFrames
using PkgBenchmark
using BenchmarkTools

versions = [ # The files containing the benchmark results to be plotted.
]

flatres(br::BenchmarkResults) = flatres(br.benchmarkgroup)
flatres(bg::BenchmarkGroup) = vcat([[(key*"⋅"*k, f) for (k,f) in flatres(subgroup)] for (key,subgroup) in bg]...)
flatres(t::BenchmarkTools.Trial) = [("",t)]

df = DataFrame([Dict(:version=>v,:benchmark=>b,:result=>r) for v in versions for (b,r) in flatres(readresults(v))])
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
        Scale.y_continuous(minvalue=0),
        Scale.x_discrete(levels=unique(versions)),
        Guide.xlabel("version"),
        Guide.ylabel("time"),
    )
    img = PNG("all_bench.png", 20cm, 5cm*length(unique(df[:,:benchmark])), dpi=200)
    draw(img, p)
end

# Geometric mean
version_time_nomiss = stack(dropmissing(unstack(df, :benchmark, :version, :time)), Not(:benchmark))
df_geomean = combine(groupby(version_time_nomiss, :variable), :value => (x->exp(sum(log.(x))/length(x))))
p = Gadfly.plot(df_geomean, x=:variable, y=:value_function,
    Geom.line(), Guide.xlabel("version"), Guide.ylabel("time"), Scale.y_continuous(minvalue=0),
)
img = PNG("geommean_bench.png", dpi=200)
draw(img, p)