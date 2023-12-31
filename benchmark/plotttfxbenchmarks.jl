using CairoMakie
using AlgebraOfGraphics
using QuantumClifford
using PkgBenchmark
using DataFrames
using DataFramesMeta
using Statistics
using Glob

##
cd("QuantumCliffordBenchmarksLog/ttfxlogs")
##
files = glob("using_*.ttfxresults")
records = []
for f in files
    m = match(r"using_tag=(v[\d\.]*)-nthreads=(\d)-julia=([\d\.\w]*).ttfxresults", f)
    if isnothing(m)
        @warn "problem with $(f)"
        continue
    end
    tag = m.captures[1]
    version = VersionNumber(tag[2:end])
    threads = m.captures[2]
    juliastr = m.captures[3]
    julia = juliastr=="nightly" ? v"1.999" : VersionNumber(m.captures[3])
    pkgimportsamples = parse.(Float64, readlines(f))
    pkgimport = isempty(pkgimportsamples) ? missing : minimum(parse.(Float64, readlines(f)))
    tasksamples = parse.(Float64, readlines("task_tag=$(tag)-nthreads=$(threads)-julia=$(juliastr).ttfxresults"))
    task = isempty(tasksamples) ? missing : minimum(tasksamples)
    push!(records, (;version, threads, julia, pkgimport, task))
end

##

df = DataFrame(records)
maxjversion = maximum(df[df.julia .!= v"1.999",:julia])
df[df.julia .== v"1.999", :julia] .= VersionNumber(maxjversion.major, maxjversion.minor+1, 0)
df[!,:total] .= df[:,:task] .+ df[:,:pkgimport]
df = DataFrames.stack(df, [:total, :task, :pkgimport], variable_name=:x, value_name=:time)
sort!(df) # by first column (version)

df = df[df.time .!== missing,:] # otherwise the ticks are messed up, probably a bug

layers = visual(Stairs, step=:center)
maxy = ceil(maximum(skipmissing(df.time)))
axis = (;width = 800, height = 400, xticklabelrotation=-pi/2, yticks=(0:0.5:maxy))
facet = (;
    linkyaxes=:minimal,
    )

##

benchmark_data = data(df)
mappings = mapping(:version, :time, color=:julia, row=:x, col=:threads)
f = draw(benchmark_data * layers * mappings; facet, axis)
save("../benchmarks_ttfx.png",f)
