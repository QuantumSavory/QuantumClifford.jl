using CairoMakie
using AlgebraOfGraphics
using QuantumClifford
using PkgBenchmark
using DataFrames
using DataFramesMeta
using Statistics
using Glob

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
    threads = m.captures[2]
    julia = m.captures[3]
    pkgimport = minimum(parse.(Float64, readlines(f)))
    task = minimum(parse.(Float64, readlines("task_tag=$(tag)-nthreads=$(threads)-julia=$(julia).ttfxresults")))
    push!(records, (;tag, threads, julia, pkgimport, task))
end

##

df = DataFrame(records)
df[!,:total] .= df[:,:task] .+ df[:,:pkgimport]
df = stack(df, [:total, :task, :pkgimport], variable_name=:x, value_name=:time)

layers = visual(Stairs, step=:post)
axis = (;width = 500, height = 400, xticklabelrotation=-pi/3)
facet = (;linkyaxes=:minimal)

##

benchmark_data = data(df)
mappings = mapping(:tag, :time, color=:julia, row=:x, col=:threads)
f = draw(benchmark_data * layers * mappings; facet, axis)
save("../benchmarks_ttfx.png",f)
