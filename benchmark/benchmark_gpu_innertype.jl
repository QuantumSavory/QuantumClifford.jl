using Revise
using QuantumClifford
using QuantumClifford: to_cpu, to_gpu
using CUDA
using BenchmarkTools
using Plots
using ProgressBars
using DataFrames
using CSV

include("surface-code-circuit.jl")

ints = [UInt8, UInt16, UInt32, UInt64, UInt128]
names = [:UInt8, :UInt16, :UInt32, :UInt64, :UInt128]

function get_stats()
    println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

    # Powers of 2 to benchmark
    powers = reverse(3:18) # start from the hardest one first! (fail fast)
    n_values = 2 .^ powers

    times = Dict([
        sm=>[]
        for (sm, tp) in zip(names, ints)
    ])

    ccircuit = if eltype(circ) <: QuantumClifford.CompactifiedGate
        circ
    else
        compactify_circuit(circ)
    end

    CUDA.allowscalar(false)


    for n in ProgressBar(n_values)    
        for (sm, T) in zip(names, ints)
            gpu_frames() = to_gpu(QuantumClifford._create_pauliframe(ccircuit; trajectories=n), T)
            gpu_column() = pftrajectories(fastcolumn(gpu_frames()), ccircuit), synchronize()
            @show n, T
            @elapsed gpu_column() # to eliminate compile time
            push!(times[sm], @belapsed $gpu_column())
        end
        # note. this is also measuring the time it takes to move data to gpu
    end
    df = DataFrame(times)
    df.n_values=n_values
    return df
end

function load_stats()
    if isfile("gpu_innertype_comparison.csv")
        DataFrame(CSV.File("gpu_innertype_comparison.csv"))
    else
        get_stats()
    end
end

function save_stats(df)
    CSV.write("gpu_innertype_comparison.csv", df)
end

function plot_result(df)
    my_y_axis(arrs) = begin
        arr = []
        for a in arrs
            append!(arr, a)
        end
        logarr = log.(arr)
        exp.(range(min(logarr...), stop=max(logarr...), length=8))
    end

    my_plots = []
    for col in names
        push!(my_plots, plot(df.n_values, df[!, col]))
    end
    plot(df.n_values, [df[!, col] for col in names], marker=:circle, label=reshape(["times $col" for col in names], (1, length(names))), xlabel="n", ylabel="Execution Time (s)", title="gpu inner type comparison (fast column)", xticks=df.n_values[1:2:end], xscale=:log2, yticks=my_y_axis([df[!, col] for col in names]), yscale=:log10, legend=:topleft, size=(1200, 800))
end

function main()
    df = load_stats()
    save_stats(df)
    plot_result(df)
    savefig("benchmark_gpu_innertype.png")
end
