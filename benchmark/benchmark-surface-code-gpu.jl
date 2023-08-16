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

function get_stats()
    # Powers of 2 to benchmark
    powers = reverse(3:19) # start from the hardest one first! (fail fast)
    n_values = 2 .^ powers

    cpu_row_times = Float64[]
    gpu_row_times = Float64[]
    cpu_column_times = Float64[]
    gpu_column_times = Float64[]

    ccircuit = if eltype(circ) <: QuantumClifford.CompactifiedGate
        circ
    else
        compactify_circuit(circ)
    end

    CUDA.allowscalar(false)

    for n in ProgressBar(n_values)

        cpu_frames() = to_cpu(QuantumClifford._create_pauliframe(ccircuit; trajectories=n))
        gpu_frames() = to_gpu(QuantumClifford._create_pauliframe(ccircuit; trajectories=n))
        
        cpu_row() = pftrajectories(fastrow(cpu_frames()), ccircuit)
        cpu_column() = pftrajectories(fastcolumn(cpu_frames()), ccircuit)
        gpu_row() = pftrajectories(fastrow(gpu_frames()), ccircuit)
        gpu_column() = pftrajectories(fastcolumn(gpu_frames()), ccircuit)

        push!(cpu_row_times, @belapsed $cpu_row())
        push!(cpu_column_times, @belapsed $cpu_column())
        push!(gpu_row_times, @belapsed $gpu_row())
        push!(gpu_column_times, @belapsed $gpu_column())

        # note. this is also measuring the time it takes to move data to gpu
    end
    return DataFrame(
        n_values=n_values,
        cpu_row_times=cpu_row_times,
        cpu_column_times=cpu_column_times,
        gpu_row_times=gpu_row_times,
        gpu_column_times=gpu_column_times
    )
end

function load_stats()
    if isfile("gpu_cpu_comparison.csv")
        DataFrame(CSV.File("gpu_cpu_comparison.csv"))
    else
        get_stats()
    end
end

function save_stats(df)
    CSV.write("gpu_cpu_comparison.csv", df)
end

function plot_result(df)
    my_y_axis(a, b) = begin
        arr = [a..., b...]
        logarr = log.(arr)
        exp.(range(min(logarr...), stop=max(logarr...), length=8))
    end

    cpu_row_col = plot(df.n_values, [df.cpu_row_times, df.cpu_column_times], marker=:circle, label=["row" "column"], xlabel="n", ylabel="Execution Time (s)", title="cpu: row vs column", xticks=df.n_values[1:2:end], xscale=:log2, yticks=my_y_axis(df.cpu_row_times, df.cpu_column_times), yscale=:log10, legend=:topleft)
    gpu_row_col = plot(df.n_values, [df.gpu_row_times, df.gpu_column_times], marker=:circle, label=["row" "column"], xlabel="n", ylabel="Execution Time (s)", title="gpu: row vs column", xticks=df.n_values[1:2:end], xscale=:log2, yticks=my_y_axis(df.gpu_row_times, df.gpu_column_times), yscale=:log10, legend=:topleft)
    gpu_vs_cpu_row = plot(df.n_values, [df.gpu_row_times, df.cpu_row_times], marker=:circle, label=["gpu row" "cpu row"], xlabel="n", ylabel="Execution Time (s)", title="gpu vs cpu (row)", xticks=df.n_values[1:2:end], xscale=:log2, yticks=my_y_axis(df.gpu_row_times, df.cpu_row_times), yscale=:log10, legend=:topleft)
    gpu_vs_cpu_column = plot(df.n_values, [df.gpu_column_times, df.cpu_column_times], marker=:circle, label=["gpu column" "cpu column"], xlabel="n", ylabel="Execution Time (s)", title="gpu vs cpu (column)", xticks=df.n_values[1:2:end], xscale=:log2, yticks=my_y_axis(df.gpu_column_times, df.cpu_column_times), yscale=:log10, legend=:topleft)

    plot(cpu_row_col, gpu_row_col, gpu_vs_cpu_row, gpu_vs_cpu_column, layout=(2, 2), size=(1200, 800))
end

function main()
    df = load_stats()
    save_stats(df)
    plot_result(df)
    savefig("benchmark_surface_code_d5.png")
end
