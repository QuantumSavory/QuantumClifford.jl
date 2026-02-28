
#=============================================================================#
function run_benchmark(
    benchmark_step!, synchronize, AT, path, title, xlabel; kwargs...
    )::Nothing

    parameters = collect(Iterators.product(block_sizes, batch_sizes))
    host_time = zeros(Float64, length(sizes_MiB))
    extrapolate_host = [false]
    host_fit = Vector{Vector}(undef, 1)
    device_time = zeros(Float64, length(sizes_MiB), size(parameters)...)
    extrapolate_device = similar(extrapolate_host, size(parameters)...)
    fill!(extrapolate_device, false)
    device_fit = similar(host_fit, size(parameters)...)

    error_flag = false
    for index in Base.OneTo(length(sizes_MiB))
        # This is intentional for short-circuit evaluation.
        error_flag = error_flag || benchmark_step!(
            synchronize, AT,
            index, parameters,
            host_time, extrapolate_host, host_fit,
            device_time, extrapolate_device, device_fit;
            kwargs...
            )
    end

    if error_flag
        error_string =
        "Unable to conclude benchmark ($title) within the current thresholds. \
        Consider relaxing the constraints of the provided configuration."
        @error error_string
        return nothing
    end

    mkpath(path)

    device_cat = [device_time[:, i] for i in CartesianIndices(parameters)]
    marked_labels = [x ? "Device*" : "Device" for x in extrapolate_device]
    label = hcat(
        [
        "$(marked_labels[i]) - (block, batch) = $x"
            for (i, x) in enumerate(parameters)
            ]...,
        extrapolate_host[1] ? "Host*" : "Host"
        )

    plot(
        sizes_MiB, 10^3 .* hcat(device_cat..., host_time);
        plot_style..., yscale = :log10, title = title, label = label,
        xlabel = xlabel, ylabel = "Runtime (ms)"
        )
    savefig("$path/runtime.$file_format")

    plot(
        sizes_MiB, hcat(map(x -> host_time ./ x, device_cat)...);
        plot_style..., title = title, label = hcat(label[1 : (end - 1)]...),
        xlabel = xlabel, ylabel = "Ratio (host/device)"
        )
    savefig("$path/ratio.$file_format")

    return nothing

end
#=============================================================================#
