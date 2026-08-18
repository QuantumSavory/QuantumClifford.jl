using Printf
using Statistics

const RAW_COLUMNS = [
    "comparison", "label", "checkout", "commit", "build", "sample", "scenario",
    "build_seconds", "cache_bytes", "import_seconds", "first_seconds",
    "first_compile_seconds", "first_recompile_seconds", "total_seconds",
    "warm_seconds", "warm_compile_seconds", "warm_recompile_seconds",
]
const TIMING_COLUMNS = [
    "build_seconds", "import_seconds", "first_seconds", "first_compile_seconds",
    "first_recompile_seconds", "total_seconds", "warm_seconds",
    "warm_compile_seconds", "warm_recompile_seconds",
]
const TOKEN_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.-]*$"

length(ARGS) == 4 || error("usage: summarize.jl RAW_TSV SUMMARY_TSV BUILD_SUMMARY_TSV SUMMARY_MD")
raw_path, summary_path, build_summary_path, markdown_path = ARGS

lines = readlines(raw_path)
isempty(lines) && error("raw result file is empty: $raw_path")
header = String.(split(first(lines), '\t'; keepempty=true))
header == RAW_COLUMNS || error("raw result header does not match the expected schema")
rows = map(Iterators.drop(lines, 1)) do line
    fields = split(line, '\t'; keepempty=true)
    length(fields) == length(header) || error("raw result row has $(length(fields)) fields; expected $(length(header))")
    Dict(zip(header, String.(fields)))
end
isempty(rows) && error("raw result file contains no measurements: $raw_path")

number(row, name) = parse(Float64, row[name])
median_iqr(values) = (median(values), quantile(values, 0.75) - quantile(values, 0.25))

function positive_integer(value, description)
    parsed = tryparse(Int, value)
    !isnothing(parsed) && parsed > 0 || error("$description must be a positive integer")
    return parsed
end

comma_list(value) = isempty(value) ? String[] : String.(split(value, ','))

function mapping_list(value, description)
    mappings = Pair{String,String}[]
    for specification in comma_list(value)
        fields = split(specification, '='; limit=2)
        length(fields) == 2 && all(!isempty, fields) || error("invalid $description entry: $specification")
        push!(mappings, String(first(fields)) => String(last(fields)))
    end
    return mappings
end

function read_design(raw_path)
    metadata_path = joinpath(dirname(raw_path), "metadata.txt")
    isfile(metadata_path) || error("benchmark metadata is missing: $metadata_path")
    metadata = Dict{String,String}()
    variants = String[]
    for (line_number, line) in enumerate(eachline(metadata_path))
        fields = split(line, '='; limit=2)
        length(fields) == 2 || error("malformed metadata line $line_number")
        key, value = String.(fields)
        haskey(metadata, key) && error("duplicate metadata key: $key")
        metadata[key] = value
        variant = match(r"^variant\.([A-Za-z0-9][A-Za-z0-9_.-]*)\.commit$", key)
        isnothing(variant) || push!(variants, only(variant.captures))
    end

    builds = positive_integer(get(metadata, "builds", ""), "metadata builds")
    samples = positive_integer(
        get(metadata, "recorded_samples_per_build", ""),
        "metadata recorded_samples_per_build",
    )
    scenarios = comma_list(get(metadata, "scenarios", ""))
    isempty(scenarios) && error("metadata scenarios must not be empty")
    all(scenario -> occursin(TOKEN_PATTERN, scenario), scenarios) || error("metadata contains an invalid scenario token")
    allunique(scenarios) || error("metadata contains duplicate scenarios")
    length(variants) >= 2 || error("metadata must describe at least two variants")
    allunique(variants) || error("metadata contains duplicate variants")

    checkouts = Dict{String,String}()
    commits = Dict{String,String}()
    for label in variants
        checkouts[label] = get(metadata, "variant.$label.checkout", "")
        commits[label] = get(metadata, "variant.$label.commit", "")
        isempty(checkouts[label]) && error("metadata has no checkout for variant $label")
        isempty(commits[label]) && error("metadata has no commit for variant $label")
    end

    baselines = Dict{String,String}()
    for (candidate, baseline) in mapping_list(get(metadata, "candidate_baselines", ""), "candidate baseline")
        candidate in variants[2:end] || error("baseline map has unknown candidate: $candidate")
        baseline in variants || error("baseline map has unknown baseline: $baseline")
        candidate != baseline || error("candidate cannot be its own baseline: $candidate")
        haskey(baselines, candidate) && error("duplicate baseline map for candidate: $candidate")
        baselines[candidate] = baseline
    end

    return (; builds, samples, scenarios, variants, checkouts, commits, baselines)
end

function validate_rows(rows, design)
    expected_keys = Set{NTuple{5,String}}()
    default_baseline = first(design.variants)
    for comparison in design.variants[2:end]
        baseline = get(design.baselines, comparison, default_baseline)
        for label in (baseline, comparison), build in 1:design.builds,
            sample in 1:design.samples, scenario in design.scenarios
            push!(expected_keys, (comparison, label, string(build), string(sample), scenario))
        end
    end

    actual_keys = Set{NTuple{5,String}}()
    build_values = Dict{NTuple{3,String},Tuple{String,String}}()
    for (row_number, row) in enumerate(rows)
        key = (
            row["comparison"], row["label"], row["build"], row["sample"], row["scenario"],
        )
        key in actual_keys && error("duplicate raw sample key at data row $row_number: $key")
        push!(actual_keys, key)
        positive_integer(row["build"], "raw build at data row $row_number")
        positive_integer(row["sample"], "raw sample at data row $row_number")
        positive_integer(row["cache_bytes"], "raw cache_bytes at data row $row_number")

        label = row["label"]
        haskey(design.checkouts, label) || error("unknown raw variant label at data row $row_number: $label")
        row["checkout"] == design.checkouts[label] || error("checkout mismatch at data row $row_number")
        row["commit"] == design.commits[label] || error("commit mismatch at data row $row_number")
        for field in TIMING_COLUMNS
            value = number(row, field)
            isfinite(value) || error("non-finite $field at data row $row_number")
            value >= 0 || error("negative $field at data row $row_number")
        end

        import_time = number(row, "import_seconds")
        first_time = number(row, "first_seconds")
        total_time = number(row, "total_seconds")
        minimum_total = import_time + first_time
        total_tolerance = 1e-9 * max(1.0, total_time, minimum_total)
        total_time + total_tolerance >= minimum_total ||
            error("total_seconds is shorter than import_seconds plus first_seconds at data row $row_number")
        tolerance = 1e-6
        number(row, "first_compile_seconds") <= first_time + tolerance ||
            error("first compile time exceeds first task time at data row $row_number")
        number(row, "first_recompile_seconds") <= number(row, "first_compile_seconds") + tolerance ||
            error("first recompile time exceeds first compile time at data row $row_number")
        warm_time = number(row, "warm_seconds")
        number(row, "warm_compile_seconds") <= warm_time + tolerance ||
            error("warm compile time exceeds warm task time at data row $row_number")
        number(row, "warm_recompile_seconds") <= number(row, "warm_compile_seconds") + tolerance ||
            error("warm recompile time exceeds warm compile time at data row $row_number")

        build_key = (row["comparison"], label, row["build"])
        value = (row["build_seconds"], row["cache_bytes"])
        if haskey(build_values, build_key) && build_values[build_key] != value
            error("cache build values vary within $build_key")
        end
        build_values[build_key] = value
    end

    missing = setdiff(expected_keys, actual_keys)
    unexpected = setdiff(actual_keys, expected_keys)
    isempty(missing) && isempty(unexpected) || error(
        "raw result design does not match metadata: $(length(missing)) missing and $(length(unexpected)) unexpected sample keys"
    )
end

design = read_design(raw_path)
validate_rows(rows, design)
comparisons = unique(row["comparison"] for row in rows)

function comparison_labels(comparison)
    labels = unique(row["label"] for row in rows if row["comparison"] == comparison)
    length(labels) == 2 || error("comparison $comparison must have one baseline and one candidate")
    comparison in labels || error("comparison $comparison has no matching candidate label")
    baseline_labels = filter(!=(comparison), labels)
    length(baseline_labels) == 1 || error("comparison $comparison must have one distinct baseline")
    return (only(baseline_labels), comparison)
end

function comparison_scenarios(comparison)
    return unique(row["scenario"] for row in rows if row["comparison"] == comparison)
end

function sample_stats(comparison, label, scenario, field)
    values = [
        number(row, field) for row in rows
        if row["comparison"] == comparison && row["label"] == label && row["scenario"] == scenario
    ]
    isempty(values) && error("no $field samples for $comparison: $label/$scenario")
    return median_iqr(values)
end

function build_stats(comparison, label, field)
    values_by_build = Dict{String,Float64}()
    for row in rows
        row["comparison"] == comparison && row["label"] == label || continue
        build = row["build"]
        value = number(row, field)
        if haskey(values_by_build, build) && values_by_build[build] != value
            error("$field varies within $comparison: $label build $build")
        end
        values_by_build[build] = value
    end
    isempty(values_by_build) && error("no $field build results for $comparison: $label")
    return median_iqr(collect(values(values_by_build)))
end

function sample_medians_by_build(comparison, label, scenario, field)
    values_by_build = Dict{String,Vector{Float64}}()
    for row in rows
        row["comparison"] == comparison && row["label"] == label && row["scenario"] == scenario || continue
        push!(get!(Vector{Float64}, values_by_build, row["build"]), number(row, field))
    end
    return Dict(build => median(values) for (build, values) in values_by_build)
end

function material_counts(comparison, label, scenario)
    baseline_label = first(comparison_labels(comparison))
    baseline = sample_medians_by_build(comparison, baseline_label, scenario, "total_seconds")
    variant = sample_medians_by_build(comparison, label, scenario, "total_seconds")
    Set(keys(baseline)) == Set(keys(variant)) || error(
        "baseline and candidate build sets differ for $comparison/$scenario"
    )
    builds = sort!(collect(keys(baseline)); by=x -> parse(Int, x))
    improved = count(builds) do build
        threshold = max(0.050, 0.05 * baseline[build])
        variant[build] <= baseline[build] - threshold
    end
    regressed = count(builds) do build
        threshold = max(0.050, 0.05 * baseline[build])
        variant[build] >= baseline[build] + threshold
    end
    return improved, regressed, length(builds)
end

delta(value, baseline) = value - baseline
percent_delta(value, baseline) = iszero(baseline) ? NaN : 100 * delta(value, baseline) / baseline
signed(value; digits=2) = @sprintf("%+.*f", digits, value)
measurement(value, iqr; scale=1.0, digits=2) = @sprintf("%.*f [%.*f]", digits, value * scale, digits, iqr * scale)

columns = [
    "comparison", "label", "scenario", "builds", "samples_per_build",
    "build_s_median", "build_s_iqr", "build_delta_s", "build_delta_pct",
    "cache_mib_median", "cache_mib_iqr", "cache_delta_mib", "cache_delta_pct",
    "import_ms_median", "import_ms_iqr", "import_delta_ms", "import_delta_pct",
    "first_ms_median", "first_ms_iqr", "first_delta_ms", "first_delta_pct",
    "compile_ms_median", "compile_ms_iqr", "recompile_ms_median", "recompile_ms_iqr",
    "total_ms_median", "total_ms_iqr", "total_delta_ms", "total_delta_pct",
    "warm_ms_median", "warm_ms_iqr", "warm_delta_ms", "warm_delta_pct",
    "warm_compile_ms_median", "warm_compile_ms_iqr",
    "warm_recompile_ms_median", "warm_recompile_ms_iqr",
    "total_material_improvements", "total_material_regressions", "compared_builds",
]

open(build_summary_path, "w") do io
    println(io, join((
        "comparison", "label", "scenario", "build", "samples",
        "build_seconds", "cache_bytes",
        "import_ms_median", "first_ms_median", "compile_ms_median",
        "recompile_ms_median", "total_ms_median", "warm_ms_median",
        "warm_compile_ms_median", "warm_recompile_ms_median",
        "baseline_total_ms_median", "total_delta_ms", "total_delta_pct",
        "material_threshold_ms", "material_improvement", "material_regression",
    ), '\t'))
    for comparison in comparisons
        labels = comparison_labels(comparison)
        baseline_label = first(labels)
        for label in labels, scenario in comparison_scenarios(comparison)
            total_by_build = sample_medians_by_build(comparison, label, scenario, "total_seconds")
            isempty(total_by_build) && continue
            baseline_total_by_build = sample_medians_by_build(
                comparison, baseline_label, scenario, "total_seconds"
            )
            for build in sort!(collect(keys(total_by_build)); by=x -> parse(Int, x))
                haskey(baseline_total_by_build, build) || error("baseline has no $scenario build $build")
                matching = [
                    row for row in rows
                    if row["comparison"] == comparison && row["label"] == label &&
                        row["scenario"] == scenario && row["build"] == build
                ]
                med(field) = median(number(row, field) for row in matching)
                total = total_by_build[build]
                baseline_total = baseline_total_by_build[build]
                threshold = max(0.050, 0.05 * baseline_total)
                println(io, join(Any[
                    comparison, label, scenario, build, length(matching),
                    first(matching)["build_seconds"], first(matching)["cache_bytes"],
                    med("import_seconds") * 1e3,
                    med("first_seconds") * 1e3,
                    med("first_compile_seconds") * 1e3,
                    med("first_recompile_seconds") * 1e3,
                    total * 1e3,
                    med("warm_seconds") * 1e3,
                    med("warm_compile_seconds") * 1e3,
                    med("warm_recompile_seconds") * 1e3,
                    baseline_total * 1e3,
                    delta(total, baseline_total) * 1e3,
                    percent_delta(total, baseline_total),
                    threshold * 1e3,
                    total <= baseline_total - threshold,
                    total >= baseline_total + threshold,
                ], '\t'))
            end
        end
    end
end

open(summary_path, "w") do io
    println(io, join(columns, '\t'))
    for comparison in comparisons
        labels = comparison_labels(comparison)
        baseline_label = first(labels)
        for label in labels, scenario in comparison_scenarios(comparison)
            matching = [
                row for row in rows
                if row["comparison"] == comparison && row["label"] == label && row["scenario"] == scenario
            ]
            isempty(matching) && continue

            builds = length(unique(row["build"] for row in matching))
            samples_per_build = unique(
                count(row -> row["build"] == build, matching)
                for build in unique(row["build"] for row in matching)
            )
            length(samples_per_build) == 1 || error("sample count varies by build for $label/$scenario")
            build, build_iqr = build_stats(comparison, label, "build_seconds")
            cache, cache_iqr = build_stats(comparison, label, "cache_bytes")
            import_time, import_iqr = sample_stats(comparison, label, scenario, "import_seconds")
            first_time, first_iqr = sample_stats(comparison, label, scenario, "first_seconds")
            compile_time, compile_iqr = sample_stats(comparison, label, scenario, "first_compile_seconds")
            recompile_time, recompile_iqr = sample_stats(comparison, label, scenario, "first_recompile_seconds")
            total_time, total_iqr = sample_stats(comparison, label, scenario, "total_seconds")
            warm_time, warm_iqr = sample_stats(comparison, label, scenario, "warm_seconds")
            warm_compile, warm_compile_iqr = sample_stats(comparison, label, scenario, "warm_compile_seconds")
            warm_recompile, warm_recompile_iqr = sample_stats(comparison, label, scenario, "warm_recompile_seconds")

            baseline_build = first(build_stats(comparison, baseline_label, "build_seconds"))
            baseline_cache = first(build_stats(comparison, baseline_label, "cache_bytes"))
            baseline_import = first(sample_stats(comparison, baseline_label, scenario, "import_seconds"))
            baseline_first = first(sample_stats(comparison, baseline_label, scenario, "first_seconds"))
            baseline_total = first(sample_stats(comparison, baseline_label, scenario, "total_seconds"))
            baseline_warm = first(sample_stats(comparison, baseline_label, scenario, "warm_seconds"))
            improved, regressed, compared = material_counts(comparison, label, scenario)

            values = Any[
                comparison, label, scenario, builds, only(samples_per_build),
                build, build_iqr, delta(build, baseline_build), percent_delta(build, baseline_build),
                cache / 2.0^20, cache_iqr / 2.0^20, delta(cache, baseline_cache) / 2.0^20, percent_delta(cache, baseline_cache),
                import_time * 1e3, import_iqr * 1e3, delta(import_time, baseline_import) * 1e3, percent_delta(import_time, baseline_import),
                first_time * 1e3, first_iqr * 1e3, delta(first_time, baseline_first) * 1e3, percent_delta(first_time, baseline_first),
                compile_time * 1e3, compile_iqr * 1e3, recompile_time * 1e3, recompile_iqr * 1e3,
                total_time * 1e3, total_iqr * 1e3, delta(total_time, baseline_total) * 1e3, percent_delta(total_time, baseline_total),
                warm_time * 1e3, warm_iqr * 1e3, delta(warm_time, baseline_warm) * 1e3, percent_delta(warm_time, baseline_warm),
                warm_compile * 1e3, warm_compile_iqr * 1e3,
                warm_recompile * 1e3, warm_recompile_iqr * 1e3,
                improved, regressed, compared,
            ]
            println(io, join(values, '\t'))
        end
    end
end

open(markdown_path, "w") do io
    println(io, "# QuantumClifford cold-start summary")
    println(io)
    println(io, "Medians are followed by interquartile ranges in brackets. Each candidate has its own adjacent baseline builds.")
    println(io)
    println(io, "| Comparison | Variant | Scenario | Cache build (s) | Build Δ | Cache (MiB) | Cache Δ | Import (ms) | Import Δ | First task (ms) | First Δ | Compile / recompile (ms) | Total (ms) | Total Δ | Material total Δ builds | Warm task (ms) | Warm compile / recompile (ms) |")
    println(io, "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for comparison in comparisons
        labels = comparison_labels(comparison)
        baseline_label = first(labels)
        for label in labels, scenario in comparison_scenarios(comparison)
            any(row -> row["comparison"] == comparison && row["label"] == label && row["scenario"] == scenario, rows) || continue
            build, build_iqr = build_stats(comparison, label, "build_seconds")
            cache, cache_iqr = build_stats(comparison, label, "cache_bytes")
            import_time, import_iqr = sample_stats(comparison, label, scenario, "import_seconds")
            first_time, first_iqr = sample_stats(comparison, label, scenario, "first_seconds")
            compile_time, compile_iqr = sample_stats(comparison, label, scenario, "first_compile_seconds")
            recompile_time, recompile_iqr = sample_stats(comparison, label, scenario, "first_recompile_seconds")
            total_time, total_iqr = sample_stats(comparison, label, scenario, "total_seconds")
            warm_time, warm_iqr = sample_stats(comparison, label, scenario, "warm_seconds")
            warm_compile, warm_compile_iqr = sample_stats(comparison, label, scenario, "warm_compile_seconds")
            warm_recompile, warm_recompile_iqr = sample_stats(comparison, label, scenario, "warm_recompile_seconds")

            baseline_build = first(build_stats(comparison, baseline_label, "build_seconds"))
            baseline_cache = first(build_stats(comparison, baseline_label, "cache_bytes"))
            baseline_import = first(sample_stats(comparison, baseline_label, scenario, "import_seconds"))
            baseline_first = first(sample_stats(comparison, baseline_label, scenario, "first_seconds"))
            baseline_total = first(sample_stats(comparison, baseline_label, scenario, "total_seconds"))

            build_change = "$(signed(delta(build, baseline_build))) s ($(signed(percent_delta(build, baseline_build)))%)"
            cache_change = "$(signed(delta(cache, baseline_cache) / 2.0^20)) MiB ($(signed(percent_delta(cache, baseline_cache)))%)"
            import_change = "$(signed(delta(import_time, baseline_import) * 1e3)) ms ($(signed(percent_delta(import_time, baseline_import)))%)"
            first_change = "$(signed(delta(first_time, baseline_first) * 1e3)) ms ($(signed(percent_delta(first_time, baseline_first)))%)"
            total_change = "$(signed(delta(total_time, baseline_total) * 1e3)) ms ($(signed(percent_delta(total_time, baseline_total)))%)"
            improved, regressed, compared = material_counts(comparison, label, scenario)

            first_compiler = "$(measurement(compile_time, compile_iqr; scale=1e3)) / $(measurement(recompile_time, recompile_iqr; scale=1e3))"
            warm_compiler = "$(measurement(warm_compile, warm_compile_iqr; scale=1e3)) / $(measurement(warm_recompile, warm_recompile_iqr; scale=1e3))"
            println(io, "| $comparison | $label | $scenario | $(measurement(build, build_iqr)) | $build_change | $(measurement(cache, cache_iqr; scale=1 / 2.0^20)) | $cache_change | $(measurement(import_time, import_iqr; scale=1e3)) | $import_change | $(measurement(first_time, first_iqr; scale=1e3)) | $first_change | $first_compiler | $(measurement(total_time, total_iqr; scale=1e3)) | $total_change | $improved better, $regressed worse / $compared | $(measurement(warm_time, warm_iqr; scale=1e3)) | $warm_compiler |")
        end
    end
end
