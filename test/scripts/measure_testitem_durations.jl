#!/usr/bin/env julia

using Dates
using Printf
using Pkg
using TestItemRunner

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TEST_ROOT = normpath(joinpath(REPO_ROOT, "test"))
const QECCORE_TEST_ROOT = normpath(joinpath(REPO_ROOT, "lib", "QECCore", "test"))

Base.@kwdef struct Config
    group::Symbol = :all
    precompile::Bool = false
    warmup::Bool = false
    top::Int = 20
    outdir::String = normpath(joinpath(TEST_ROOT, "timings"))
end

function print_usage()
    println("""
    Usage:
      julia --project=test test/scripts/measure_testitem_durations.jl [--group=all|general|ecc] [--precompile] [--warmup] [--top=N] [--outdir=PATH]

    Notes:
      - Run from QuantumClifford.jl repository root.
      - This script uses TestItemRunner directly and measures one @testitem at a time.
    """)
end

function parse_args(args)
    group = :all
    precompile = false
    warmup = false
    top = 20
    outdir = normpath(joinpath(TEST_ROOT, "timings"))

    for arg in args
        if arg in ("-h", "--help")
            print_usage()
            exit(0)
        elseif startswith(arg, "--group=")
            v = Symbol(split(arg, "=", limit=2)[2])
            v in (:all, :general, :ecc) || error("Invalid --group value: $v")
            group = v
        elseif arg == "--precompile"
            precompile = true
        elseif arg == "--warmup"
            warmup = true
        elseif startswith(arg, "--top=")
            n = parse(Int, split(arg, "=", limit=2)[2])
            n > 0 || error("--top must be positive.")
            top = n
        elseif startswith(arg, "--outdir=")
            d = split(arg, "=", limit=2)[2]
            outdir = normpath(d)
        else
            error("Unknown argument: $arg")
        end
    end

    return Config(;
        group=group,
        precompile=precompile,
        warmup=warmup,
        top=top,
        outdir=outdir,
    )
end

is_pkg_available(name::String) = !isnothing(Base.find_package(name))
is_oscar_supported() = !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11" && is_pkg_available("Oscar")
is_tesseract_supported() = !Sys.iswindows() && is_pkg_available("PyTesseractDecoder")
is_bitpack_supported() = Base.Sys.islinux() && (Int === Int64)

function is_ecc_item(ti)
    p = normpath(ti.filename)
    return (:ecc in ti.tags) || startswith(p, QECCORE_TEST_ROOT)
end

function has_any_tag(tags, badtags)
    for t in tags
        t in badtags && return true
    end
    return false
end

function include_item(ti, group::Symbol)
    p = normpath(ti.filename)
    in_package_tests = startswith(p, TEST_ROOT) || startswith(p, QECCORE_TEST_ROOT)
    in_package_tests || return false

    excluded = Set{Symbol}([:jet, :aqua, :doctests, :cuda, :rocm, :opencl])
    !is_oscar_supported() && push!(excluded, :oscar_required)
    !is_tesseract_supported() && push!(excluded, :tesseract_required)
    !is_bitpack_supported() && push!(excluded, :bitpack)

    has_any_tag(ti.tags, excluded) && return false

    if group == :general
        return !is_ecc_item(ti)
    elseif group == :ecc
        return is_ecc_item(ti)
    else
        return true
    end
end

function collect_testitems(group::Symbol)
    items = NamedTuple{(:filename, :name, :tags),Tuple{String,String,Vector{Symbol}}}[]
    TestItemRunner.run_tests(REPO_ROOT; filter=ti -> begin
        include_item(ti, group) && push!(items, ti)
        return false
    end)

    sort!(items; by=ti -> (ti.filename, ti.name))
    return items
end

function run_testitem_once(item)
    # Measure a single @testitem by filtering on filename+name.
    elapsed = @elapsed begin
        TestItemRunner.run_tests(REPO_ROOT; filter=ti -> (
            ti.filename == item.filename && ti.name == item.name
        ))
    end
    return elapsed
end

function write_results(cfg::Config, rows)
    mkpath(cfg.outdir)
    stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    csv_path = joinpath(cfg.outdir, "testitem_durations_$(cfg.group)_$(stamp).csv")
    md_path = joinpath(cfg.outdir, "testitem_durations_$(cfg.group)_$(stamp).md")

    open(csv_path, "w") do io
        println(io, "rank,seconds,status,filename,name,tags,error")
        for (i, row) in enumerate(rows)
            tags = join(String.(row.tags), "|")
            filename = replace(row.filename, ',' => ';')
            name = replace(row.name, ',' => ';')
            err = replace(row.error, ',' => ';')
            @printf(io, "%d,%.6f,%s,%s,%s,%s,%s\n", i, row.seconds, row.status, filename, name, tags, err)
        end
    end

    open(md_path, "w") do io
        println(io, "# Test Item Durations ($(cfg.group))")
        println(io)
        timestamp = Dates.format(now(), Dates.DateFormat("yyyy-mm-dd HH:MM:SS"))
        println(io, "- generated: $timestamp")
        println(io, "- measured items: $(length(rows))")
        println(io)
        println(io, "| rank | seconds | status | test item | file | tags |")
        println(io, "|---:|---:|---|---|---|---|")
        for (i, row) in enumerate(rows[1:min(end, cfg.top)])
            tags = isempty(row.tags) ? "" : join(String.(row.tags), ", ")
            file_rel = relpath(row.filename, REPO_ROOT)
            @printf(
                io,
                "| %d | %.3f | %s | %s | `%s` | %s |\n",
                i,
                row.seconds,
                row.status,
                row.name,
                file_rel,
                tags
            )
        end
    end

    return csv_path, md_path
end

function main(args)
    cfg = parse_args(args)

    println("Repository root: $REPO_ROOT")
    println("Group: $(cfg.group)")
    println("Precompile: $(cfg.precompile)")
    println("Warmup: $(cfg.warmup)")
    println("Output directory: $(cfg.outdir)")

    if cfg.precompile
        println("Precompiling test environment...")
        Pkg.instantiate()
        Pkg.precompile()
    end

    items = collect_testitems(cfg.group)
    println("Collected $(length(items)) test items.")

    rows = NamedTuple{(:seconds, :status, :error, :filename, :name, :tags),Tuple{Float64,String,String,String,String,Vector{Symbol}}}[]
    for (i, item) in enumerate(items)
        print(@sprintf("[%3d/%3d] %s ... ", i, length(items), item.name))
        flush(stdout)

        t = 0.0
        status = "ok"
        errtxt = ""
        try
            if cfg.warmup
                run_testitem_once(item)
            end
            t = run_testitem_once(item)
            println(@sprintf("%.3fs", t))
        catch err
            t = NaN
            status = "error"
            errtxt = sprint(showerror, err)
            println("error")
        end

        push!(rows, (
            seconds=t,
            status=status,
            error=errtxt,
            filename=item.filename,
            name=item.name,
            tags=item.tags,
        ))
    end

    sort!(rows; by=r -> (r.status == "ok" ? 0 : 1, -coalesce(r.seconds, -Inf)))
    csv_path, md_path = write_results(cfg, rows)

    println()
    println("Top $(min(cfg.top, length(rows))) slowest test items:")
    ranked_ok = filter(r -> r.status == "ok", rows)
    for (rank, row) in enumerate(ranked_ok[1:min(cfg.top, end)])
        @printf("%3d. %8.3fs  %-60s  (%s)\n", rank, row.seconds, row.name, relpath(row.filename, REPO_ROOT))
    end

    println()
    println("Saved CSV: $csv_path")
    println("Saved report: $md_path")
end

main(ARGS)
