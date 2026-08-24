using InteractiveUtils
versioninfo(; verbose=true)

using ParallelTestRunner

const TEST_PROJECTS = Dict(
    "oscar" => normpath(joinpath(@__DIR__, "projects", "oscar")),
    "python_decoders" => normpath(joinpath(@__DIR__, "projects", "python_decoders")),
    "cuda" => normpath(joinpath(@__DIR__, "projects", "cuda")),
    "rocm" => normpath(joinpath(@__DIR__, "projects", "rocm")),
    "opencl" => normpath(joinpath(@__DIR__, "projects", "opencl")),
    "jet" => normpath(joinpath(@__DIR__, "projects", "jet")),
)

const BACKEND_PREFIXES = (
    "KernelAbstractions/CUDA/",
    "KernelAbstractions/ROCm/",
    "KernelAbstractions/OpenCL/",
)
const PYTHON_PREFIX = "general/ecc/python_decoders/"
const OSCAR_PREFIXES = ("general/oscar/", "general/ecc/oscar/")
const OSCAR_ONLY_TESTS = (
    "general/oscar/closure_boxes_tests",
    "general/ecc/oscar/bivaraite_bicycle_via_quotient_ring_tests",
    "general/ecc/oscar/d_dimensional_codes_tests",
    "general/ecc/oscar/dihedral2bga_tests",
    "general/ecc/oscar/generalized_toric_codes_tests",
    "general/ecc/oscar/trivariate_tricycle_tests",
)

const DOWNGRADE_TEST = get(ENV, "QUANTUMSAVORY_DOWNGRADE_TEST", "") == "true"
const OSCAR_SUPPORTED =
    !DOWNGRADE_TEST &&
    !Sys.iswindows() &&
    Sys.ARCH == :x86_64 &&
    VERSION >= v"1.11"
const PYTHON_SUPPORTED = !DOWNGRADE_TEST && !Sys.iswindows()

args = parse_args(isempty(ARGS) ? ["general"] : ARGS)
isempty(args.positionals) && push!(args.positionals, "general")
jet_only =
    isnothing(args.list) &&
    length(args.positionals) == 1 &&
    startswith(only(args.positionals), "jet")

if jet_only
    using Pkg
    Pkg.activate(TEST_PROJECTS["jet"])
    include(joinpath(@__DIR__, "jet_tests.jl"))
    include(joinpath(@__DIR__, "jet_opt_tests.jl"))
else
    testsuite = find_tests(@__DIR__)
    filter!(testsuite) do (name, _)
        endswith(name, "_tests") ||
            any(prefix -> startswith(name, prefix), BACKEND_PREFIXES)
    end

    if !OSCAR_SUPPORTED
        foreach(name -> delete!(testsuite, name), OSCAR_ONLY_TESTS)
    end
    if !PYTHON_SUPPORTED
        filter!(testsuite) do (name, _)
            !startswith(name, PYTHON_PREFIX)
        end
    end
    if DOWNGRADE_TEST
        delete!(testsuite, "general/aqua_tests")
        delete!(testsuite, "general/oscar/doctests_tests")
    end
    if !(Sys.islinux() && Int === Int64)
        delete!(testsuite, "general/bitpack_tests")
    end

    function test_project(name)
        if PYTHON_SUPPORTED && startswith(name, PYTHON_PREFIX)
            return TEST_PROJECTS["python_decoders"]
        elseif OSCAR_SUPPORTED &&
               any(prefix -> startswith(name, prefix), OSCAR_PREFIXES)
            return TEST_PROJECTS["oscar"]
        elseif startswith(name, "KernelAbstractions/CUDA/")
            return TEST_PROJECTS["cuda"]
        elseif startswith(name, "KernelAbstractions/ROCm/")
            return TEST_PROJECTS["rocm"]
        elseif startswith(name, "KernelAbstractions/OpenCL/")
            return TEST_PROJECTS["opencl"]
        end
        return nothing
    end

    function test_worker(name)
        project = test_project(name)
        project === nothing && return nothing

        init_worker_code = if project == TEST_PROJECTS["oscar"]
            quote
                using Pkg
                Pkg.activate($project)
                using QuantumClifford
                import Oscar
            end
        else
            quote
                using Pkg
                Pkg.activate($project)
            end
        end

        return addworker(; init_worker_code)
    end

    using QuantumClifford
    runtests(
        QuantumClifford,
        args;
        testsuite,
        test_worker,
        init_code = :(using QuantumClifford),
    )
end
