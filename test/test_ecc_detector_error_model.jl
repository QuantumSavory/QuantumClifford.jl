@testitem "ECC detector error model export" tags=[:ecc, :ecc_base] begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    function expected_detectors(code, error)
        [i - 1 for i in findall(!iszero, comm(error, parity_checks(code)))]
    end

    function expected_observables(code, q, pauli)
        n = code_n(code)
        fm = faults_matrix(code)
        bits =
            if pauli === :X
                @view fm[:, q]
            elseif pauli === :Z
                @view fm[:, n + q]
            else
                xor.(@view(fm[:, q]), @view(fm[:, n + q]))
            end
        [i - 1 for i in findall(!iszero, bits)]
    end

    code = Bitflip3()
    dem = detector_error_model(code; px=0.1, py=0.2, pz=0.3)

    @test dem.num_detectors == code_s(code)
    @test dem.num_observables == 2code_k(code)
    @test length(dem) == 3code_n(code)

    first_x = dem.errors[1]
    @test first_x.probability == 0.1
    @test first_x.detectors == expected_detectors(code, single_x(code_n(code), 1))
    @test first_x.observables == expected_observables(code, 1, :X)

    third_y = dem.errors[8]
    @test third_y.probability == 0.2
    @test third_y.detectors == expected_detectors(code, single_y(code_n(code), 3))
    @test third_y.observables == expected_observables(code, 3, :Y)

    io = IOBuffer()
    write_detector_error_model(io, dem)
    dem_text = String(take!(io))

    @test startswith(dem_text, "detector D0\ndetector D1\nlogical_observable L0\nlogical_observable L1\n")
    @test occursin("error(0.1) D0\n", dem_text)
    @test occursin("error(0.2) D1 L0 L1\n", dem_text)
    @test occursin("error(0.3) L0\n", dem_text)

    sparse_dem = detector_error_model(Steane7(); px=1e-3, pz=1e-3)
    @test all(term.probability == 1e-3 for term in sparse_dem.errors)
    @test length(sparse_dem) <= 2code_n(Steane7())

    @test_throws DomainError detector_error_model(code; px=-0.1)
    @test_throws DomainError detector_error_model(code; py=1.1)

    mktempdir() do dir
        path = joinpath(dir, "bitflip.dem")
        write_detector_error_model(path, dem)
        @test read(path, String) == dem_text

        python = Sys.which("python")
        python === nothing && (python = Sys.which("python3"))
        if python !== nothing
            script = """
import sys
try:
    import stim
except Exception:
    sys.exit(0)
stim.DetectorErrorModel(open(sys.argv[1], encoding='utf-8').read())
"""
            run(`$python -c $script $path`)
        end
    end
end

@testitem "Detector error model coverage smoke" begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    # Keep this small smoke test in the default CI set so coverage exercises the
    # public DEM API while the heavier ECC assertions stay under ECC_TEST=base.
    code = Bitflip3()
    empty_dem = detector_error_model(code)
    @test isempty(empty_dem)
    @test length(empty_dem) == 0
    @test empty_dem.num_detectors == code_s(code)
    @test empty_dem.num_observables == 2code_k(code)

    dem = detector_error_model(code; px=0.125, py=0.25, pz=0.5)
    @test !isempty(dem)
    @test length(dem) == 3code_n(code)
    @test all(term.probability in (0.125, 0.25, 0.5) for term in dem.errors)
    @test any(!isempty(term.detectors) for term in dem.errors)
    @test any(!isempty(term.observables) for term in dem.errors)
    @test_throws DomainError detector_error_model(code; pz=1.01)

    io = IOBuffer()
    @test write_detector_error_model(io, dem) === nothing
    text = String(take!(io))
    @test occursin("detector D0", text)
    @test occursin("logical_observable L0", text)
    @test occursin("error(0.125)", text)
    @test occursin("error(0.25)", text)
    @test occursin("error(0.5)", text)

    code_io = IOBuffer()
    @test write_detector_error_model(code_io, code; px=0.125) === nothing
    @test occursin("error(0.125)", String(take!(code_io)))

    mktempdir() do dir
        dem_path = joinpath(dir, "bitflip.dem")
        @test write_detector_error_model(dem_path, dem) === nothing
        @test read(dem_path, String) == text

        code_path = joinpath(dir, "bitflip-from-code.dem")
        @test write_detector_error_model(code_path, code; pz=0.5) === nothing
        @test occursin("error(0.5)", read(code_path, String))
    end
end
