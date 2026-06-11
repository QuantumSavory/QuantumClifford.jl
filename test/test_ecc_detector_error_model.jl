@testitem "ECC detector error model export" tags=[:ecc, :ecc_base] begin
    using Test
    using QuantumClifford
    using QuantumClifford: stab_to_gf2, comm, single_x, single_y, single_z
    using QuantumClifford.ECC

    """Recompute the expected mechanisms from first principles: detectors from the
    anticommutation of the fault with each parity check row (acceptance criterion:
    "detector targets match the stabilizers flipped by the corresponding Pauli error")
    and observables literally as `faults_matrix(code) * stab_to_gf2(fault)` (criterion:
    "logical targets match `faults_matrix(code)`")."""
    function expected_mechanisms(H, O, px, py, pz)
        s, n = size(H)
        expected = Vector{Tuple{Float64,Vector{Int},Vector{Int}}}()
        for q in 1:n
            for (p, fault) in ((px, single_x(n, q)), (py, single_y(n, q)), (pz, single_z(n, q)))
                iszero(p) && continue
                dets = Int[i-1 for i in 1:s if comm(fault, H, i) == 0x1]
                flips = O * stab_to_gf2(fault)
                obs = Int[j-1 for j in 1:size(O, 1) if isodd(flips[j])]
                isempty(dets) && isempty(obs) && continue
                push!(expected, (Float64(p), dets, obs))
            end
        end
        return expected
    end

    dem_text(dem) = sprint(write_detector_error_model, dem)

    # Codes with full-rank parity checks (`faults_matrix` warns on redundant rows;
    # such codes still export fine, but we keep the test output clean).
    testable_codes() = [Steane7(), Shor9(), Perfect5(), Cleve8(), Bitflip3(), Gottesman(3)]

    @testset "detector targets match flipped stabilizers; logical targets match faults_matrix" begin
        for code in testable_codes()
            H = parity_checks(code)
            O = faults_matrix(code)
            px, py, pz = 0.001, 0.002, 0.003
            dem = detector_error_model(code; px, py, pz)
            @test dem.num_detectors == size(H, 1)
            @test dem.num_observables == size(O, 1) == 2 * code_k(code)
            expected = expected_mechanisms(H, O, px, py, pz)
            @test length(dem.errors) == length(expected)
            for (e, (p, dets, obs)) in zip(dem.errors, expected)
                @test e.probability == p
                @test e.detectors == dets       # also pins the ascending order
                @test e.observables == obs
                @test issorted(e.detectors) && issorted(e.observables)
            end
        end
    end

    @testset "deterministic output" begin
        for code in testable_codes()
            a = dem_text(detector_error_model(code; px=0.001, py=0.002, pz=0.003))
            b = dem_text(detector_error_model(code; px=0.001, py=0.002, pz=0.003))
            @test a == b
            # the `Stabilizer` and the code-object methods agree
            @test detector_error_model(parity_checks(code); px=0.001, py=0.002, pz=0.003) ==
                  detector_error_model(code; px=0.001, py=0.002, pz=0.003)
        end
    end

    @testset "golden file: Steane7 (CSS), the issue's flagship call" begin
        dem = detector_error_model(Steane7(); px=1e-3, py=0.0, pz=1e-3)
        @test dem_text(dem) == """detector D0
detector D1
detector D2
detector D3
detector D4
detector D5
logical_observable L0
logical_observable L1
error(0.001) D5
error(0.001) D2
error(0.001) D4 L1
error(0.001) D1
error(0.001) D4 D5
error(0.001) D1 D2 L0
error(0.001) D3 L1
error(0.001) D0
error(0.001) D3 D5
error(0.001) D0 D2 L0
error(0.001) D3 D4 L1
error(0.001) D0 D1 L0
error(0.001) D3 D4 D5
error(0.001) D0 D1 D2
"""
    end

    @testset "golden file: Perfect5 (non-CSS)" begin
        dem = detector_error_model(Perfect5(); px=0.001, py=0.002, pz=0.003)
        @test dem_text(dem) == """detector D0
detector D1
detector D2
detector D3
logical_observable L0
logical_observable L1
error(0.001) D3 L0 L1
error(0.002) D0 D2 D3 L0 L1
error(0.003) D0 D2
error(0.001) D0 L1
error(0.002) D0 D1 D3 L1
error(0.003) D1 D3
error(0.001) D0 D1 L1
error(0.002) D0 D1 D2 L1
error(0.003) D2
error(0.001) D1 D2 L0 L1
error(0.002) D0 D1 D2 D3 L0 L1
error(0.003) D0 D3
error(0.001) D2 D3 L1
error(0.002) D1 D2 D3 L0 L1
error(0.003) D1 L0
"""
    end

    @testset "writing to IO and to a file path agree" begin
        dem = detector_error_model(Steane7(); px=1e-3, pz=1e-3)
        mktempdir() do dir
            path = joinpath(dir, "steane7.dem")
            write_detector_error_model(path, dem)
            @test read(path, String) == dem_text(dem)
        end
    end

    @testset "zero probabilities are skipped" begin
        dem = detector_error_model(Steane7(); px=1e-3)
        @test length(dem.errors) == 7                       # one X mechanism per qubit
        @test all(e.probability == 1e-3 for e in dem.errors)
        dem0 = detector_error_model(Steane7())
        @test isempty(dem0.errors)                          # declarations only
        @test dem_text(dem0) == join(vcat(["detector D$(i)" for i in 0:5],
                                          ["logical_observable L$(i)" for i in 0:1]), '\n') * '\n'
    end

    @testset "undetectable logical faults are kept; duplicate mechanisms are not merged" begin
        # All three Z faults on the bit-flip code commute with every check but flip the
        # logical X observable, producing three identical `error(p) L0` mechanisms --
        # Stim's format explicitly allows repeated mechanisms with the same targets.
        dem = detector_error_model(Bitflip3(); pz=0.05)
        @test length(dem.errors) == 3
        @test all(e.detectors == Int[] && e.observables == [0] for e in dem.errors)
    end

    @testset "probability validation" begin
        @test_throws ArgumentError detector_error_model(Steane7(); px=-0.1)
        @test_throws ArgumentError detector_error_model(Steane7(); py=1.5)
        @test_throws ArgumentError detector_error_model(Steane7(); pz=Inf)
    end
end

@testitem "ECC detector error model -- Stim parser smoke test" tags=[:ecc, :ecc_base] begin
    # Optional and off by default so that the Julia test suite does not require Python:
    # enable by setting the environment variable QC_TEST_STIM=true with a `python3` on
    # PATH (override with QC_TEST_STIM_PYTHON) that can `import stim`.
    using Test
    using QuantumClifford.ECC
    if get(ENV, "QC_TEST_STIM", "") == "true"
        pyexe = get(ENV, "QC_TEST_STIM_PYTHON", "python3")
        mktempdir() do dir
            path = joinpath(dir, "steane7.dem")
            write_detector_error_model(path, detector_error_model(Steane7(); px=1e-3, py=0.0, pz=1e-3))
            script = """
                     import stim
                     dem = stim.DetectorErrorModel.from_file(r'$(path)')
                     errors = [inst for inst in dem if inst.type == 'error']
                     assert all(inst.args_copy() == [0.001] for inst in errors)
                     print(dem.num_detectors, dem.num_observables, len(errors))
                     """
            out = strip(read(`$(pyexe) -c $(script)`, String))
            @test out == "6 2 14"
        end
    else
        @info "Skipping the optional Stim parser smoke test -- set ENV QC_TEST_STIM=true (with python3 + stim available) to enable it."
    end
end
