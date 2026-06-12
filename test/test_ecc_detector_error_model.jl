@testitem "ECC detector error model export" tags=[:ecc, :ecc_base] begin
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: DetectorErrorModelError

    function zero_based_true_indices(bits)
        return findall(bit -> bit != 0, bits) .- 1
    end

    function y_logical_targets(faults, n, qubit)
        return [i - 1 for i in 1:size(faults, 1) if xor(faults[i, qubit], faults[i, n + qubit])]
    end

    @testset "Steane7 targets match parity checks and faults matrix" begin
        code = Steane7()
        checks = parity_checks(code)
        faults = faults_matrix(code)
        n = code_n(code)

        dem = detector_error_model(code; px=1e-3, py=2e-3, pz=3e-3)

        @test dem.num_detectors == code_s(code)
        @test dem.num_observables == size(faults, 1)
        @test length(dem.errors) == 3n

        x_error = dem.errors[1]
        @test x_error isa DetectorErrorModelError
        @test x_error.probability == 1e-3
        @test x_error.detectors == zero_based_true_indices(comm(checks, single_x(n, 1)))
        @test x_error.logical_observables == zero_based_true_indices(@view faults[:, 1])

        y_error = dem.errors[2]
        @test y_error.probability == 2e-3
        @test y_error.detectors == zero_based_true_indices(comm(checks, single_y(n, 1)))
        @test y_error.logical_observables == y_logical_targets(faults, n, 1)

        z_error = dem.errors[3]
        @test z_error.probability == 3e-3
        @test z_error.detectors == zero_based_true_indices(comm(checks, single_z(n, 1)))
        @test z_error.logical_observables == zero_based_true_indices(@view faults[:, n + 1])
    end

    @testset "writer emits deterministic Stim DEM text" begin
        dem = detector_error_model(Steane7(); px=1e-3, pz=1e-3)

        io = IOBuffer()
        write_detector_error_model(io, dem)
        text = String(take!(io))

        io_again = IOBuffer()
        write_detector_error_model(io_again, dem)
        @test text == String(take!(io_again))

        lines = split(chomp(text), '\n')
        @test lines[1:dem.num_detectors] == ["detector D$(i)" for i in 0:dem.num_detectors - 1]
        @test lines[dem.num_detectors + 1:dem.num_detectors + dem.num_observables] ==
              ["logical_observable L$(i)" for i in 0:dem.num_observables - 1]
        @test all(startswith(line, "error(") for line in lines[dem.num_detectors + dem.num_observables + 1:end])
        @test endswith(text, "\n")

        path = tempname() * ".dem"
        write_detector_error_model(path, dem)
        @test read(path, String) == text

        python = Sys.which("python3")
        if python !== nothing
            has_stim = readchomp(
                `$python -c "import importlib.util; print(importlib.util.find_spec('stim') is not None)"`,
            )
            if has_stim == "True"
                @test success(
                    `$python -c "import pathlib, stim, sys; stim.DetectorErrorModel(pathlib.Path(sys.argv[1]).read_text())" $path`,
                )
            else
                @info "Skipping Stim parse smoke test because Python package stim is not installed."
            end
        end
    end

    @testset "probabilities are validated and zero probabilities are omitted" begin
        @test isempty(detector_error_model(Steane7(); px=0, py=0, pz=0).errors)
        @test_throws ArgumentError detector_error_model(Steane7(); px=-0.1)
        @test_throws ArgumentError detector_error_model(Steane7(); py=1.1)
    end
end
