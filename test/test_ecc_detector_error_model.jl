@testitem "ECC detector error model export" tags=[:ecc, :ecc_base] begin
    using QuantumClifford
    using QuantumClifford.ECC

    function zero_based_true_indices(bits)
        return findall(bit -> bit != 0, bits) .- 1
    end

    function y_logical_targets(faults, n, qubit)
        return [i - 1 for i in 1:size(faults, 1) if xor(faults[i, qubit], faults[i, n + qubit])]
    end

    function dem_error_line(probability, detectors, logical_observables)
        targets = vcat(["D$(i)" for i in detectors], ["L$(i)" for i in logical_observables])
        return isempty(targets) ? "error($probability)" : "error($probability) " * join(targets, " ")
    end

    @testset "Steane7 targets match parity checks and faults matrix" begin
        code = Steane7()
        checks = parity_checks(code)
        faults = faults_matrix(code)
        n = code_n(code)

        dem = detector_error_model(code; px=1e-3, py=2e-3, pz=3e-3)
        lines = split(chomp(dem), '\n')
        first_error_line = code_s(code) + size(faults, 1) + 1

        @test dem isa String
        @test lines[1:code_s(code)] == ["detector D$(i)" for i in 0:code_s(code) - 1]
        @test lines[code_s(code) + 1:code_s(code) + size(faults, 1)] ==
              ["logical_observable L$(i)" for i in 0:size(faults, 1) - 1]
        @test length(lines) == code_s(code) + size(faults, 1) + 3n

        @test lines[first_error_line] == dem_error_line(
            0.001,
            zero_based_true_indices(comm(checks, single_x(n, 1))),
            zero_based_true_indices(@view faults[:, 1]),
        )

        @test lines[first_error_line + 1] == dem_error_line(
            0.002,
            zero_based_true_indices(comm(checks, single_y(n, 1))),
            y_logical_targets(faults, n, 1),
        )

        @test lines[first_error_line + 2] == dem_error_line(
            0.003,
            zero_based_true_indices(comm(checks, single_z(n, 1))),
            zero_based_true_indices(@view faults[:, n + 1]),
        )
    end

    @testset "writer emits deterministic Stim DEM text" begin
        dem_text = detector_error_model(Steane7(); px=1e-3, pz=1e-3)

        io = IOBuffer()
        write_detector_error_model(io, dem_text)
        text = String(take!(io))

        io_again = IOBuffer()
        write_detector_error_model(io_again, dem_text)
        @test text == String(take!(io_again))
        @test text == dem_text

        lines = split(chomp(text), '\n')
        code = Steane7()
        detector_count = code_s(code)
        observable_count = size(faults_matrix(code), 1)
        @test lines[1:detector_count] == ["detector D$(i)" for i in 0:detector_count - 1]
        @test lines[detector_count + 1:detector_count + observable_count] ==
              ["logical_observable L$(i)" for i in 0:observable_count - 1]
        @test all(startswith(line, "error(") for line in lines[detector_count + observable_count + 1:end])
        @test endswith(text, "\n")

        path = tempname() * ".dem"
        write_detector_error_model(path, dem_text)
        @test read(path, String) == text
    end

    @testset "probabilities are validated and zero probabilities are omitted" begin
        dem = detector_error_model(Steane7(); px=0, py=0, pz=0)
        @test !occursin("error(", dem)
        @test_throws "`px`, `py`, and `pz` must be between 0 and 1." detector_error_model(
            Steane7();
            px=-0.1,
        )
        @test_throws "`px`, `py`, and `pz` must be between 0 and 1." detector_error_model(
            Steane7();
            py=1.1,
        )
    end
end
