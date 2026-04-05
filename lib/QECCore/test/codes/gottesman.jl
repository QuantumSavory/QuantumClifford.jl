@testitem "Gottesman codes should correct all single-qubit errors" begin
    using QECCore
    using Test
    import QuantumClifford: single_x, single_y, single_z, comm, nqubits
    import QuantumClifford.ECC: parity_checks

    @testset "Gottesman codes should correct all single-qubit errors" begin
        for j in 3:12
            H = parity_checks(Gottesman(j))
            syndromes = Set([]) # the set automatically removes repeated entries
            for error_type in (single_x, single_y, single_z)
                for bit_index in 1:nqubits(H)
                    syndrome = comm(H, error_type(nqubits(H), bit_index))
                    @test any(==(0x1), syndrome) # checking the syndrome is not trivially zero
                    push!(syndromes, syndrome)
                end
            end
            @test length(syndromes) == 3*nqubits(H)
        end
    end
end
