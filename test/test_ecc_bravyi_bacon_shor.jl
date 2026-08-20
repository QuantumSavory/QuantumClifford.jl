@testitem "Bravyi-Bacon-Shor Subsystem Code" tags=[:ecc] begin
    using QuantumClifford.ECC
    using QuantumClifford
    using Nemo: GF, matrix, Nemo

    @testset "BBS Code from Hamming [7,4,3]" begin
        # As given in arXiv:2002.06257 Section 2.2
        A = [0 0 1 0 0 1 1;
             0 1 0 1 0 1 0;
             1 0 0 0 1 1 0;
             0 1 0 0 1 0 1;
             0 0 1 1 1 0 0;
             1 1 1 0 0 0 0;
             1 0 0 1 0 0 1]

        bbs = BravyiBaconShor(A)

        # N = number of nonzero entries in A
        @test code_n(bbs) == 21

        # K = rank(A) over GF(2)
        @test code_k(bbs) == 4

        # Test CSS property
        @test iscss(bbs)

        # Stabilizer generators: left kernel of A has dim 3 → 3 X stabilizers
        # Right kernel of A has dim 3 → 3 Z stabilizers
        @test length(parity_checks(bbs)) == 6
        @test size(parity_matrix_x(bbs), 1) == 3
        @test size(parity_matrix_z(bbs), 1) == 3

        # Gauge generators exist
        gg = gauge_generators(bbs)
        @test length(gg) > 0

        # Gauge qubits: g = n - k - s
        @test code_g(bbs) == 21 - 4 - 6  # = 11

        # Stabilizers must commute with each other
        stab = parity_checks(bbs)
        for i in 1:length(stab)
            for j in i+1:length(stab)
                @test comm(stab[i], stab[j]) == 0x0
            end
        end

        # Stabilizers must commute with all gauge generators (key subsystem code property)
        for s in stab
            for g in gg
                @test comm(s, g) == 0x0
            end
        end

        # BBS gauge generators are all weight 2 (each pairs exactly two qubits)
        for g in gg
            @test count(i -> g[i] != (false, false), 1:g.nqubits) == 2
        end
    end
end
