@testitem "Subsystem Codes: BBS, SHP, and SHYPS" tags=[:ecc] begin
    using QuantumClifford.ECC
    using QuantumClifford
    using Nemo: GF, matrix, Nemo
    using QECCore: Simplex, parity_matrix, Hamming

    @testset "Simplex Classical Code" begin
        # Simplex(r) = dual of Hamming(r), parameters [2^r-1, r, 2^(r-1)]
        for r in 2:5
            s = Simplex(r)
            n = 2^r - 1
            k = r
            d = 2^(r-1)
            @test code_n(s) == n
            @test code_k(s) == k
            @test distance(s) == d

            H = Matrix{Int}(parity_matrix(s))
            @test size(H) == (n - k, n)

            # Every row of the PCM should have nonzero weight
            for i in 1:size(H, 1)
                @test sum(H[i, :]) > 0
            end

            # H has full rank (n-k) over GF(2)
            F = GF(2)
            @test Nemo.rank(matrix(F, H)) == n - k

            # Right nullspace of H has dimension k (the information dimension of the simplex code)
            K = Nemo.kernel(transpose(matrix(F, H)))
            @test Nemo.nrows(K) == k
        end

        @test_throws ArgumentError Simplex(1)
    end

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
    end

    @testset "SHP Code from Hamming [7,4,3]" begin
        # As given in arXiv:2002.06257 Section 3.3
        H = [1 1 0 1 1 0 0;
             1 0 1 1 0 1 0;
             0 1 1 1 0 0 1]

        shp = SubsystemHypergraphProduct(H, H)

        # N = n1 * n2 = 7 * 7 = 49
        @test code_n(shp) == 49

        # K = nullity(H1) * nullity(H2) = 4 * 4 = 16
        @test code_k(shp) == 16

        # Test CSS
        @test iscss(shp)

        # S_X = kron(H, nullspace(H)): 3×4 = 12 rows
        # S_Z = kron(nullspace(H), H): 4×3 = 12 rows
        @test length(parity_checks(shp)) == 24
        @test size(parity_matrix_x(shp), 1) == 12
        @test size(parity_matrix_z(shp), 1) == 12

        # Gauge generators exist
        gg = gauge_generators(shp)
        @test length(gg) > 0

        # Gauge qubits: g = n - k - rank(stabilizer)
        @test code_g(shp) == 49 - 16 - 24  # = 9

        # Stabilizers must commute with each other
        stab = parity_checks(shp)
        for i in 1:length(stab)
            for j in i+1:length(stab)
                @test comm(stab[i], stab[j]) == 0x0
            end
        end
    end

    @testset "SHYPS Codes" begin
        # SHYPS(r) = SHP(H_simplex(r), H_simplex(r))
        # Parameters: [(2^r-1)^2, r^2, 2^(r-1)]

        @testset "SHYPS(2)" begin
            shyps = SubsystemHypergraphProductSimplex(2)
            # [9, 4, 2]
            @test code_n(shyps) == 9
            @test code_k(shyps) == 4
            @test distance(shyps) == 2
            @test iscss(shyps)

            # Gauge generators exist
            gg = gauge_generators(shyps)
            @test length(gg) > 0

            # n = k + g + s
            g = code_g(shyps)
            s_rank = code_n(shyps) - code_k(shyps) - g
            @test code_n(shyps) == code_k(shyps) + g + s_rank

            # Stabilizers commute
            stab = parity_checks(shyps)
            for i in 1:length(stab)
                for j in i+1:length(stab)
                    @test comm(stab[i], stab[j]) == 0x0
                end
            end
        end

        @testset "SHYPS(3)" begin
            shyps = SubsystemHypergraphProductSimplex(3)
            # [49, 9, 4]
            @test code_n(shyps) == 49
            @test code_k(shyps) == 9
            @test distance(shyps) == 4
            @test iscss(shyps)

            gg = gauge_generators(shyps)
            @test length(gg) > 0

            g = code_g(shyps)
            s_rank = code_n(shyps) - code_k(shyps) - g
            @test code_n(shyps) == code_k(shyps) + g + s_rank

            stab = parity_checks(shyps)
            for i in 1:length(stab)
                for j in i+1:length(stab)
                    @test comm(stab[i], stab[j]) == 0x0
                end
            end
        end

        @test_throws ArgumentError SubsystemHypergraphProductSimplex(1)
    end

    @testset "Gauge generators do not all commute (non-abelian gauge group)" begin
        # For a true subsystem code, the gauge group is non-abelian
        H = [1 1 0 1 1 0 0;
             1 0 1 1 0 1 0;
             0 1 1 1 0 0 1]
        shp = SubsystemHypergraphProduct(H, H)
        gg = gauge_generators(shp)
        found_noncommuting = false
        for i in 1:length(gg)
            for j in i+1:length(gg)
                if comm(gg[i], gg[j]) != 0x0
                    found_noncommuting = true
                    break
                end
            end
            found_noncommuting && break
        end
        @test found_noncommuting
    end
end
