@testitem "Subsystem Codes: BBS, SHP, and SHYPS" tags=[:ecc] begin
    using QuantumClifford.ECC
    using QuantumClifford
    using Nemo: GF, matrix, Nemo
    using QECCore

    @testset "Simplex Classical Code" begin
        # Simplex(r) = dual of Hamming(r), parameters [2^r-1, r, 2^(r-1)]
        for r in 2:5
            s = QECCore.Simplex(r)
            n = 2^r - 1
            k = r
            d = 2^(r-1)
            @test code_n(s) == n
            @test code_k(s) == k
            @test distance(s) == d

            H = Matrix{Int}(QECCore.parity_matrix(s))
            @test size(H) == (n - k, n)

            # every PCM row is a nonzero codeword of the dual (Hamming) code
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

        @test_throws ArgumentError QECCore.Simplex(1)
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

        # Stabilizers must commute with all gauge generators (key subsystem code property)
        for s in stab
            for g in gg
                @test comm(s, g) == 0x0
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

            # gauge qubits and stabilizer count
            @test code_g(shyps) == 1
            @test length(parity_checks(shyps)) == 4

            # Stabilizers commute with each other
            stab = parity_checks(shyps)
            for i in 1:length(stab)
                for j in i+1:length(stab)
                    @test comm(stab[i], stab[j]) == 0x0
                end
            end

            # Stabilizers commute with all gauge generators
            for s in stab
                for g in gg
                    @test comm(s, g) == 0x0
                end
            end

            # SHYPS is QLDPC: gauge weight is at most r+1 (= 3 for r=2)
            # rows of the simplex PCM are Hamming codewords; the Gaussian basis gives weight <= r+1
            for g in gg
                @test count(i -> g[i] != (false, false), 1:g.nqubits) <= 3
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

            # gauge qubits and stabilizer count
            @test code_g(shyps) == 16
            @test length(parity_checks(shyps)) == 24

            stab = parity_checks(shyps)
            for i in 1:length(stab)
                for j in i+1:length(stab)
                    @test comm(stab[i], stab[j]) == 0x0
                end
            end

            # Stabilizers commute with all gauge generators
            for s in stab
                for g in gg
                    @test comm(s, g) == 0x0
                end
            end

            # SHYPS is QLDPC: gauge weight is at most r+1 (= 4 for r=3)
            for g in gg
                @test count(i -> g[i] != (false, false), 1:g.nqubits) <= 4
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
