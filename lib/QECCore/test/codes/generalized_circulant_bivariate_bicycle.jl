@testitem "CirculantBivariateBicycle" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using QECCore
        using QECCore.LinearAlgebra
        using QuantumClifford
        using QuantumClifford: stab_looks_good
        using QuantumClifford: stab_to_gf2
        using QuantumClifford.ECC
        using Oscar

        @testset "Generalized CirculantBivariateBicycle codes from [bravyi2024high](@cite)" begin
            test_cases = [
            (
                name = "BB [[108, 8, 10]]",
                l = 9,
                m = 6,
                A = [(:x, 3), (:y, 1), (:y, 2)], # A = x³ + y + y²
                B = [(:y, 3), (:x, 1), (:x, 2)], # B = y³ + x + x²
                expected_n = 108,
                expected_k = 8
            ),
            (
                name = "BB [[90, 8, 10]]",
                l = 15,
                m = 3,
                A = [(:x, 9), (:y, 1), (:y, 2)], # A = x⁹ + y + y²
                B = [(:y, 0), (:x, 2), (:x, 7)], # B = 1 + x² + x⁷
                expected_n = 90,
                expected_k = 8
            ),
            (
                name = "BB [[288, 12, 18]]",
                l = 12,
                m = 12,
                A = [(:x, 3), (:y, 2), (:y, 7)], # A = x³ + y² + y⁷
                B = [(:y, 3), (:x, 1), (:x, 2)], # B = y³ + x + x²
                expected_n = 288,
                expected_k = 12
            ),
            (
                name = "BB [[144, 12, 12]]",
                l = 12,
                m = 6,
                A = [(:x, 3), (:y, 1), (:y, 2)], # A = x³ + y + y²
                B = [(:y, 3), (:x, 1), (:x, 2)], # B = y³ + x + x²
                expected_n = 144,
                expected_k = 12
            ),
            (
                name = "BB [[72, 12, 6]]",
                l = 6,
                m = 6,
                A = [(:x, 3), (:y, 1), (:y, 2)], # A = x³ + y + y²
                B = [(:y, 3), (:x, 1), (:x, 2)], # B = y³ + x + x²
                expected_n = 72,
                expected_k = 12
            ),
            (
                name = "BB [[360, 12, ≤24]]",
                l = 30,
                m = 6,
                A = [(:x, 9), (:y, 1), (:y, 2)],  # A = x⁹ + y + y²
                B = [(:y, 3), (:x, 25), (:x, 26)],# B = y³ + x²⁵ + x²⁶
                expected_n = 360,
                expected_k = 12
            ),
            (
                name = "BB [[756, 16, ≤34]]",
                l = 21,
                m = 18,
                A = [(:x, 3), (:y, 10), (:y, 17)], # A = x³ + y¹⁰ + y¹⁷
                B = [(:y, 5), (:x, 3), (:x, 19)],  # B = y⁵ + x³ + x¹⁹
                expected_n = 756,
                expected_k = 16
            ),
        ]
            @testset "$(case.name): l=$(case.l), m=$(case.m)" for case in test_cases
                c = GeneralizedCirculantBivariateBicycle(case.l, case.m, case.A, case.B)
                stab = parity_checks(c)
                n, k = code_n(c), code_k(c)
                nₛ, kₛ = code_n(stab), code_k(stab)
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k && computed_rank == nₛ - kₛ
                @test n == nₛ == case.expected_n && k == kₛ == case.expected_k
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test bivariate_bicycle_code_k(c) == case.expected_k == code_k(c)
                # Both A and B are matrices where each row and column has exactly three non-zero entries when using 3-term polynomial representation.
                Hx = parity_matrix_x(c)
                n = size(Hx,2)÷2
                A = Hx[:,1:n]
                B = Hx[:,n+1:end]
                @test all(sum(A, dims=2) .== 3)
                @test all(sum(B, dims=2) .== 3)
                @test all(sum(A, dims=1) .== 3)
                @test all(sum(B, dims=1) .== 3)
                @test A*B == B*A
            end
        end

        @testset "Verify number of logical qubits `k` from Table 1: [berthusen2024toward](@cite)" begin
            test_cases = [
            (
                name = "BB [[72 , 8,6]]",
                l = 12,
                m = 3,
                A = [(:x, 9), (:y, 1), (:y, 2)],
                B = [(:y, 0), (:x, 1), (:x, 11)],
                expected_n = 72,
                expected_k = 8
            ),
            (
                name = "[[90 , 8,6]]",
                l = 9,
                m = 5,
                A = [(:x, 8), (:y, 4), (:y, 1)],
                B = [(:y, 5), (:x, 8), (:x, 7)],
                expected_n = 90,
                expected_k = 8
            ),
            (
                name = "BB [[120, 8, 8]]",
                l = 12,
                m = 5,
                A = [(:x, 10), (:y, 4), (:y, 1)],
                B = [(:y, 0), (:x, 1), (:x, 2)],
                expected_n = 120,
                expected_k = 8
            ),
            (
                name = "BB [[150, 8,8]]",
                l = 15,
                m = 5,
                A = [(:x, 5), (:y, 2), (:y, 3)],
                B = [(:y, 2), (:x, 7), (:x, 6)],
                expected_n = 150,
                expected_k = 8
            ),
            (
                name = "BB [[196,12,8]]",
                l = 14,
                m = 7,
                A = [(:x, 5), (:y, 6), (:y, 5)],
                B = [(:y, 0), (:x, 4), (:x, 13)],
                expected_n = 196,
                expected_k = 12
            )
        ]
            @testset "$(case.name): a=$(case.l), b=$(case.m)" for case in test_cases
                c = GeneralizedCirculantBivariateBicycle(case.l, case.m, case.A, case.B)
                stab = parity_checks(c)
                n, k = code_n(c), code_k(c)
                nₛ, kₛ = code_n(stab), code_k(stab)
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k && computed_rank == nₛ - kₛ
                @test n == nₛ == case.expected_n && k == kₛ == case.expected_k
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test bivariate_bicycle_code_k(c) == case.expected_k == code_k(c)
                # Both A and B are matrices where each row and column has exactly three non-zero entries when using 3-term polynomial representation.
                Hx = parity_matrix_x(c)
                n = size(Hx,2)÷2
                A = Hx[:,1:n]
                B = Hx[:,n+1:end]
                @test all(sum(A, dims=2) .== 3)
                @test all(sum(B, dims=2) .== 3)
                @test all(sum(A, dims=1) .== 3)
                @test all(sum(B, dims=1) .== 3)
                @test A*B == B*A
            end
        end

        @testset "Verify number of logical qubits `k` from Table 1: wang2024coprime" begin
            # Test cases from [wang2024coprime](@cite)
            test_cases = [
            (
                name = "BB [[54, 8, 6]]",
                l = 3,
                m = 9,
                A = [(:x, 0), (:y, 2), (:y, 4)],
                B = [(:y, 3), (:x, 1), (:x, 2)],
                expected_n = 54,
                expected_k = 8
            ),
            (
                name = "BB [[98, 6, 12]]",
                l = 7,
                m = 7,
                A = [(:x, 3), (:y, 5), (:y, 6)],
                B = [(:y, 2), (:x, 3), (:x, 5)],
                expected_n = 98,
                expected_k = 6
            ),
            (
                name = "BB [[126, 8, 10]]",
                l = 3,
                m = 21,
                A = [(:x, 0), (:y, 2), (:y, 10)],
                B = [(:y, 3), (:x, 1), (:x, 2)],
                expected_n = 126,
                expected_k = 8
            ),
            (
                name = "BB [[150, 16, 8]]",
                l = 5,
                m = 15,
                A = [(:x, 0), (:y, 6), (:y, 8)],
                B = [(:y, 5), (:x, 1), (:x, 4)],
                expected_n = 150,
                expected_k = 16
            ),
            (
                name = "BB [[162, 8, 14]]",
                l = 3,
                m = 27,
                A = [(:x, 0), (:y, 10), (:y, 14)],
                B = [(:y, 12), (:x, 1), (:x, 2)],
                expected_n = 162,
                expected_k = 8
            ),
            (
                name = "BB [[180, 8, 16]]",
                l = 6,
                m = 15,
                A = [(:x, 3), (:y, 1), (:y, 2)],
                B = [(:y, 6), (:x, 4), (:x, 5)],
                expected_n = 180,
                expected_k = 8
            )
        ]
            @testset "$(case.name): l=$(case.l), m=$(case.m)" for case in test_cases
                c = GeneralizedCirculantBivariateBicycle(case.l, case.m, case.A, case.B)
                stab = parity_checks(c)
                n, k = code_n(c), code_k(c)
                nₛ, kₛ = code_n(stab), code_k(stab)
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test bivariate_bicycle_code_k(c) == case.expected_k == code_k(c)
                # Both A and B are matrices where each row and column has exactly three non-zero entries when using 3-term polynomial representation.
                Hx = parity_matrix_x(c)
                n = size(Hx,2)÷2
                A = Hx[:,1:n]
                B = Hx[:,n+1:end]
                @test all(sum(A, dims=2) .== 3)
                @test all(sum(B, dims=2) .== 3)
                @test all(sum(A, dims=1) .== 3)
                @test all(sum(B, dims=1) .== 3)
                @test A*B == B*A
            end
        end

        @testset "4-term polynomial BB [[128, 14, 12]] code from [eberhardt2024logical](@cite)" begin
            l, m = 8, 8
            A = [(:x,2), (:y,1), (:y,3), (:y,4)]
            B = [(:y,2), (:x,1), (:x,3), (:x,4)]
            c = GeneralizedCirculantBivariateBicycle(l, m, A, B)
            stab = parity_checks(c)
            n, k = code_n(c), code_k(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            H = stab_to_gf2(stab)
            mat = matrix(GF(2), H)
            computed_rank = rank(mat)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test bivariate_bicycle_code_k(c) == code_k(c)
            Hx = parity_matrix_x(c)
            n = size(Hx,2)÷2
            A = Hx[:,1:n]
            B = Hx[:,n+1:end]
            @test all(sum(A, dims=2) .== 4)
            @test all(sum(B, dims=2) .== 4)
            @test all(sum(A, dims=1) .== 4)
            @test all(sum(B, dims=1) .== 4)
            @test A*B == B*A
        end
    end
end
