@testitem "CirculantBivariateBicycle" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using QECCore
        using QECCore.LinearAlgebra
        using QuantumClifford
        using QuantumClifford: stab_looks_good
        using QuantumClifford: stab_to_gf2
        using QuantumClifford.ECC
        using Oscar

        @testset "CirculantBivariateBicycle codes from [bravyi2024high](@cite)" begin
            test_cases = [
                (name = "BB [[108, 8,  10]]", l =  9, m =  6, A = [3, 1, 2], B = [3, 1, 2], expected_n = 108, expected_k = 8),
                (name = "BB [[90,  8,  10]]", l = 15, m =  3, A = [9, 1, 2], B = [0, 2, 7], expected_n =  90, expected_k = 8),
                (name = "BB [[288, 12, 18]]", l = 12, m = 12, A = [3, 2, 7], B = [3, 1, 2], expected_n = 288, expected_k = 12),
                (name = "BB [[144, 12, 12]]", l = 12, m =  6, A = [3, 1, 2], B = [3, 1, 2], expected_n = 144, expected_k = 12),
                (name = "BB [[72,  12,  6]]", l =  6, m =  6, A = [3, 1, 2], B = [3, 1, 2], expected_n =  72, expected_k = 12),
                (name = "BB [[360, 12,≤24]]", l = 30, m =  6, A = [9, 1, 2], B = [3,25,26], expected_n = 360, expected_k = 12),
                (name = "BB [[756, 18,≤34]]", l = 21, m = 18, A = [3,10,17], B = [5, 3,19], expected_n = 756, expected_k = 16),
                (name = "BB [[784, 24, d ]]", l = 28, m = 14, A = [26,6, 8], B = [7, 9,20], expected_n = 784, expected_k = 24),
            ]
            @testset "$(case.name): l=$(case.l), m=$(case.m)" for case in test_cases
                c = CirculantBivariateBicycle(case.l, case.m, case.A, case.B)
                @test code_k(c) == case.expected_k
                stab = parity_checks(c)
                n, k = code_n(c), code_k(c)
                nₛ, kₛ = code_n(stab), code_k(stab)
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k && computed_rank == nₛ - kₛ
                @test n == nₛ == case.expected_n && k == kₛ == case.expected_k
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test bivariate_bicycle_code_k(c) == case.expected_k == QECCore.code_k(c)
            end
        end

        @testset "Verify number of logical qubits `k` from Table 1: [berthusen2024toward](@cite)" begin
            test_cases = [
                (name = "BB [[72 , 8,6]]", l = 12, m = 3, A = [9, 1, 2], B = [0, 1,11], expected_n = 72, expected_k = 8),
                (name = "BB [[90 , 8,6]]", l =  9, m = 5, A = [8, 4, 1], B = [5, 8, 7], expected_n = 90, expected_k = 8),
                (name = "BB [[120, 8,8]]", l = 12, m = 5, A = [10,4, 1], B = [0, 1, 2], expected_n = 120, expected_k = 8),
                (name = "BB [[150, 8,8]]", l = 15, m = 5, A = [5, 2, 3], B = [2, 7, 6], expected_n = 150, expected_k = 8),
                (name = "BB [[196,12,8]]", l = 14, m = 7, A = [6, 5, 6], B = [0, 4,13], expected_n = 196, expected_k = 12)
            ]
            @testset "$(case.name): a=$(case.l), b=$(case.m)" for case in test_cases
                c = CirculantBivariateBicycle(case.l, case.m, case.A, case.B)
                @test code_k(c) == case.expected_k
                stab = parity_checks(c)
                n, k = code_n(c), code_k(c)
                nₛ, kₛ = code_n(stab), code_k(stab)
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k && computed_rank == nₛ - kₛ
                @test n == nₛ == case.expected_n && k == kₛ == case.expected_k
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test bivariate_bicycle_code_k(c) == case.expected_k == QECCore.code_k(c)
            end
        end

        @testset "Verify number of logical qubits `k` from Table 1: wang2024coprime" begin
            # Test cases from [wang2024coprime](@cite)
            test_cases = [
                (name = "BB [[54 ,  8,  6]]", a = 3,  b = 9,  ax = [0, 2, 4],  bx = [3, 1, 2],  expected_k = 8),
                (name = "BB [[98 ,  6, 12]]", a = 7,  b = 7,  ax = [3, 5, 6],  bx = [2, 3, 5],  expected_k = 6),
                (name = "BB [[126,  8, 10]]", a = 3,  b = 21, ax = [0, 2,10],  bx = [3, 1, 2],  expected_k = 8),
                (name = "BB [[150, 16,  8]]", a = 5,  b = 15, ax = [0, 6, 8],  bx = [5, 1, 4],  expected_k = 16),
                (name = "BB [[162,  8, 14]]", a = 3,  b = 27, ax = [0,10,14],  bx = [12,1, 2],  expected_k = 8),
                (name = "BB [[180,  8, 16]]", a = 6,  b = 15, ax = [3, 1, 2],  bx  =[6, 4, 5],  expected_k = 8)
            ]
            @testset "$(case.name): a=$(case.a), b=$(case.b)" for case in test_cases
                c = CirculantBivariateBicycle(case.a, case.b, case.ax, case.bx)
                @test code_k(c) == case.expected_k
                stab = parity_checks(c)
                n, k = code_n(c), code_k(c)
                nₛ, kₛ = code_n(stab), code_k(stab)
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test bivariate_bicycle_code_k(c) == case.expected_k == QECCore.code_k(c)
            end
        end
    end
end
