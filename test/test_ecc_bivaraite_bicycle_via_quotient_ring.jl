@testitem "ECC Bivaraite Bicycle" tags=[:ecc] begin
    using Hecke
    using HiGHS
    using JuMP
    using Hecke: group_algebra, GF, abelian_group, gens, one
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, distance, DistanceMIPAlgorithm

    @testset "Reproduce Table 3 bravyi2024high" begin
        # [[72, 12, 6]]
        l=6; m=6
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y + y^2)
        B = S(y^3 + x + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 72 && code_k(c) == 12
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[90, 8, 10]]
        l=15; m=3
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^9 + y   + y^2)
        B = S(1   + x^2 + x^7)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 90 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[108, 8, 10]]
        l=9; m=6
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y + y^2)
        B = S(y^3 + x + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 108 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[144, 12, 12]]
        l=12; m=6
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y + y^2)
        B = S(y^3 + x + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 144 && code_k(c) == 12
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[288, 12, 12]]
        l=12; m=12
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y^2 + y^7)
        B = S(y^3 + x   + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 288 && code_k(c) == 12

        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[360, 12, ≤ 24]]
        l=30; m=6
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^9 + y    + y^2)
        B = S(y^3 + x^25 + x^26)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 360 && code_k(c) == 12
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[756, 16, ≤ 34]]
        l=21; m=18
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y^10 + y^17)
        B = S(y^5 + x^3  + x^19)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 756 && code_k(c) == 16
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)
    end

    @testset "Reproduce Table 1 berthusen2024toward" begin
        # [[72, 8, 6]]
        l=12; m=3
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^9 + y + y^2)
        B = S(1   + x + x^11)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 72 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[90, 8, 6]]
        l=9; m=5
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^8 + y^4 + y)
        B = S(y^5 + x^8 + x^7)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 90 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[120, 8, 8]]
        l=12; m=5
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^10 + y^4 + y)
        B = S(1    + x   + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 120 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[150, 8, 8]]
        l=15; m=5
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^5 + y^2 + y^3)
        B = S(y^2 + x^7 + x^6)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 150 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[196, 12, 8]]
        l=14; m=7
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^6 + y^5 + y^6)
        B = S(1   + x^4 + x^13)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 196 && code_k(c) == 12
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)
    end

    @testset "Reproduce Table 1 wang2024coprime" begin
        # [[54, 8, 6]]
        l=3; m=9
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1   + y^2 + y^4)
        B = S(y^3 + x   + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 54 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[98, 6, 12]]
        l=7; m=7
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y^5 + y^6)
        B = S(y^2 + x^3 + x^5)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 98 && code_k(c) == 6
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[126, 8, 10]]
        l=3; m=21
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1   + y^2 + y^10)
        B = S(y^3 + x  +  x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 126 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[150, 16, 8]]
        l=5; m=15
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1   + y^6 + y^8)
        B = S(y^5 + x   + x^4)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 150 && code_k(c) == 16
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[162, 8, 14]]
        l=3; m=27
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1    + y^10 + y^14)
        B = S(y^12 + x    + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 162 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[180, 8, 16]]
        l=6; m=15
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y   + y^2)
        B = S(y^6 + x^4 + x^5)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 180 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)
    end

    @testset "Reproduce Table 1 eberhardt2024logical" begin
        # [[108, 16, 6]]
        l=6; m=9
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1 +   y   + y^2)
        B = S(y^3 + x^2 + x^4)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 108 && code_k(c) == 16
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[128, 14, 12]]
        l=8; m=8
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^2 + y + y^3 + y^4)
        B = S(y^2 + x + x^3 + x^4)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 128 && code_k(c) == 14
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 4)
        @test all(sum(B, dims=2) .== 4)
        @test all(sum(A, dims=1) .== 4)
        @test all(sum(B, dims=1) .== 4)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[162, 4, 16]]
        l=9; m=9
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1   + x + y)
        B = S(x^3 + y + y^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 162 && code_k(c) == 4
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[162, 12, 8]]
        l=9; m=9
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1   + x   + y^6)
        B = S(y^3 + x^2 + x^3)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 162 && code_k(c) == 12
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[162, 24, 6]]
        l=9; m=9
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(1   + y   + y^2)
        B = S(y^3 + x^3 + x^6)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 162 && code_k(c) == 24
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[270, 8, 18]]
        l=9; m=15
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y + y^2)
        B = S(y^3 + x + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 270 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[98, 6, 12]]
        l=7; m=7
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x + y^3 + y^4)
        B = S(y + x^3 + x^4)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 98 && code_k(c) == 6
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)

        # [[162, 8, 12]]
        l=9; m=9
        R, (x, y) = polynomial_ring(GF(2), [:x, :y])
        I = ideal(R, [x^l-1, y^m-1])
        S, _ = quo(R, I)
        A = S(x^3 + y + y^2)
        B = S(y^3 + x + x^2)
        c = BivariateBicycleCode(l, m, A, B)
        @test code_n(c) == 162 && code_k(c) == 8
        Hx = parity_matrix_x(c)
        n = size(Hx,2)÷2
        A = Hx[:,1:n]
        B = Hx[:,n+1:end]
        @test all(sum(A, dims=2) .== 3)
        @test all(sum(B, dims=2) .== 3)
        @test all(sum(A, dims=1) .== 3)
        @test all(sum(B, dims=1) .== 3)
        @test A*B == B*A
        @test iszero(A*B+B*A) == iszero(2*A*B)
    end
end
