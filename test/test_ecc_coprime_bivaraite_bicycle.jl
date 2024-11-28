@testitem "ECC coprime Bivaraite Bicycle" begin
    using Nemo
    using Nemo: gcd
    using Hecke
    using JuMP
    using GLPK
    using Hecke: group_algebra, GF, abelian_group, gens
    using QuantumClifford
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, parity_checks

    include("test_ecc_util.jl") # minimum_distance

    function get_hx_lx(c)
        hx = stab_to_gf2(Stabilizer(parity_checks(c))[1:end÷2,:])
        lx = stab_to_gf2(logicalxview(canonicalize!(MixedDestabilizer(parity_checks(c)))))
        return hx,lx
    end

    function code_distance(hx, lx, k)
        w_values = []
        for i in 1:k
            w = minimum_distance(hx, lx[i, :])
            push!(w_values, w)
        end
        return round(Int, sum(w_values)/k)
    end

    @testset "Reproduce Table 2 wang2024coprime" begin
        # [[30,4,6]]
        l=3; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋   + 𝜋^2
        B = 𝜋 + 𝜋^3 + 𝜋^8
        c = two_block_group_algebra_codes(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 30 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 6

        # [[42,6,6]]
        l=3; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋^2 + 𝜋^3
        B = 𝜋 + 𝜋^3 + 𝜋^11
        c = two_block_group_algebra_codes(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 42 && code_k(c) == 6
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 6

        # [[70,6,8]]
        l=5; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^5;
        B = 1 + 𝜋 + 𝜋^12;
        c = two_block_group_algebra_codes(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 70 && code_k(c) == 6
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 8

        # [[108,12,6]]
        l=2; m=27;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 𝜋^2 + 𝜋^5  + 𝜋^44
        B = 𝜋^8 + 𝜋^14 + 𝜋^47
        c = two_block_group_algebra_codes(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 108 && code_k(c) == 12

        # [[126,12,10]]
        l=7; m=9
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1   + 𝜋    + 𝜋^58
        B = 𝜋^3 + 𝜋^16 + 𝜋^44
        c = two_block_group_algebra_codes(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 126 && code_k(c) == 12
    end
end
