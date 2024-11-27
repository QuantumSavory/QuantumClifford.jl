@testitem "ECC MILP Code Distance Solver" begin
    using Hecke
    using JuMP
    using GLPK
    using Hecke: group_algebra, GF, abelian_group, gens, one
    using QuantumClifford
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, parity_checks

    # Computing minimum distance for quantum LDPC codes is a NP-Hard problem,
    # but we can solve mixed integer linear program (MILP) for small instance sizes.
    function minimum_distance(hx,lx)
        n = size(hx, 2) # number of qubits
        m = size(hx, 1) # number of stabilizers
        # Maximum stabilizer weight
        whx = maximum(sum(hx[i, :] for i in 1:m))  # Maximum row sum of hx
        # Weight of the logical operator
        wlx = count(!iszero, lx)  # Count non-zero elements in lx
        # Number of slack variables needed to express orthogonality constraints modulo two
        num_anc_hx = ceil(Int, log2(whx))
        num_anc_logical = ceil(Int, log2(wlx))
        num_var = n + m * num_anc_hx + num_anc_logical
        model = Model(GLPK.Optimizer) # model
        set_silent(model)
        @variable(model, x[1:num_var], Bin) # binary variables
        @objective(model, Min, sum(x[i] for i in 1:n)) # objective function
        # Orthogonality to rows of hx constraints
        for row in 1:m
            weight = zeros(Int, num_var)
            supp = findall(hx[row, :] .!= 0)  # Non-zero indices in hx[row, :]
            for q in supp
                weight[q] = 1
            end
            cnt = 1
            for q in 1:num_anc_hx
                weight[n + (row - 1) * num_anc_hx + q] = -(1 << cnt)
                cnt += 1
            end
            @constraint(model, sum(weight[i] * x[i] for i in 1:num_var) == 0)
        end
        # Odd overlap with lx constraint
        supp = findall(lx .!= 0)  # Non-zero indices in lx
        weight = zeros(Int, num_var)
        for q in supp
            weight[q] = 1
        end
        cnt = 1
        for q in 1:num_anc_logical
            weight[n + m * num_anc_hx + q] = -(1 << cnt)
            cnt += 1
        end
        @constraint(model, sum(weight[i] * x[i] for i in 1:num_var) == 1)
        optimize!(model)
        opt_val = sum(value(x[i]) for i in 1:n)
        return Int(opt_val)
    end

    function get_hx_lx(c)
        hx = stab_to_gf2(Stabilizer(parity_checks(c))[1:end÷2,:])
        lx = stab_to_gf2(logicalxview(canonicalize!(MixedDestabilizer(parity_checks(c)))))
        return hx, lx
    end

    @testset "Reproduce Table 3 bravyi2024high" begin
        # [[72, 12, 6]]
        l=6; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[90, 8, 10]]
        l=15; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y   + y^2
        B = 1   + x^2 + x^7
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 10

        # [[108, 8, 10]]
        l=9; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 10

        # [[144, 12, 12]]
        l=12; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 12
    end

    @testset "Reproduce Table 1 berthusen2024toward" begin
        # [[72, 8, 6]]
        l=12; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y + y^2
        B = 1   + x + x^11
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[90, 8, 6]]
        l=9; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^8 + y^4 + y
        B = y^5 + x^8 + x^7
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[120, 8, 8]]
        l=12; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^10 + y^4 + y
        B = 1    + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 8

        # [[150, 8, 8]]
        l=15; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^5 + y^2 + y^3
        B = y^2 + x^7 + x^6
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 8

        # [[196, 12, 8]]
        l=14; m=7
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^6 + y^5 + y^6
        B = 1   + x^4 + x^13
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 8
    end

    @testset "Reproduce Table 1 wang2024coprime" begin
        # [[54, 8, 6]]
        l=3; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^2 + y^4
        B = y^3 + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[98, 6, 12]]
        l=7; m=7
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y^5 + y^6
        B = y^2 + x^3 + x^5
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 12

        # [[126, 8, 10]]
        l=3; m=21
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^2 + y^10
        B = y^3 + x  +  x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 10

        # [[150, 16, 8]]
        l=5; m=15
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^6 + y^8
        B = y^5 + x   + x^4
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 8

        # [[162, 8, 14]]
        l=3; m=27
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1    + y^10 + y^14
        B = y^12 + x    + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 14

        # [[180, 8, 16]]
        l=6; m=15
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y   + y^2
        B = y^6 + x^4 + x^5
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 16
    end
end
