@testitem "Random" begin
    using QuantumClifford
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Random sampling of operators" begin
        for n in [1, test_sizes..., 200,500]
            p = random_pauli(n)
            s = random_stabilizer(n)
            ss = random_stabilizer(rand(1:n),n)
            ms = MixedDestabilizer(ss)
            d = random_destabilizer(n)
            c = random_clifford(n)
            sq = random_clifford1(n÷2+1)
            @test stab_looks_good(s)
            @test stab_looks_good(ss)
            @test destab_looks_good(d)
            @test mixed_destab_looks_good(ms)
            @test stab_looks_good(c*s)
            @test stab_looks_good(c*ss)
            @test destab_looks_good(c*d)
            @test mixed_destab_looks_good(c*ms)
            @test stab_looks_good(p*s)
            @test stab_looks_good(p*ss)
            @test destab_looks_good(p*d)
            @test mixed_destab_looks_good(p*ms)
            @test stab_looks_good(apply!(s,sq,phases=false))
            @test stab_looks_good(apply!(ss,sq,phases=false))
            @test destab_looks_good(apply!(d,sq,phases=false))
            @test mixed_destab_looks_good(apply!(ms,sq,phases=false))
        end
    end

    @testset "Random Paulis" begin
        for n in [1, test_sizes..., 200,500]
            @test all((random_pauli(n).phase[] == 0  for _ in 1:100))
            @test all((random_pauli(n, 0.1).phase[] == 0  for _ in 1:100))
            @test any((random_pauli(n; nophase=false, realphase=false).phase[] == 1  for _ in 1:100))
            @test any((random_pauli(n, 0.1; nophase=false, realphase=false).phase[] == 1  for _ in 1:100))
            @test any((random_pauli(n; nophase=false).phase[] ∈ [0,2]  for _ in 1:100))
            @test any((random_pauli(n, 0.1; nophase=false).phase[] ∈ [0,2]  for _ in 1:100))
        end
        for i in 1:10
            e = 0.2
            n = 10000
            expected = 2/3*e * 2 * n
            bound = 1/sqrt(n)
            @test expected * (1-10bound) <= sum(count_ones.(random_pauli(10000,0.2).xz)) <= expected * (1+10bound)
            e = 0.75
            n = 10000
            expected = 2/3*e * 2 * n
            bound = 1/sqrt(n)
            @test expected * (1-10bound) <= sum(count_ones.(random_pauli(10000).xz)) <= expected * (1+10bound)

        end
    end
end
