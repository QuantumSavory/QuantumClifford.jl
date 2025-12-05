@testitem "Alternative bit packing" tags=[:bitpack] begin
    using Random
    using QuantumClifford: Tableau, mul_left!

    @testset "alternative bit packing" begin
        for n in [1,3] # can not go higher than 4 (limitation from SIMD acting on transposed/strided arrays)
            N = 64*n-2
            s64 = random_stabilizer(N,N);
            _phases = phases(s64);
            xzs64 = tab(s64).xzs;
            xzs64T = collect(xzs64')';
            p64 = random_pauli(N;nophase=true);
            c64_stab = tab(random_destabilizer(N;phases=false));

            _after_p = p64*s64
            after_p = stab_to_gf2(_after_p);
            after_p_phases = phases(_after_p);
            after_can = stab_to_gf2(canonicalize!(copy(s64)));
            _after_clif = apply!(copy(s64),CliffordOperator(c64_stab));
            after_cliff = stab_to_gf2(_after_clif);
            after_cliff_phases = phases(_after_clif);

            after_mul = stab_to_gf2(mul_left!(copy(s64), p64))

            for int in [UInt8, UInt16, UInt32, UInt64]
                p = PauliOperator(p64.phase, N, collect(reinterpret(int,p64.xz)));
                xzs = collect(reinterpret(int, collect(xzs64)));
                xzsT = collect(xzs')';
                _s = Stabilizer(_phases,N,xzs);
                _sT = Stabilizer(_phases,N,xzsT);

                for trans in (true, false)
                    s = trans ? _sT : _s
                    apply_pauli = p*s
                    @test after_p_phases == phases(apply_pauli)
                    canon = canonicalize!(deepcopy(s))
                    @test after_can == stab_to_gf2(canon)

                    cxzs = collect(reinterpret(int, collect(c64_stab.xzs)));
                    cxzsT = collect(cxzs')';
                    _c = CliffordOperator(Tableau(zeros(UInt8,2N),N,cxzs));
                    _cT = CliffordOperator(Tableau(zeros(UInt8,2N),N,cxzsT));

                    for ctrans in (true, false)
                        c = ctrans ? _c : _cT
                        after_clifford = apply!(deepcopy(s),c)
                        @test after_cliff == stab_to_gf2(after_clifford)
                        @test after_cliff_phases == phases(after_clifford)
                    end

                    @test after_mul == stab_to_gf2(mul_left!(copy(s), p))
                end
            end
        end
    end

    @testset "fast column and fast row mul_left correctness" begin
        reinterpret_stab(s) = Stabilizer(Tableau(copy(phases(s)), nqubits(s), collect(reinterpret(UInt8, collect(s.tab.xzs)))))
        reinterpret_p(p) = PauliOperator(p.phase, nqubits(p), collect(reinterpret(UInt8, p.xz)))
        for N in [7,8,9,33,61,62,63,64,65,66]

            s0 = random_stabilizer(2,N)
            p = random_pauli(N)
            s = copy(s0)
            sr = QuantumClifford.fastrow(copy(s))
            sc = QuantumClifford.fastcolumn(copy(s))

            mul_left!(s, p)
            mul_left!(sr, p)
            mul_left!(sc, p)

            s8 = reinterpret_stab(copy(s0))
            s8r = QuantumClifford.fastrow(copy(s8))
            s8c = QuantumClifford.fastcolumn(copy(s8))
            p8 = reinterpret_p(copy(p))
            mul_left!(s8, p8)
            mul_left!(s8r, p8)
            mul_left!(s8c, p8)

            @test stab_to_gf2(s) == stab_to_gf2(sr) == stab_to_gf2(sc) == stab_to_gf2(s8) == stab_to_gf2(s8r) == stab_to_gf2(s8c)
        end
    end

    @testset "memory layout performance comparison" begin
        # fastrow should be faster than fastcolumn for canonicalization
        s_row = fastrow(random_stabilizer(100, 128))
        s_col = fastcolumn(copy(s_row))
        
        # Both layouts should produce identical results
        result_row = canonicalize!(copy(s_row); phases=true)
        result_col = canonicalize!(copy(s_col); phases=true)
        @test stab_to_gf2(result_row) == stab_to_gf2(result_col)
        
        # Test sparse gate application
        s_row_gates = fastrow(random_stabilizer(50, 64))
        s_col_gates = fastcolumn(copy(s_row_gates))
        gate = sCNOT(1, 2)
        
        # Apply sparse gates and verify identity
        s_row_after = apply!(copy(s_row_gates), gate)
        s_col_after = apply!(copy(s_col_gates), gate)
        @test stab_to_gf2(s_row_after) == stab_to_gf2(s_col_after)
        
        # Test dense clifford operator application
        s_row_clif = fastrow(random_stabilizer(50, 64))
        s_col_clif = fastcolumn(copy(s_row_clif))
        c = CliffordOperator(random_destabilizer(64; phases=false))
        
        # Apply dense gates and verify identity
        s_row_clif_after = apply!(copy(s_row_clif), c)
        s_col_clif_after = apply!(copy(s_col_clif), c)
        @test stab_to_gf2(s_row_clif_after) == stab_to_gf2(s_col_clif_after)
    end
end
