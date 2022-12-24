using Random
using QuantumClifford
using QuantumClifford: Tableau

@testset "Alternative bit packing" begin
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
            end
        end
    end
end
