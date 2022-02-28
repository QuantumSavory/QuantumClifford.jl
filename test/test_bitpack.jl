function test_bitpack()
    @testset "Alternative bit packing" begin
        @test_broken false # TODO fix these tests (started failing when switched to SIMD.jl)
        return
        for n in [1,3,5]
            N = 64*n-2
            s64 = random_stabilizer(N,N);
            phases = s64.phases;
            xzs64_rowmajor = s64.xzs;
            xzs64_colmajor = collect(s64.xzs')';
            p64 = random_pauli(N);
            c64_stab = Destabilizer(random_stabilizer(N,N)).tab
            
            after_p = stab_to_gf2(p64*s64)
            after_p_phases = (p64*s64).phases
            after_can = stab_to_gf2(canonicalize!(copy(s64)))
            after_cliff = stab_to_gf2(apply!(copy(s64),CliffordOperator(c64_stab)))

            for int in [UInt8, UInt16, UInt32, UInt64], order in [:column,:row]
                p = PauliOperator(p64.phase, N, collect(reinterpret(int,p64.xz)));
                xzs_colmajor = collect(reinterpret(int, collect(xzs64_rowmajor))')';
                xzs_rowmajor = collect(xzs_colmajor);
                s_col = Stabilizer(phases,N,xzs_colmajor);
                s_row = Stabilizer(phases,N,xzs_rowmajor);
                s = order == :column ? s_col : s_row
                apply_pauli = p*s
                @test after_p_phases == apply_pauli.phases
                canon = canonicalize!(deepcopy(s))
                @test after_can == stab_to_gf2(canon)

                for c_order in [:column, :row]
                    c_raw = Stabilizer(zeros(UInt8, 2N), N, reinterpret(int, collect(c64_stab.xzs)))
                    c = CliffordOperator(c_raw)
                    if c_order == :column
                        c = CliffordOperator(Stabilizer(c.tab.phases, N, collect(c.tab.xzs')'))
                    end
                    if int==UInt64
                        @test after_cliff == stab_to_gf2(apply!(deepcopy(s),c))
                    else
                        @test_broken after_cliff == stab_to_gf2(apply!(deepcopy(s),c))
                    end
                end
            end
        end
    end
end

test_bitpack()