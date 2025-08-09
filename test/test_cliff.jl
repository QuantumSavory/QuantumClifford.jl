@testitem "Clifford" begin
    using Random
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
    using QuantumClifford: mul_left!

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Low-level tableaux ops" begin
        for n in test_sizes
            p1 = random_pauli(n)
            p11 = copy(p1)
            p2 = random_pauli(n)
            p3 = p2*p1
            s = Stabilizer([p1,p2])
            @test QuantumClifford._mul_left_nonvec!(copy(p1).xz,p2.xz)&3 == mul_left!(copy(p1).xz,p2.xz)&3
            @test prodphase(p2,p1) == mul_left!(p1,p2).phase[]
            mul_left!(p11,s,2)
            mul_left!(s,1,2)
            @test p1 == p11 == p3 == s[1]
        end
    end
    @testset "Clifford Operators" begin
        @testset "Constructors" begin
            @test_throws DimensionMismatch CliffordOperator(T"X")
        end
        @testset "Constructor from PauliOperator" begin
            for n in test_sizes
                l = random_clifford(n)
                pauli = random_pauli(n)
                @test apply!(copy(l), pauli; phases=true) == apply!(l, CliffordOperator(pauli); phases=true)
            end
        end
        @testset "Permutations of qubits" begin
            for c in [tCNOT, tId1⊗tHadamard, tCNOT⊗tCNOT, tensor_pow(tCNOT,6), tensor_pow(tCNOT,7), tensor_pow(tCNOT,6)⊗tPhase, tensor_pow(tCNOT,7)⊗tPhase]
                for rep in 1:5
                    p = randperm(nqubits(c))
                    s = random_stabilizer(nqubits(c))
                    @test permutesystems(c,p)*s[:,p] == (c*s)[:,p]
                end
            end
            for i in 1:5
                p = randperm(125)
                c = rand([tId1, tHadamard, tPhase], 125)
                @test ⊗(c[p]...) == permutesystems(⊗(c...), p)
            end
        end
        @testset "Tensor products" begin
            for n in test_sizes
                for np in [2,3,4]
                    for pow in [1,2,10]
                        s1 = random_stabilizer(n)
                        sps = [random_stabilizer(np) for i in 1:pow]
                        ss = [s1, sps...]
                        c1 = random_clifford(n)
                        cps = repeat([random_clifford(np)],pow)
                        cs = [c1, cps...]
                        res1 = ⊗([c*s for (c,s) in zip(cs,ss)]...)
                        res2 = ⊗(cs...)*⊗(ss...)
                        res3 = (c1*s1) ⊗ (⊗(tensor_pow(cps[1],pow)) * ⊗(sps...))
                        @test res1==res2==res3
                    end
                end
            end
        end
        @testset "Clifford acting on Stabilizer" begin
            for size in test_sizes
                size < 5 && continue
                s = random_stabilizer(size)
                gates = vcat([tCNOT, tHadamard, tPhase], repeat([tId1],size-4))
                gates_perm = randperm(size-1)
                gates = gates[gates_perm]
                big_gate = reduce(⊗,gates)

                s1 = apply!(copy(s),big_gate)
                @test stab_looks_good(s1)

                igates_perm = invperm(gates_perm)
                s2 = copy(s)
                canonicalize!(s2)
                s2 = apply!(s2, tCNOT, [igates_perm[1],igates_perm[1]+1])
                canonicalize!(s2)
                s2 = apply!(s2, tHadamard, [igates_perm[2]+(igates_perm[1]<igates_perm[2])])
                canonicalize!(s2)
                s2 = apply!(s2, tPhase, [igates_perm[3]+(igates_perm[1]<igates_perm[3])])

                @test canonicalize!(s1) == canonicalize!(s2)
            end
        end
        @testset "Clifford acting on (Mixed)(De)Stabilizer" begin
            for size in test_sizes
                size < 2 && continue # TODO remove this line
                d = random_destabilizer(size)
                md = MixedDestabilizer(copy(d),size÷2)
                s = copy(stabilizerview(md))
                c = random_clifford(size)
                @test QuantumClifford._apply_nonthread!(s,c) == stabilizerview(apply!(md,c)) == stabilizerview(MixedDestabilizer(apply!(d,c),size÷2))
                newsize = min(size, 5)
                indices = randperm(size)[1:newsize]
                cn = random_clifford(newsize)
                @test QuantumClifford._apply_nonthread!(s,cn,indices) == stabilizerview(apply!(md,cn,indices)) == stabilizerview(MixedDestabilizer(apply!(d,cn,indices),size÷2))
            end
        end
        @testset "Inversions" begin
            for n in test_sizes
                c = random_clifford(n)
                ci = inv(c; phases=false)
                cip = inv(c)
                id = one(c)
                @test c*cip == cip*c == id
                @test (ci*c).tab.xzs == (c*ci).tab.xzs == id.tab.xzs
            end
        end
    end
end
