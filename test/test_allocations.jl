@testitem "Allocation checks" begin
    using QuantumClifford: mul_left!, RandDestabMemory, Tableau
    using Random, InteractiveUtils
    n = Threads.nthreads()
    allocated(f::F) where {F} = @allocations f()
    @testset "apply! mul_left! canonicalize!" begin
        p1 = random_pauli(500)
        p2 = random_pauli(500)
        f1() = mul_left!(p1,p2)
        f1()
        allocated(f1)
        @test allocated(f1) == 0
        s = random_stabilizer(500)
        f2() = canonicalize!(s)
        f2()
        allocated(f2)
        @test allocated(f2) <= 1
        f2a() = QuantumClifford._canonicalize!(s)
        f2a()
        allocated(f2a)
        @test allocated(f2a) <= 1
        c = random_clifford(500)
        f3() = apply!(s,c)
        f3()
        allocated(f3)
        #@test allocated(f3) <= 1 # TODO lower it by making apply! more efficient
        @test_broken false # the test above does not always work on julia 1.11+, depending on whether it runs in CI or not
        f4() = apply!(s,tCNOT,[5,20])
        f4()
        allocated(f4)
        #@test allocated(f4) <= 3 # TODO lower it by making apply! more efficient
        @test_broken false # the test above does not always work on julia 1.11+, depending on whether it runs in CI or not
        for phases in [(false,false),(false,true),(true,false),(true,true)], i in 1:6
            g = enumerate_single_qubit_gates(i,qubit=10,phases=phases)
            f5() = apply!(s,g)
            f5()
            allocated(f5)
            @test allocated(f5) <= 2
        end
        for g in [sSWAP(10,200), sCNOT(10,200)]
            f6() = apply!(s,g)
            f6()
            allocated(f6)
            @test allocated(f6) <= 2
        end
    end
    @testset "random_destabilizer" begin
        N = 100
        memory = RandDestabMemory(N)
        f1() = RandDestabMemory(N)
        f1()
        f2() = random_destabilizer(Random.GLOBAL_RNG, memory)
        f2()
        f3() = Destabilizer(Tableau(memory.phasesarray, memory.xzs))
        f3()
        @test allocated(f1) < 12.5 * N^2 + 50 * N + 1000
        if VERSION >= v"1.11"
            @test abs(allocated(f2) - allocated(f3)) / allocated(f3) < 0.05
        end
    end
    @testset "project!" begin
        N = 100
        d = random_destabilizer(N)
        md = MixedDestabilizer(random_destabilizer(N))
        md.rank = 50
        s = copy(stabilizerview(d))
        ms = MixedStabilizer(s)
        ms.rank = 50
        p = s[end];
        f1() = project!(s,p)
        f1()
        f2() = project!(ms,p)
        f2()
        f3() = project!(d,p)
        f3()
        f4() = project!(md,p)
        f4()
        allocated(f1)
        allocated(f2)
        allocated(f3)
        allocated(f4)
        @test allocated(f1) <= 15
        @test allocated(f2) <= 12
        @test allocated(f3) <= 6
        @test allocated(f4) <= 8
        for p! in [projectX!, projectY!, projectZ!]
            md = MixedDestabilizer(random_destabilizer(N))
            md.rank = 50
            f5() = p!(md,40)
            f5()
            allocated(f5)
            @test allocated(f5) <= 7
        end
    end
    @testset "tensor product" begin
        stabs = [s[1:5] for s in [random_stabilizer(n) for n in [63,64,65,127,128,129]]]
        f1() = ⊗(stabs...)
        f1()
        allocated(f1)
        @test allocated(f1) <= 18
    end

    test_sizes = [2,63,64,65,127,128,129,511,512,513]
    @testset "apply_right! symbolic" begin
        for q in test_sizes
            q1 = rand(1:q)
            q2 = rand(setdiff(1:q, [q1]))
            for _gate in subtypes(AbstractSingleQubitOperator)
                _gate == SingleQubitOperator && continue

                l = random_clifford(q)
                gate = _gate(q1)
                f1() = apply_right!(l, gate)
                f1()
                allocated(f1)
                @test allocated(f1) == 0
            end
            for _gate in subtypes(AbstractTwoQubitOperator)
                l = random_clifford(q)
                gate = _gate(q1, q2)
                f2() = apply_right!(l, gate)
                f2()
                allocated(f2)
                @test allocated(f2) == 0
            end
        end
    end
    @testset "apply_right! pauli" begin
        for q in test_sizes
            l = random_clifford(q)
            pauli = random_pauli(q)
            f1() = apply_right!(l, pauli)
            f1()
            allocated(f1)
            @test allocated(f1) == 0
        end
    end
    @testset "apply_right! dense" begin
        for q in test_sizes
            l = random_clifford(q)
            r = random_clifford(q)
            f1() = apply_right!(l, r)
            f1()
            allocated(f1)
            @test_broken allocated(f1) == 0
        end
    end
end
