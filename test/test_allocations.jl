@testitem "Allocation checks" begin
    using QuantumClifford
    using QuantumClifford: mul_left!
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
        @test allocated(f3) <= 1 # TODO lower it by making apply! more efficient
        f4() = apply!(s,tCNOT,[5,20])
        f4()
        allocated(f4)
        @test allocated(f4) <= 3 # TODO lower it by making apply! more efficient
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
        @test allocated(f1) <= 11
        @test allocated(f2) <= 12
        @test allocated(f3) <= 6
        @test allocated(f4) <= 5
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
end
