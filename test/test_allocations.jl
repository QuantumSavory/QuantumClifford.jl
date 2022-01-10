function test_allocations()
    VERSION >= v"1.7" || return
    @testset "Allocation checks" begin
        n = Threads.nthreads()
        allocated(f::F) where {F} = @allocated f()
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
            @test allocated(f2) < 50
            c = random_clifford(500)
            f3() = apply!(s,c)
            f3()
            @test allocated(f3) < 1500*n # TODO lower it by making apply! more efficient
            f4() = apply!(s,CNOT,[5,20])
            f4()
            @test allocated(f4) < 1500*n # TODO lower it by making apply! more efficient
            for phases in [(0,0),(0,1),(1,0),(1,1)], i in 1:6
                g = enumerate_single_qubit_gates(i,qubit=10,phases=phases)
                f5() = apply!(s,g)
                f5()
                @test allocated(f5) < 100*n
            end
            for g in [sSWAP(10,200), sCNOT(10,200)]
                f6() = apply!(s,g)
                f6()
                @test allocated(f6) < 150*n
            end 
        end
        # TODO project! and the other canonicalize           
        @testset "tensor product" begin
            stabs = [s[1:5] for s in [random_stabilizer(n) for n in [63,64,65,127,128,129]]]
            f1() = ⊗(stabs...)
            f1()
            @test allocated(f1) < 6000
        end
    end
end

test_allocations()