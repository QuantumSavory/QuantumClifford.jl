@testitem "Symbolic Measurements" begin
    using Random
    using QuantumClifford
    using Test
    test_sizes = [10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Symbolic Measurements" begin
        for n in test_sizes
            for i in [rand(1:n), 1, n, n÷2+1, n÷2-1, n÷2]
                for T in [sMRZ,sMZ,sMRX,sMX] # TODO sMY sMRY
                    for _ in 1:10
                        apply!(PauliFrame(n,n,n), T(i,1))
                    end
                end
            end
            for i in [rand(1:n), 1, n, n÷2+1, n÷2-1, n÷2]
                for T in [sMRZ,sMZ,sMRX,sMX,sMY] # TODO sMRY
                    for _ in 1:10
                        s = random_stabilizer(n)
                        apply!(Register(copy(s),1), T(i,1))
                    end
                end
            end
        end
        @test_throws DimensionMismatch SingleQubitOperator(tCNOT,1)
        @test_throws DimensionMismatch CliffordOperator(sHadamard(5),2)
        @test_throws ArgumentError CliffordOperator(sHadamard(5),6,compact=true)
        for T in [sMRZ,sMZ,sMRX,sMX,sMRY,sMY]
            @test_throws ArgumentError T(-1)
            @test_throws ArgumentError T(-1,0)
            @test_throws ArgumentError T(-1,1)
            @test T(1,0) == T(1)
        end
    end
end
