function test_paulis()
    @testset "Pauli Operators" begin
        @testset "Parsing, constructors, and properties" begin
            @test P"-iXYZ" == PauliOperator(0x3, 3, vcat(BitArray([1,1,0]).chunks, BitArray([0,1,1]).chunks))
            @test P"-iXYZ" == PauliOperator(0x3, Bool[1,1,0], Bool[0,1,1])
            @test xbit(P"-iXYZ") == Bool[1,1,0]
            @test zbit(P"-iXYZ") == Bool[0,1,1]
            @test P"-iXYZ".xz == UInt64[0x03, 0x06]
            @test P"-iXYZ".phase[] == 0x03
            @test P"-iXYZ".nqubits == 3
            @test size(P"-iXYZ") == (3,)
        end
        @testset "Indexing" begin
            @test eachindex(P"IXYZ") == 1:4
            @test P"IXYZ"[3] == (true, true)
            p = P"IXYZ"
            @test p[[3,2]] == P"YX"
            p[4] = (true,false)
            @test p == P"IXYX"
        end
        @testset "Elementary operations" begin
            @test P"X"*P"Z" == P"-iY"
            @test comm(P"XX",P"YY") == 0x0
            @test comm(P"XZ",P"YZ") == 0x1
            @test prodphase(P"XX",P"YY") == 0x2
            @test prodphase(P"ZZZ",P"XXX") == 0x3
        end
        @testset "Commutation implies real phase" begin
            for i in 1:10
                for n in test_sizes
                    p1,p2 = random_pauli(n; nophase=true), random_pauli(n; nophase=true)
                    com = comm(p1,p2)==0x0
                    p = prodphase(p1,p2)
                    rea = p==0x0 || p==0x2
                    @test (com && rea) || (!com && !rea)
                end
            end
        end
        @testset "Prodphase" begin
            for i in 1:10
                for n in test_sizes
                    p1,p2 = random_pauli(n; nophase=true), random_pauli(n; nophase=true)
                    p = prodphase(p1.xz,p2.xz)
                    @test p == QuantumClifford._stim_prodphase(p1.xz,p2.xz)&0x3
                end
            end
        end
        @testset "Single qubit Paulis and their action" begin
            for i in 1:3
                for n in test_sizes
                    for t in [Stabilizer, Destabilizer, MixedDestabilizer]
                        ix, iy, iz = rand(1:n), rand(1:n), rand(1:n)
                        px = single_x(n,ix)
                        py = single_y(n,iy)
                        pz = single_z(n,iz)
                        s1 = t(random_stabilizer(n))
                        s2 = copy(s1)
                        apply!(s1,px)
                        apply_single_x!(s2,ix)
                        @test s1==s2
                        apply!(s1,py)
                        apply_single_y!(s2,iy)
                        @test s1==s2
                        apply!(s1,pz)
                        apply_single_z!(s2,iz)
                        @test s1==s2
                    end
                end
            end
        end
    end
end

test_paulis()