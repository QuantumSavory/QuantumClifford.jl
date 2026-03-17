@testitem "Pauli Operators" begin
    using QuantumClifford: apply_single_x!, apply_single_y!, apply_single_z!
    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Parsing, constructors, and properties" begin
        @test P"-iXYZ" == PauliOperator(0x3, 3, vcat(BitArray([1,1,0]).chunks, BitArray([0,1,1]).chunks))
        @test P"-iXYZ" == PauliOperator(0x3, Bool[1,1,0], Bool[0,1,1])
        @test P"XYZ" == PauliOperator(Bool[1,1,0],Bool[0,1,1])
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
        @test firstindex(p):lastindex(p) == eachindex(p) == 1:length(p)
        @test p[[3,2]] == P"YX"
        p[4] = (true,false)
        @test p == P"IXYX"
    end
    @testset "Elementary operations" begin
        @test P"X"*P"Z" == P"-iY"
        @test P"X"⊗P"Z" == P"XZ"
        @test -P"X" == P"-X"
        @test +P"X" == P"X"
        @test 1*P"X" == P"X"
        @test 1.0*P"X" == P"X"
        @test -1*P"X" == P"-X"
        @test im*P"X" == P"iX"
        @test -im*P"X" == P"-iX"
        @test_throws DomainError 0.5*P"X"
        @test comm(P"XX",P"YY") == 0x0
        @test comm(P"XZ",P"YZ") == comm(S"II XZ",P"YZ", 2) == comm(P"XZ",S"II YZ", 2) == comm(S"II XZ",P"YZ")[2] == comm(P"XZ",S"II YZ")[2] == 0x1
        @test prodphase(P"XX",P"YY") == 0x2
        @test prodphase(P"ZZZ",P"XXX") == prodphase(S"III ZZZ",P"XXX",2) == prodphase(P"ZZZ",S"III XXX",2) == prodphase(S"III ZZZ",S"III XXX",2,2) == 0x3
    end

    for Pop in [P"X", P"iX", P"-iXYZ", random_pauli(100; nophase=false, realphase=false)]
        @test Pop * inv(Pop) == zero(Pop)
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
                    x, y, z = rand(1:n), rand(1:n), rand(1:n)
                    px = single_x(n,x)
                    py = single_y(n,y)
                    pz = single_z(n,z)
                    rstab = random_stabilizer(n)
                    s1 = t(rstab)
                    s2 = copy(s1)
                    s3 = copy(s1)
                    apply!(s1,px)
                    apply_single_x!(s2,x)
                    apply!(s3,P"X",[x])
                    @test s1==s2==s3
                    apply!(s1,py)
                    apply_single_y!(s2,y)
                    apply!(s3,P"Y",[y])
                    @test s1==s2==s3
                    apply!(s1,pz)
                    apply_single_z!(s2,z)
                    apply!(s3,P"Z",[z])
                    @test s1==s2==s3
                end
            end
        end
    end
end
