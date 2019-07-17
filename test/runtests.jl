using SimpleClifford, Test, Random

function stab_looks_good(s)
    c = canonicalize!(s)
    phasesok = all((c.phases .== 0x0) .| (c.phases .== 0x2))
    H = stab_to_gf2(c)
    good_indices = reduce(|,H,dims=(1,))
    good_indices = good_indices[1:end÷2] .| good_indices[end÷2+1:end]
    rowsok = all(good_indices)
    good_indices = reduce(|,H,dims=(2,))
    colsok = all(good_indices)
    return phasesok && rowsok && colsok && check_allrowscommute(c)
end

function destab_looks_good(destabilizer)
    s = destabilizer.stabilizer
    d = destabilizer.destabilizer
    good = true
    for i in eachindex(s)
        good &= comm(s[i],d[i])==0x1
        for j in 1:lastindex(s)
            j==i && continue
            good &= comm(s[i],d[j])==0x0
        end
    end
    good
end

function tests()

Random.seed!(42)

@testset "Pauli Operators" begin
    @testset "Parsing, constructors, and properties" begin
        @test P"-iXYZ" == PauliOperator(0x3, 3, vcat(BitArray([1,1,0]).chunks, BitArray([0,1,1]).chunks))
        @test P"-iXYZ" == PauliOperator(0x3, Bool[1,1,0], Bool[0,1,1])
        @test P"-iXYZ".xbit == Bool[1,1,0]
        @test P"-iXYZ".xz == UInt64[0x03, 0x06]
        @test P"-iXYZ".phase[] == 0x03 # TODO why is this failing?
        @test P"-iXYZ".nqbits == 3
        @test size(P"-iXYZ") == (3,)
    end
    
    @testset "Elementary operations" begin
        @test P"X"*P"Z" == P"-iY"
        @test comm(P"XX",P"YY") == 0x0
        @test comm(P"XZ",P"YZ") == 0x1
        @test prodphase(P"XX",P"YY") == 0x2
        @test prodphase(P"ZZZ",P"XXX") == 0x3
    end
    
    @testset "Commutation implies real phase" begin
        for i in 10
            for n in [63,64,65,127,128,129]
                p1,p2 = random_pauli(n; nophase=true), random_pauli(n; nophase=true)
                com = comm(p1,p2)==0x0
                p = prodphase(p1,p2)
                rea = p==0x0 || p==0x2
                @test (com && rea) || (!com && !rea)
            end
        end
    end
end

@testset "Stabilizer canonicalization" begin
    s = S"- XZZZZ_____
          - _YZY___YX_
          - __XXZ__YX_
          + Z_Z_Y__YXZ
          + _____Z____
          + __________
          + __________
          + ______YYX_
          + __________
          + __________"
    canonicalize!(s)
    t = S"- XZZZZ_____
          - _YZY___YX_
          - __XXZ__YX_
          + Z_Z_Y__YXZ
          + ______YYX_
          + _____Z____
          + __________
          + __________
          + __________
          + __________"
    @test s == t
    for n in [10,63,64,65,127,128,129]
        @test stab_looks_good(random_stabilizer(n))
    end
end

@testset "Projective measurements" begin
    s = S"XXX
          ZZI
          IZZ"
    ps, anticom, res = project!(copy(s), P"ZII")
    ps = canonicalize!(ps)
    @test anticom==1 && isnothing(res) && ps == S"ZII
                                                  IZI
                                                  IIZ"
    @test stab_looks_good(ps)

    ps, anticom, res = project!(copy(s), P"-XXX")
    @test anticom==0 && res[]==0x2 && ps == canonicalize!(copy(s))
    @test stab_looks_good(ps)

    ps, anticom, res = project!(copy(s), P"-XXX"; keep_result=false)
    @test anticom==0 && isnothing(res) && ps == s
    @test stab_looks_good(ps)

    for n in [10,63,64,65,127,128,129]
        s = random_stabilizer(n)
        m = random_pauli(n;nophase=true)
        ps, anticom, res = project!(copy(s),m)
        @test anticom==0x0 || ps[anticom]==m
        @test stab_looks_good(ps)
        m = single_z(n,1)
        ps, anticom, res = project!(copy(s),m)
        @test anticom==0x0 || ps[anticom]==m
        @test stab_looks_good(ps)
        m = single_x(n,1)
        ps, anticom, res = project!(copy(s),m)
        @test anticom==0x0 || ps[anticom]==m
        @test stab_looks_good(ps)
    end
end

@testset "Destabilizer initialization" begin
    for n in [10,63,64,65,127,128,129]
        @test destab_looks_good(calculate_destabilizer(random_stabilizer(n)))
    end
end

@testset "Projective measurements with destabilizers" begin
    for n in [10,63,64,65,127,128,129]
        s = canonicalize!(random_stabilizer(n))
        m = random_pauli(n;nophase=true)
        ps, anticom, res = project!(copy(s),m)
        dps, danticom, dres = project!(calculate_destabilizer(copy(s)),m)
        @test destab_looks_good(dps)
        @test anticom==danticom && res==dres && canonicalize!(ps)==canonicalize!(dps.stabilizer)
        m = single_z(n,1)
        ps, anticom, res = project!(copy(s),m)
        dps, danticom, dres = project!(calculate_destabilizer(copy(s)),m)
        @test destab_looks_good(dps)
        @test anticom==danticom && res==dres && canonicalize!(ps)==canonicalize!(dps.stabilizer)
        m = single_x(n,1)
        ps, anticom, res = project!(copy(s),m)
        dps, danticom, dres = project!(calculate_destabilizer(copy(s)),m)
        @test destab_looks_good(dps)
        @test anticom==danticom && res==dres && canonicalize!(ps)==canonicalize!(dps.stabilizer)
    end
end

end

tests()
