using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
using Test


@testset "Stabilizer to Clifford Conversion" begin
    for i in 1:10
        stab = random_stabilizer(4)
        cl = CliffordOperator(stab)
        cltab = tab(cl)
        @test Stabilizer(cltab[end÷2+1:end]) == (stab);
    end


    for i in 1:10

        gnd4 = S"ZIII
                 IZII
                 IIZI
                 IIIZ"

        stab = random_stabilizer(4)
        cl = CliffordOperator(stab)
        cltab = tab(cl)
        @test apply!(gnd4, cl) == (stab);
    end
    ## TODO: More tests??
end

@testset "Basic Stabs to Clifford" begin
    stabSWAP = S"IZ
                 ZI"
    @test CliffordOperator(stabSWAP) == tSWAP
    ## TODO: More tests??
end
