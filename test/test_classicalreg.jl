@testitem "Classical" begin
    using Random
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good

    n=5
    stab = random_stabilizer(n)
    mdstab = MixedDestabilizer(stab)
    _ = Register(stab)
    reg = Register(stab, [0,0,0,0])
    regmd = Register(mdstab, [0,0,0,0])
    @test reg==regmd
    @test stabilizerview(reg) == stabilizerview(mdstab)
    @test destabilizerview(reg) == destabilizerview(mdstab)
    @test logicalxview(reg) == logicalxview(mdstab)
    @test logicalzview(reg) == logicalzview(mdstab)
    @test bitview(reg) == bitview(regmd)
    @test quantumstate(reg) == mdstab
    for state in [mdstab,reg,regmd]
        for op in [sMX(1,1),sMY(2,2),sMZ(3,3),PauliMeasurement(P"XYZZZ",4),sCNOT(1,2),sCPHASE(2,3),sCNOT(3,4),NoiseOpAll(UnbiasedUncorrelatedNoise(0.1))]
            apply!(state,op)
        end
        for (i,proj) in enumerate([projectXrand!, projectYrand!, projectZrand!])
            proj(state, i)
        end
    end
    @test tab(canonicalize!(stabilizerview(reg))).xzs == tab(canonicalize!(stabilizerview(mdstab))).xzs
end
