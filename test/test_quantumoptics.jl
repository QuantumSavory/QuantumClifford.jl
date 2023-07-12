using Test
using QuantumClifford
using QuantumOpticsBase
using QuantumClifford: @S_str, random_stabilizer
using LinearAlgebra
#using QuantumCliffordQOpticsExt: _l0, _l1, _s₊, _s₋, _i₊, _i₋
const qo = Base.get_extension(QuantumClifford, :QuantumCliffordQOpticsExt)
const _l0 = qo._l0
const _l1 = qo._l1
const _s₊ = qo._s₊
const _s₋ = qo._s₋
const _i₊ = qo._i₊
const _i₋ = qo._i₋

@testset "correct conversion to Ket" begin
    for n in 1:5
        stabs = [random_stabilizer(1) for _ in 1:n]
        stab = tensor(stabs...)
        translate = Dict(S"X"=>_s₊,S"-X"=>_s₋,S"Z"=>_l0,S"-Z"=>_l1,S"Y"=>_i₊,S"-Y"=>_i₋)
        kets = [translate[s] for s in stabs]
        ket = tensor(kets...)
        @test ket.data ≈ Ket(stab).data

        rstab = random_stabilizer(n)
        lstab = random_stabilizer(n)
        lket = Ket(rstab)
        rket = Ket(lstab)
        dotket = abs(lket'*rket)
        dotstab = abs(dot(lstab,rstab))
        @test (dotket<=1e-10 && dotstab==0) || dotket≈dotstab
    end
end

@testset "correct conversion from StabMixture to Operator" begin
    for n in 1:5
        stab = random_stabilizer(n)
        @test dm(Ket(stab)) == Operator(StabMixture(stab))
    end
end
