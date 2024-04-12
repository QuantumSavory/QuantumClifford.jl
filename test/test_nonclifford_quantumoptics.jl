using QuantumClifford
using QuantumClifford: GeneralizedStabilizer, rowdecompose, PauliChannel, mul_left!, mul_right!
using QuantumClifford: @S_str, random_stabilizer
using QuantumOpticsBase
using LinearAlgebra
using Test
using InteractiveUtils
using Random

#using QuantumCliffordQOpticsExt: _l0, _l1, _s₊, _s₋, _i₊, _i₋
const qo = Base.get_extension(QuantumClifford, :QuantumCliffordQOpticsExt)
const _l0 = qo._l0
const _l1 = qo._l1
const _s₊ = qo._s₊
const _s₋ = qo._s₋
const _i₊ = qo._i₊
const _i₋ = qo._i₋

qo_basis = SpinBasis(1//2)
qo_tgate = sparse(identityoperator(qo_basis))
qo_tgate.data[2,2] = exp(im*pi/4)

##

for n in 1:5
    i = rand(1:n)
    stab = random_stabilizer(n)
    genstab = GeneralizedStabilizer(stab)
    ket = Ket(stab)
    @test dm(ket) ≈ Operator(stab)
    @test dm(ket) ≈ Operator(genstab)

    pauli = random_pauli(n; nophase=false, realphase=true)
    qo_pauli = Operator(pauli)

    qo_bigtgate = n==1 ? qo_tgate : embed(qo_basis^n, i, qo_tgate)
    bigtgate = embed(n,i, pcT)
    @test qo_bigtgate ≈ Operator(bigtgate)

    for step in 1:10
        # apply!(ket, qo_bigtgate) TODO implement this API
        ket = qo_bigtgate*ket
        apply!(genstab, bigtgate)
        @test dm(ket) ≈ Operator(genstab)
        @test isapprox(expect(qo_pauli, ket), expect(pauli, genstab); atol=1e-5)
    end
end
