using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
using QuantumClifford: projectremoverand!

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

function alter_expect(p::PauliOperator, s::Stabilizer)
    nqubits(p) == nqubits(s) || error("The number of qubits does not match")
    n = nqubits(p)
    s_anc = s ⊗ one(Stabilizer, n)
    p_anc = zero(PauliOperator, 2n)
    for i = 1:n
        if p[i] != (false, false)
            if p[i][1]
                apply!(s_anc, sHadamard(i))
                p[i][2] && apply!(s_anc, sInvPhase(i))
            end
            apply!(s_anc, sCNOT(i, i+n))
            p_anc[i+n] = (false, true)
        end
    end
    p_anc.phase[] = p.phase[]
    _, _, result = project!(s_anc, p_anc)
    result === nothing && return 0
    result === 0x00 && return 1
    result === 0x01 && return im
    result === 0x02 && return -1
    result === 0x03 && return -im
end

@testset "Expectation of Pauli strings on stabilizer states" begin
    st = bell()
    apply!(st, sX(2))
    @test expect(P"XX", st) == alter_expect(P"XX", st) == 1
    @test expect(P"ZZ", st) == alter_expect(P"ZZ", st) == -1
    @test expect(P"XI", st) == alter_expect(P"XI", st) == 0
    @test expect(P"ZI", st) == alter_expect(P"ZI", st) == 0
    @test expect(P"II", st) == alter_expect(P"II", st) == 1
    @test expect(P"-XX", st) == alter_expect(P"-XX", st) == -1
    @test expect(P"-iXX", st) == alter_expect(P"-iXX", st) == -im
end
