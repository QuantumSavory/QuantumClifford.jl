using Test
using QuantumClifford

function alter_expect(s::Stabilizer, p::PauliOperator)
    nqubits(p) == nqubits(s) || error("The number of qubits does not match")
    n = nqubits(p)
    s_anc = s ⊗ one(Stabilizer, n)
    p_anc = ['I' for _ = 1:2n]
    for i = 1:n
        if p[i] != (false, false)
            if p[i][1]
                apply!(s_anc, sHadamard(i))
                p[i][2] && apply!(s_anc, sInvPhase(i))
            end
            apply!(s_anc, sCNOT(i, i+n))
            p_anc[i+n] = 'Z'
        end
    end
    traceout!(s_anc, [i for i = 1:n])
    p_anc = QuantumClifford._P_str(string(p_anc...))
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
    @test expect(st, P"XX") == alter_expect(st, P"XX") == 1
    @test expect(st, P"ZZ") == alter_expect(st, P"ZZ") == -1
    @test expect(st, P"XI") == alter_expect(st, P"XI") == 0
    @test expect(st, P"ZI") == alter_expect(st, P"ZI") == 0
    @test expect(st, P"II") == alter_expect(st, P"II") == 1
    @test expect(st, P"-XX") == alter_expect(st, P"-XX") == -1        
    @test expect(st, P"-iXX") == alter_expect(st, P"-iXX") == -im
end