module QuantumCliffordLDPCDecodersExt

import LDPCDecoders
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, parity_checks

struct BeliefPropDecoder <: AbstractSyndromeDecoder
    """Stabilizer tableau defining the code"""
    H
    """Faults matrix corresponding to the code"""
    faults_matrix
    """The number of qubits in the code"""
    n
    """The depth of the code"""
    s
    """The number of encoded qubits"""
    k
    """Number of X checks"""
    cx
    """Number of X checks"""
    cz
    """The classical BP decoder for Hx"""
    bpdecoderx
    """The classical BP decoder for Hz"""
    bpdecoderz
end

function BeliefPropDecoder(c, p_init=0, max_iters=10)
    Hx = parity_checks_x(c)
    Hz = parity_checks_z(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    cx = size(Hx, 1)
    cz = size(Hx, 1)
    fm = faults_matrix(H)

    bpx = LDPCDecoders.BeliefPropagationDecoder(Hx, p_init, max_iters)
    bpz = LDPCDecoders.BeliefPropagationDecoder(Hz, p_init, max_iters)

    return BeliefPropDecoder(H, fm, n, s, k, cx, cz, bpx, bpz)
end

parity_checks(d::BeliefPropDecoder) = d.H

function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.cx]
    row_z = syndrome_sample[d.cx+1:d.cx+d.cz]

    guess_x = falses(d.n)
    guess_z = falses(d.n)

    success = LDPCDecoders.decode!(d.bpx, row_x, guess_x)
    success = LDPCDecoders.decode!(d.bpz, row_z, guess_z)
    return vcat(guess_z, guess_x)
end

end
