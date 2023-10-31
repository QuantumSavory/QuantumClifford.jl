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
    """Empty array to hold temporary values in belief decoding"""
    log_probabs
    """Error probabilities of each channel"""
    channel_probs
    """Number of X checks, used to get the syndrome corresponding to the X checks"""
    numchecks_X
    """Empty matrix used to hold error probabilities for the X channels"""
    b2c_X
    """Empty matrix used to temporary belief propagation values"""
    c2b_X
    """Number of X checks, used to get the syndrome corresponding to the X checks"""
    numchecks_Z
    """Empty matrix used to hold error probabilities for the Z channels"""
    b2c_Z
    """Empty matrix used to temporary belief propagation values"""
    c2b_Z
    """The measured error syndrome"""
    err
    """Sparse array of Cx matrix"""
    sparse_Cx
    """Sparse array of the transpose of the Cx matrix"""
    sparse_CxT
    """Sparse array of Cz matrix"""
    sparse_Cz
    """Sparse array of the transpose of the Cx matrix"""
    sparse_CzT
    """Maximum number of iterations before giving up"""
    max_iters
end

function BeliefPropDecoder(c, p_init=0, max_iters=10)
    Hx = parity_checks_x(c)
    Hz = parity_checks_z(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    fm = faults_matrix(H)
    log_probabs = zeros(n)
    channel_probs = fill(p_init, n)

    numchecks_X = size(Hx, 1)
    b2c_X = zeros(numchecks_X, n)
    c2b_X = zeros(numchecks_X, n)

    numchecks_Z = size(Hz, 1)
    b2c_Z = zeros(numchecks_Z, n)
    c2b_Z = zeros(numchecks_Z, n)
    err = zeros(n)

    sparse_Cx = sparse(Hx)
    sparse_CxT = sparse(Hx')
    sparse_Cz = sparse(Hz)
    sparse_CzT = sparse(Hz')
    return BeliefPropDecoder(H, fm, n, s, k, log_probabs, channel_probs, numchecks_X, b2c_X, c2b_X, numchecks_Z, b2c_Z, c2b_Z, err, sparse_Cx, sparse_CxT, sparse_Cz, sparse_CzT, max_iters)
end

parity_checks(d::BeliefPropDecoder) = d.H

function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.numchecks_X]
    row_z = syndrome_sample[d.numchecks_X+1:d.numchecks_X+d.numchecks_Z]

    KguessX, success = LDPCDecoders.syndrome_decode(d.sparse_Cx, d.sparse_CxT, row_x, d.max_iters, d.channel_probs, d.b2c_X, d.c2b_X, d.log_probabs, Base.copy(d.err))
    KguessZ, success = LDPCDecoders.syndrome_decode(d.sparse_Cz, d.sparse_CzT, row_z, d.max_iters, d.channel_probs, d.b2c_Z, d.c2b_Z, d.log_probabs, Base.copy(d.err))
    guess = vcat(KguessZ, KguessX)
end

end
