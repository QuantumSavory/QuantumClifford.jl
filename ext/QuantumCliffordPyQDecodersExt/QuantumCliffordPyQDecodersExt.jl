module QuantumCliffordPyQDecodersExt

using PyQDecoders: np, sps, ldpc, pm
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, parity_checks


struct PyBeliefPropDecoder <: AbstractSyndromeDecoder
    code
    H
    Hx
    Hz
    nx
    nz
    faults_matrix
    pyx
    pyz
end

function PyBeliefPropDecoder(c)
    Hx = parity_checks_x(c) |> collect
    Hz = parity_checks_z(c) |> collect
    H = parity_checks(c)
    fm = faults_matrix(c)
    pyx = ldpc.bp_decoder(Hx)
    pyz = ldpc.bp_decoder(Hz)
    return PyBeliefPropDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::PyBeliefPropDecoder) = d.H

function decode(d::PyBeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.nx] # TODO This copy is costly!
    row_z = syndrome_sample[d.nx+1:end]
    guess_x = d.pyx.decode(row_x)
    guess_z = d.pyz.decode(row_z)
    vcat(guess_z, guess_x)
end


struct PyMatchingDecoder <: AbstractSyndromeDecoder
    code
    H
    Hx
    Hz
    nx
    nz
    faults_matrix
    pyx
    pyz
end

function PyMatchingDecoder(c; weights=nothing)
    Hx = parity_checks_x(c) |> collect # TODO keep these sparse
    Hz = parity_checks_z(c) |> collect
    H = parity_checks(c)
    fm = faults_matrix(c)
    if isnothing(weights)
        pyx = pm.Matching.from_check_matrix(Hx)
        pyz = pm.Matching.from_check_matrix(Hz)
    else
        pyx = pm.Matching.from_check_matrix(Hx, weights=weights)
        pyz = pm.Matching.from_check_matrix(Hz, weights=weights)
    end
    return PyMatchingDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::PyMatchingDecoder) = d.H

function decode(d::PyMatchingDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.nx] # TODO This copy is costly!
    row_z = syndrome_sample[d.nx+1:end]
    guess_z_errors = d.pyx.decode(row_x)
    guess_x_errors = d.pyz.decode(row_z)
    vcat(guess_x_errors, guess_z_errors)
end

end
