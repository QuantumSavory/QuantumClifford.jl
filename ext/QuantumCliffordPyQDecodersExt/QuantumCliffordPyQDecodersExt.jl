module QuantumCliffordPyQDecodersExt

using PyQDecoders: np, sps, ldpc, pm, PythonCall
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, parity_checks

abstract type PyBP <: AbstractSyndromeDecoder end

struct PyBeliefPropDecoder <: PyBP
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

struct PyBeliefPropOSDecoder <: PyBP
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

function PyBeliefPropDecoder(c; maxiter=nothing)
    Hx = parity_checks_x(c) |> collect # TODO should be sparse
    Hz = parity_checks_z(c) |> collect # TODO should be sparse
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter=isnothing(maxiter) ? 0 : maxiter
    pyx = ldpc.bp_decoder(np.array(Hx); max_iter) # TODO should be sparse
    pyz = ldpc.bp_decoder(np.array(Hz); max_iter) # TODO should be sparse
    return PyBeliefPropDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

function PyBeliefPropOSDecoder(c; maxiter=nothing)
    Hx = parity_checks_x(c) |> collect # TODO should be sparse
    Hz = parity_checks_z(c) |> collect # TODO should be sparse
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter=isnothing(maxiter) ? 0 : maxiter
    pyx = ldpc.bposd_decoder(np.array(Hx); max_iter) # TODO should be sparse
    pyz = ldpc.bposd_decoder(np.array(Hz); max_iter) # TODO should be sparse
    return PyBeliefPropOSDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::PyBP) = d.H

function decode(d::PyBP, syndrome_sample)
    row_x = syndrome_sample[1:d.nx] # TODO These copies and indirections might be costly!
    row_z = syndrome_sample[d.nx+1:end]
    guess_x = PythonCall.PyArray(d.pyx.decode(np.array(row_x)))
    guess_z = PythonCall.PyArray(d.pyz.decode(np.array(row_z)))
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
    guess_z_errors = PythonCall.PyArray(d.pyx.decode(row_x))
    guess_x_errors = PythonCall.PyArray(d.pyz.decode(row_z))
    vcat(guess_x_errors, guess_z_errors)
end

end
