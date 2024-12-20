module QuantumCliffordPyQDecodersExt

using PyQDecoders: np, sps, ldpc, pm, PythonCall
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, batchdecode, parity_checks

abstract type PyBP <: AbstractSyndromeDecoder end

# A common structure to hold shared fields for different decoders
mutable struct GenericLDPCDecoder{D} <: PyBP
    code
    H
    Hx
    Hz
    nx
    nz
    faults_matrix
    pyx::D
    pyz::D
end

# Common function to initialize PyBeliefPropDecoder or PyBeliefPropOSDecoder
function initialize_decoder(c, maxiter, bpmethod, errorrate, osdmethod, osdorder, decoder_type)
    Hx = parity_checks_x(c) |> collect # TODO keep these sparse
    Hz = parity_checks_z(c) |> collect # TODO keep these sparse
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter = isnothing(maxiter) ? 0 : maxiter
    bpmethod ∈ (nothing, :productsum, :minsum, :minsumlog) || error("Unknown bpmethod")
    bp_method = get(Dict(:productsum => 0, :minsum => 1, :minsumlog => 2), bpmethod, 0)
    error_rate = isnothing(errorrate) ? PythonCall.Py(nothing) : errorrate
    osd_method = isnothing(osdmethod) ? "osd0" : osdmethod
    if decoder_type == :beliefprop
        pyx = ldpc.bp_decoder(np.array(Hx); max_iter, bp_method, error_rate) # TODO keep these sparse
        pyz = ldpc.bp_decoder(np.array(Hz); max_iter, bp_method, error_rate) # TODO keep these sparse
    elseif decoder_type == :beliefprop_os
        pyx = ldpc.bposd_decoder(np.array(Hx); max_iter, bp_method, error_rate, osd_method, osdorder)
        pyz = ldpc.bposd_decoder(np.array(Hz); max_iter, bp_method, error_rate, osd_method, osdorder)
    else
        error("Unknown decoder type.")
    end
    return GenericPyLDPCDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

function PyBeliefPropDecoder(c; maxiter=nothing, bpmethod=nothing, errorrate=nothing)
    return initialize_decoder(c, maxiter, bpmethod, errorrate, nothing, 0, :beliefprop)
end

function PyBeliefPropOSDecoder(c; maxiter=nothing, bpmethod=nothing, errorrate=nothing, osdmethod=nothing, osdorder=0)
    return initialize_decoder(c, maxiter, bpmethod, errorrate, osdmethod, osdorder, :beliefprop_os)
end

function PyMatchingDecoder(c; weights=nothing)
    Hx = parity_checks_x(c) |> collect # TODO keep these sparse
    Hz = parity_checks_z(c) |> collect # TODO keep these sparse
    H = parity_checks(c)
    fm = faults_matrix(c)
    if isnothing(weights)
        pyx = pm.Matching.from_check_matrix(Hx)
        pyz = pm.Matching.from_check_matrix(Hz)
    else
        pyx = pm.Matching.from_check_matrix(Hx, weights=weights)
        pyz = pm.Matching.from_check_matrix(Hz, weights=weights)
    end
    return GenericPyLDPCDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::GenericPyLDPCDecoder) = d.H

function decode(d::GenericPyLDPCDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.nx]
    row_z = @view syndrome_sample[d.nx+1:end]
    guess_z_errors = PythonCall.PyArray(d.pyx.decode(np.array(row_x)))
    guess_x_errors = PythonCall.PyArray(d.pyz.decode(np.array(row_z)))
    return vcat(guess_x_errors, guess_z_errors)
end

function batchdecode(d::GenericPyLDPCDecoder, syndrome_samples)
    row_x = @view syndrome_samples[:,1:d.nx]
    row_z = @view syndrome_samples[:,d.nx+1:end]
    guess_z_errors = PythonCall.PyArray(d.pyx.decode_batch(row_x))
    guess_x_errors = PythonCall.PyArray(d.pyz.decode_batch(row_z))
    return hcat(guess_x_errors, guess_z_errors)
end

end
