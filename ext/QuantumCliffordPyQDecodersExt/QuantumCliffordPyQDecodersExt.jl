module QuantumCliffordPyQDecodersExt

using PyQDecoders: np, sps, ldpc, pm, PythonCall
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, batchdecode, parity_checks

abstract type PyBP <: AbstractSyndromeDecoder end

struct PyBeliefPropDecoder <: PyBP # TODO all these decoders have the same fields, maybe we can factor out a common type
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

struct PyBeliefPropOSDecoder <: PyBP # TODO all these decoders have the same fields, maybe we can factor out a common type
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

function PyBeliefPropDecoder(c; maxiter=nothing, bpmethod=nothing, errorrate=nothing)
    Hx = reinterpret(UInt8,collect(parity_checks_x(c)))
    Hz = reinterpret(UInt8,collect(parity_checks_z(c)))
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter=isnothing(maxiter) ? 0 : maxiter
    bpmethod ∈ (nothing, :productsum, :minsum) || error(lazy"PyBeliefPropDecoder got an unknown belief propagation method argument. `bpmethod` must be one of :productsum, :minsum.")
    bp_method = get(Dict(:productsum => "product_sum", :minsum => "minimum_sum"), bpmethod, "minimum_sum")
    isnothing(errorrate) || 0≤errorrate≤1 || error(lazy"PyBeliefPropDecoder got an invalid error rate argument. `errorrate` must be in the range [0, 1].")
    error_rate = isnothing(errorrate) ? PythonCall.Py(0.0001) : errorrate
    pyx = ldpc.BpDecoder(np.array(Hx); max_iter, bp_method, error_rate) # TODO should be sparse
    pyz = ldpc.BpDecoder(np.array(Hz); max_iter, bp_method, error_rate) # TODO should be sparse
    return PyBeliefPropDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

function PyBeliefPropOSDecoder(c; maxiter=nothing, bpmethod=nothing, errorrate=nothing, osdmethod=nothing, osdorder=0)
    Hx = reinterpret(UInt8,collect(parity_checks_x(c)))
    Hz = reinterpret(UInt8,collect(parity_checks_z(c)))
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter=isnothing(maxiter) ? 0 : maxiter
    bpmethod ∈ (nothing, :productsum, :minsum) || error(lazy"PyBeliefPropOSDecoder got an unknown belief propagation method argument. `bpmethod` must be one of :productsum, :minsum.")
    bp_method = get(Dict(:productsum => "product_sum", :minsum => "minimum_sum"), bpmethod, "minimum_sum")
    isnothing(errorrate) || 0≤errorrate≤1 || error(lazy"PyBeliefPropOSDecoder got an invalid error rate argument. `errorrate` must be in the range [0, 1].")
    error_rate = isnothing(errorrate) ? PythonCall.Py(0.0001) : errorrate
    isnothing(osdmethod) || osdmethod ∈ (:zeroorder, :exhaustive, :combinationsweep) || error(lazy"PyBeliefPropOSDecoder got an unknown OSD method argument. `osdmethod` must be one of :zeroorder, :exhaustive, :combinationsweep.")
    osd_method = get(Dict(:zeroorder => "OSD_0", :exhaustive => "OSD_E", :combinationsweep => "OSD_CS"), osdmethod, 0)
    0≥osdorder || error(lazy"PyBeliefPropOSDecoder got an invalid OSD order argument. `osdorder` must be ≥0.")
    osd_order = osdorder
    pyx = ldpc.BpOsdDecoder(np.array(Hx); max_iter, bp_method, error_rate, osd_method, osd_order) # TODO should be sparse
    pyz = ldpc.BpOsdDecoder(np.array(Hz); max_iter, bp_method, error_rate, osd_method, osd_order) # TODO should be sparse
    return PyBeliefPropOSDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::PyBP) = d.H

function decode(d::PyBP, syndrome_sample)
    row_x = @view syndrome_sample[1:d.nx]
    row_z = @view syndrome_sample[d.nx+1:end]
    guess_z_errors = PythonCall.PyArray(d.pyx.decode(np.array(row_x)))
    guess_x_errors = PythonCall.PyArray(d.pyz.decode(np.array(row_z)))
    vcat(guess_x_errors, guess_z_errors)
end

struct PyMatchingDecoder <: AbstractSyndromeDecoder # TODO all these decoders have the same fields, maybe we can factor out a common type
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
    row_x = @view syndrome_sample[1:d.nx]
    row_z = @view syndrome_sample[d.nx+1:end]
    guess_z_errors = PythonCall.PyArray(d.pyx.decode(row_x))
    guess_x_errors = PythonCall.PyArray(d.pyz.decode(row_z))
    vcat(guess_x_errors, guess_z_errors)
end

function batchdecode(d::PyMatchingDecoder, syndrome_samples)
    row_x = @view syndrome_samples[:,1:d.nx]
    row_z = @view syndrome_samples[:,d.nx+1:end]
    guess_z_errors = PythonCall.PyArray(d.pyx.decode_batch(row_x))
    guess_x_errors = PythonCall.PyArray(d.pyz.decode_batch(row_z))
    n_cols_x = size(guess_x_errors, 2)
    hcat(guess_x_errors, guess_z_errors)
end

end
