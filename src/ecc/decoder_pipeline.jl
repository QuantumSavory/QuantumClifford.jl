"""An abstract type for QECC syndrome decoding algorithms.

All `AbstractSyndromeDecoder` types are expected to:
- have a `parity_checks` method giving the parity checks for the code under study
- have a `decode` method that guesses error which caused the syndrome
- have an `evaluate_decoder` method which runs a full simulation but it supports only a small number of ECC protocols"""
abstract type AbstractSyndromeDecoder end

"""Decode a syndrome using a given decoding algorithm."""
function decode end

"""Decode a batch of syndromes using a given decoding algorithm."""
function batchdecode(d::AbstractSyndromeDecoder, syndrome_samples)
    H = parity_checks(d)
    s, n = size(H)
    samples, _s = size(syndrome_samples)
    s == _s || throw(ArgumentError(lazy"The syndromes given to `batchdecode` have the wrong dimensions. The syndrome length is $(_s) while it should be $(s)"))
    results = falses(samples, 2n)
    for (i,syndrome_sample) in enumerate(eachrow(syndrome_samples))
        guess = decode(d, syndrome_sample)# TODO use `decode!`
        isnothing(guess) || (results[i,:] = guess)
    end
    return results
end


function add_two_qubit_gate_noise(g, gate_error)
    return ()
end

"""Applies gate_error to a given two-qubit gate g."""
function add_two_qubit_gate_noise(g::AbstractTwoQubitOperator, gate_error)
    return (PauliError(affectedqubits(g), gate_error), )
end


"""An abstract type mostly used by [`evaluate_decoder`](@ref) to specify in what context to evaluate an ECC."""
abstract type AbstractECCSetup end

"""A helper function that takes a parity check tableau and an `AbstractECCSetup` type and provides the circuit that needs to be simulated."""
function physical_ECC_circuit end # XXX Do not export! This might need to be refactored as we add more interesting setups!

"""Configuration for ECC evaluator that does not simulate any ECC circuits, rather it simply checks the commutation of the parity check and the Pauli error.

This is much faster than any other simulation method, but it is incapable of noisy-circuit simulations and thus useless for fault-tolerance studies.

See also: [`NaiveSyndromeECCSetup`](@ref), [`ShorSyndromeECCSetup`](@ref)"""
struct CommutationCheckECCSetup <: AbstractECCSetup
    xz_noise::Float64
    function CommutationCheckECCSetup(xz_noise)
        0<=xz_noise<=1 || throw(DomainError(xz_noise, "The independent X/Z memory noise in `CommutationCheckECCSetup` should be between 0 and 1."))
        new(xz_noise)
    end
end

"""Configuration for ECC evaluator that runs the simplest syndrome measurement circuit.

The circuit is being simulated (as opposed to doing only a quick commutation check).
This circuit would give poor performance if there is non-zero gate noise.

See also: [`CommutationCheckECCSetup`](@ref), [`ShorSyndromeECCSetup`](@ref)"""
struct NaiveSyndromeECCSetup <: AbstractECCSetup
    mem_noise::Float64
    two_qubit_gate_noise::Float64
    function NaiveSyndromeECCSetup(mem_noise, two_qubit_gate_noise)
        0<=mem_noise<=1 || throw(DomainError(mem_noise, "The memory noise in `NaiveSyndromeECCSetup` should be between 0 and 1."))
        0<=two_qubit_gate_noise<=1 || throw(DomainError(two_qubit_gate_noise, "The two-qubit gate noise in `NaiveSyndromeECCSetup` should be between 0 and 1."))
        new(mem_noise, two_qubit_gate_noise)
    end
end

"""Configuration for ECC evaluators that simulate the Shor-style syndrome measurement (without a flag qubit).

The simulated circuit includes:
- perfect noiseless encoding (encoding and its fault tolerance are not being studied here)
- one round of "memory noise" after the encoding but before the syndrome measurement
- perfect preparation of entangled ancillary qubits
- noisy Shor-style syndrome measurement (only two-qubit gate noise)
- noiseless "logical state measurement" (providing the comparison data when evaluating the decoder)

See also: [`CommutationCheckECCSetup`](@ref), [`NaiveSyndromeECCSetup`](@ref)"""
struct ShorSyndromeECCSetup <: AbstractECCSetup
    mem_noise::Float64
    two_qubit_gate_noise::Float64
    function ShorSyndromeECCSetup(mem_noise, two_qubit_gate_noise)
        0<=mem_noise<=1 || throw(DomainError(mem_noise, "The memory noise in `ShorSyndromeECCSetup` should be between 0 and 1."))
        0<=two_qubit_gate_noise<=1 || throw(DomainError(two_qubit_gate_noise, "The two-qubit gate noise in `ShorSyndromeECCSetup` should be between 0 and 1."))
        new(mem_noise, two_qubit_gate_noise)
    end
end

function physical_ECC_circuit(H, setup::NaiveSyndromeECCSetup)
    syndrome_circ, n_anc, syndrome_bits = naive_syndrome_circuit(H)
    new_circuit = []

    for op in syndrome_circ
        push!(new_circuit, op)
        
        for noise_op in add_two_qubit_gate_noise(op, setup.two_qubit_gate_noise)
            push!(new_circuit, noise_op)
        end
    end

    mem_error_circ = [PauliError(i, setup.mem_noise) for i in 1:nqubits(H)]


    circ = vcat(mem_error_circ, new_circuit)
    return circ, syndrome_bits, n_anc
end


function physical_ECC_circuit(H, setup::ShorSyndromeECCSetup)
    prep_anc, syndrome_circ, n_anc, syndrome_bits = shor_syndrome_circuit(H)
    noisy_syndrome_circ = syndrome_circ # add_two_qubit_gate_noise(syndrome_circ, gate_error)
    mem_error_circ = [PauliError(i, setup.mem_noise) for i in 1:nqubits(H)];
    circ = [prep_anc..., mem_error_circ..., noisy_syndrome_circ...]
    circ, syndrome_bits, n_anc
end

"""Evaluate the performance of a given decoder (e.g. [`TableDecoder`](@ref)) and a given style of running an ECC code (e.g. [`ShorSyndromeECCSetup`](@ref))"""
function evaluate_decoder(d::AbstractSyndromeDecoder, setup::AbstractECCSetup, nsamples::Int)
    H = parity_checks(d)
    n = code_n(H)
    k = code_k(H)
    O = faults_matrix(H)

    physical_noisy_circ, syndrome_bits, n_anc = physical_ECC_circuit(H, setup)
    encoding_circ = naive_encoding_circuit(H)
    preX = [sHadamard(i) for i in n-k+1:n]

    mdH = MixedDestabilizer(H)
    logX_circ, _, logX_bits = naive_syndrome_circuit(logicalxview(mdH), n_anc+1, last(syndrome_bits)+1)
    logZ_circ, _, logZ_bits = naive_syndrome_circuit(logicalzview(mdH), n_anc+1, last(syndrome_bits)+1)

    # Evaluate the probability for X logical error (the Z-observable part of the faults matrix is used)
    X_error = evaluate_decoder(
        d, nsamples,
        [encoding_circ..., physical_noisy_circ..., logZ_circ...],
        syndrome_bits, logZ_bits, O[end÷2+1:end,:])
    # Evaluate the probability for Z logical error (the X-observable part of the faults matrix is used)
    Z_error = evaluate_decoder(
        d, nsamples,
        [preX..., encoding_circ..., physical_noisy_circ..., logX_circ...],
        syndrome_bits, logX_bits, O[1:end÷2,:])
    return (X_error, Z_error)
end

"""Evaluate the performance of an error-correcting circuit.

This method requires you give the circuit that performs both syndrome measurements and (probably noiseless) logical state measurements.
The faults matrix that translates an error vector into corresponding logical errors is necessary as well.

This is a relatively barebones method that assumes the user prepares necessary circuits, etc.
It is a method that is used internally by more user-frienly methods providing automatic conversion of codes and noise models
to the necessary noisy circuits.
"""
function evaluate_decoder(d::AbstractSyndromeDecoder, nsamples, circuit, syndrome_bits, logical_bits, faults_submatrix)
    frames = pftrajectories(circuit;trajectories=nsamples,threads=true)

    syndromes = @view pfmeasurements(frames)[:, syndrome_bits]
    measured_faults = @view pfmeasurements(frames)[:, logical_bits]
    guesses = batchdecode(d, syndromes)
    evaluate_guesses(measured_faults, guesses, faults_submatrix)
end

function evaluate_guesses(measured_faults, guesses, faults_matrix)
    nsamples = size(guesses, 1)
    guess_faults = (faults_matrix * guesses') .% 2 # TODO this can be faster and non-allocating by turning it into a loop
    decoded = 0
    for i in 1:nsamples # TODO this can be faster and non-allocating by having the loop and the matrix multiplication on the line above work together and not store anything
        (@view guess_faults[:,i]) == (@view measured_faults[i,:]) && (decoded += 1)
    end
    return (nsamples - decoded) / nsamples
end

function evaluate_decoder(d::AbstractSyndromeDecoder, setup::CommutationCheckECCSetup, nsamples::Int)
    H = parity_checks(d)
    fm = faults_matrix(H)
    fmtab = QuantumClifford.Tableau(fm[:,end÷2+1:end],fm[:,1:end÷2]) # TODO there should be a special method for this
    k, s = size(fm)
    n = nqubits(H)
    err = zero(PauliOperator, n)
    syndromes = zeros(Bool, nsamples, length(H)) # TODO will this be better and faster if we use bitmatrices
    measured_faults = zeros(UInt8, nsamples, k)
    for i in 1:nsamples
        err = random_pauli!(err, setup.xz_noise, nophase=true)
        comm!((@view syndromes[i,:]), H,err)
        comm!((@view measured_faults[i,:]), fmtab,err)
    end
    measured_faults .%= 2
    guesses = batchdecode(d, syndromes)
    evaluate_guesses(measured_faults, guesses, fm)
end

"""A simple look-up table decoder for error correcting codes.

The lookup table contains only weight=1 errors, thus it is small,
but at best it provides only for distance=3 decoding.

The size of the lookup table would grow exponentially quickly for higher distances."""
struct TableDecoder <: AbstractSyndromeDecoder
    """Stabilizer tableau defining the code"""
    H
    """Faults matrix corresponding to the code"""
    faults_matrix
    """The number of qubits in the code"""
    n::Int
    """The depth of the code"""
    s::Int
    """The number of encoded qubits"""
    k::Int
    """The lookup table corresponding to the code, slow to create"""
    lookup_table::Dict{Vector{Bool},Vector{Bool}}
    lookup_buffer::Vector{Bool}
    TableDecoder(H, faults_matrix, n, s, k, lookup_table) = new(H, faults_matrix, n, s, k, lookup_table, fill(false, s))
end

function TableDecoder(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    lookup_table = create_lookup_table(H)
    fm = faults_matrix(H)
    return TableDecoder(H, fm, n, s, k, lookup_table)
end

parity_checks(d::TableDecoder) = d.H

function create_lookup_table(code::Stabilizer)
    lookup_table = Dict{Vector{Bool},Vector{Bool}}()
    constraints, qubits = size(code)
    # In the case of no errors
    lookup_table[ zeros(UInt8, constraints) ] = stab_to_gf2(zero(PauliOperator, qubits))
    # In the case of single bit errors
    for bit_to_be_flipped in 1:qubits
        for error_type in [single_x, single_y, single_z]
            # Generate e⃗
            error = error_type(qubits, bit_to_be_flipped)
            # Calculate s⃗
            # (check which stabilizer rows do not commute with the Pauli error)
            syndrome = comm(error, code)
            # Store s⃗ → e⃗
            lookup_table[syndrome] = stab_to_gf2(error)
        end
    end
    lookup_table
end;

function decode(d::TableDecoder, syndrome_sample)
    d.lookup_buffer .= syndrome_sample # TODO have this work without data copying, by supporting the correct types, especially in the batch decode case
    return get(d.lookup_table, d.lookup_buffer, nothing)
end

# From extensions:

"""A simple Belief Propagation decoder built around tools from `LDPCDecoders.jl`."""
function BeliefPropDecoder(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordLDPCDecodersExt)
    if isnothing(ext)
        throw("The `BeliefPropDecoder` depends on the package `LDPCDecoders` but you have not installed or imported `LDPCDecoders` yet. Immediately after you import `LDPCDecoders`, the `BeliefPropDecoder` will be available.")
    end
    return ext.BeliefPropDecoder(args...; kwargs...)
end

"""An Iterative Bitflip decoder built around tools from `LDPCDecoders.jl`."""
function BitFlipDecoder(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordLDPCDecodersExt)
    if isnothing(ext)
        throw("The `BitFlipDecoder` depends on the package `LDPCDecoders` but you have not installed or imported `LDPCDecoders` yet. Immediately after you import `LDPCDecoders`, the `BitFlipDecoder` will be available.")
    end
    return ext.BitFlipDecoder(args...; kwargs...)
end


"""A Belief Propagation decoder built around tools from the python package `ldpc` available from the julia package `PyQDecoders.jl`."""
function PyBeliefPropDecoder(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordPyQDecodersExt)
    if isnothing(ext)
        throw("The `PyBeliefPropDecoder` depends on the package `PyQDecoders` but you have not installed or imported `PyQDecoders` yet. Immediately after you import `PyQDecoders`, the `PyBeliefPropDecoder` will be available.")
    end
    return ext.PyBeliefPropDecoder(args...; kwargs...)
end

"""A Belief Propagation decoder with ordered statistics decoding, built around tools from the python package `ldpc` available from the julia package `PyQDecoders.jl`."""
function PyBeliefPropOSDecoder(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordPyQDecodersExt)
    if isnothing(ext)
        throw("The `PyBeliefPropOSDecoder` depends on the package `PyQDecoders` but you have not installed or imported `PyQDecoders` yet. Immediately after you import `PyQDecoders`, the `PyBeliefPropOSDecoder` will be available.")
    end
    return ext.PyBeliefPropOSDecoder(args...; kwargs...)
end

"""A perfect matching decoder built around tools from the python package `pymatching` available from the julia package `PyQDecoders.jl`."""
function PyMatchingDecoder(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordPyQDecodersExt)
    if isnothing(ext)
        throw("The `PyMatchingDecoder` depends on the package `PyQDecoders` but you have not installed or imported `PyMatchingDecoder` yet. Immediately after you import `PyQDecoders`, the `PyMatchingDecoder` will be available.")
    end
    return ext.PyMatchingDecoder(args...; kwargs...)
end
