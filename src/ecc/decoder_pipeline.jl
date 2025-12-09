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

function add_two_qubit_gate_noise(g, gate_error)
    return ()
end

"""Applies gate_error to a given two-qubit gate g."""
function add_two_qubit_gate_noise(g::AbstractTwoQubitOperator, gate_error)
    qubits = affectedqubits(g)
    return (PauliError(qubits, gate_error), )
end

function physical_ECC_circuit(H, setup::NaiveSyndromeECCSetup)
    syndrome_circ, n_anc, syndrome_bits = naive_syndrome_circuit(H)
    noisy_syndrome_circ = []

    for op in syndrome_circ
        push!(noisy_syndrome_circ, op)
        for noise_op in add_two_qubit_gate_noise(op, setup.two_qubit_gate_noise)
            push!(noisy_syndrome_circ, noise_op)
        end
    end

    mem_error_circ = [PauliError(i, setup.mem_noise) for i in 1:nqubits(H)]
    circ = vcat(mem_error_circ, noisy_syndrome_circ)
    return circ, syndrome_bits, n_anc
end


function physical_ECC_circuit(H, setup::ShorSyndromeECCSetup)
    prep_anc, syndrome_circ, n_anc, syndrome_bits = shor_syndrome_circuit(H)

    noisy_syndrome_circ = []
    for op in syndrome_circ
        push!(noisy_syndrome_circ, op)
        for noise_op in add_two_qubit_gate_noise(op, setup.two_qubit_gate_noise)
            push!(noisy_syndrome_circ, noise_op)
        end
    end

    mem_error_circ = [PauliError(i, setup.mem_noise) for i in 1:nqubits(H)]

    circ = vcat(prep_anc, mem_error_circ, noisy_syndrome_circ)
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
    preX = sHadamard[sHadamard(i) for i in n-k+1:n]

    mdH = MixedDestabilizer(H)
    logX_circ, _, logX_bits = naive_syndrome_circuit(logicalxview(mdH), n_anc+1, last(syndrome_bits)+1)
    logZ_circ, _, logZ_bits = naive_syndrome_circuit(logicalzview(mdH), n_anc+1, last(syndrome_bits)+1)

    # Evaluate the probability for X logical error (the Z-observable part of the faults matrix is used)
    X_error = evaluate_decoder(
        d, nsamples,
        vcat(encoding_circ, physical_noisy_circ, logZ_circ),
        syndrome_bits, logZ_bits, O[end÷2+1:end,:])
    # Evaluate the probability for Z logical error (the X-observable part of the faults matrix is used)
    Z_error = evaluate_decoder(
        d, nsamples,
        vcat(preX, encoding_circ, physical_noisy_circ, logX_circ),
        syndrome_bits, logX_bits, O[1:end÷2,:])
    return (X_error, Z_error)
end

"""Evaluate the performance of an error-correcting circuit.

This method requires you give the circuit that performs both syndrome measurements and (probably noiseless) logical state measurements.
The faults matrix that translates an error vector into corresponding logical errors is necessary as well.

This is a relatively barebones method that assumes the user prepares necessary circuits, etc.
It is a method that is used internally by more user-friendly methods providing automatic conversion of codes and noise models
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
    fails = 0
    for i in 1:nsamples
        for j in 1:size(faults_matrix, 1)
            sum_mod = 0
            @inbounds @simd for k in 1:size(faults_matrix, 2)
                sum_mod += faults_matrix[j, k] * guesses[i, k]
            end
            sum_mod %= 2
            if sum_mod != measured_faults[i, j]
                fails += 1
                break
            end
        end
    end
    return fails / nsamples
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
    """The maximum weight of errors in the lookup table"""
    error_weight::Int
    """The lookup table corresponding to the code, slow to create"""
    lookup_table::Dict{Vector{Bool},Vector{Bool}}
    lookup_buffer::Vector{Bool}
    TableDecoder(H, faults_matrix, n, s, k, error_weight, lookup_table) = new(H, faults_matrix, n, s, k, error_weight, lookup_table, fill(false, s))
end

function TableDecoder(c; error_weight=1)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    lookup_table = create_lookup_table(H; error_weight)
    fm = faults_matrix(H)
    return TableDecoder(H, fm, n, s, k, error_weight, lookup_table)
end

parity_checks(d::TableDecoder) = d.H

function create_lookup_table(code::Stabilizer; error_weight=1)
    lookup_table = Dict{Vector{Bool},Vector{Bool}}()
    constraints, qubits = size(code)
    # Process errors from highest weight to lowest
    # so that lower-weight errors overwrite higher-weight ones
    # (lower-weight errors are more probable)
    for w in error_weight:-1:1
        # Iterate over all combinations of w qubit positions
        for positions in combinations(1:qubits, w)
            # For each position, there are 3 choices: X, Y, Z
            # We iterate over all 3^w combinations using (x_bit, z_bit) tuples
            for pauli_types in Iterators.product(ntuple(_ -> ((true,false), (true,true), (false,true)), w)...)
                error = zero(PauliOperator, qubits)
                for (i, pos) in enumerate(positions)
                    error[pos] = pauli_types[i]
                end
                # Calculate syndrome (check which stabilizer rows do not commute with the error)
                syndrome = comm(error, code)
                # Store s⃗ → e⃗
                lookup_table[syndrome] = stab_to_gf2(error)
            end
        end
    end
    # In the case of no errors
    lookup_table[zeros(UInt8, constraints)] = stab_to_gf2(zero(PauliOperator, qubits))
    lookup_table
end;

function decode(d::TableDecoder, syndrome_sample)
    d.lookup_buffer .= syndrome_sample # TODO have this work without data copying, by supporting the correct types, especially in the batch decode case
    return get(d.lookup_table, d.lookup_buffer, nothing)
end

"""A simple look-up table decoder for classical codes.

Similar to [`TableDecoder`](@ref) but works with a classical parity check matrix
instead of a stabilizer tableau."""
struct ClassicalTableDecoder <: AbstractSyndromeDecoder
    """Parity check matrix defining the code"""
    H::Matrix{Bool}
    """The number of bits in the code"""
    n::Int
    """The number of parity checks"""
    s::Int
    """The maximum weight of errors in the lookup table"""
    error_weight::Int
    """The lookup table corresponding to the code"""
    lookup_table::Dict{Vector{Bool},Vector{Bool}}
    lookup_buffer::Vector{Bool}
    ClassicalTableDecoder(H, n, s, error_weight, lookup_table) = new(H, n, s, error_weight, lookup_table, fill(false, s))
end

function ClassicalTableDecoder(H::Matrix{Bool}; error_weight=1)
    s, n = size(H)
    lookup_table = create_lookup_table(H; error_weight)
    return ClassicalTableDecoder(H, n, s, error_weight, lookup_table)
end

parity_checks(d::ClassicalTableDecoder) = d.H

function create_lookup_table(H::Matrix{Bool}; error_weight=1) # TODO there is inefficient casting between Bool vectors and bitvectors here
    lookup_table = Dict{Vector{Bool},Vector{Bool}}()
    s, n = size(H)
    # Process errors from highest weight to lowest
    # so that lower-weight errors overwrite higher-weight ones
    # (lower-weight errors are more probable)
    for w in error_weight:-1:1
        for positions in combinations(1:n, w)
            error = falses(n)
            for pos in positions
                error[pos] = true
            end
            # Calculate syndrome: s = H * e (mod 2)
            syndrome = Bool[(sum(H[row, pos] for pos in positions) % 2) == 1 for row in 1:s]
            lookup_table[syndrome] = error
        end
    end
    # In the case of no errors
    lookup_table[falses(s)] = falses(n)
    lookup_table
end

function decode(d::ClassicalTableDecoder, syndrome_sample)
    d.lookup_buffer .= syndrome_sample
    return get(d.lookup_table, d.lookup_buffer, nothing)
end

"""A look-up table decoder for CSS codes.

Uses two [`ClassicalTableDecoder`](@ref) instances internally to decode
the X and Z errors separately."""
struct CSSTableDecoder <: AbstractSyndromeDecoder
    """Stabilizer tableau defining the code"""
    H
    """Faults matrix corresponding to the code"""
    faults_matrix
    """The number of qubits in the code"""
    n::Int
    """The number of parity checks"""
    s::Int
    """The number of encoded qubits"""
    k::Int
    """The number of X checks"""
    cx::Int
    """The number of Z checks"""
    cz::Int
    """The maximum weight of errors in the lookup table"""
    error_weight::Int
    """Classical decoder for X errors (decodes Z syndrome)"""
    tabledecoderx::ClassicalTableDecoder
    """Classical decoder for Z errors (decodes X syndrome)"""
    tabledecoderz::ClassicalTableDecoder
end

function CSSTableDecoder(c; error_weight=1)
    Hx = parity_matrix_x(c)
    Hz = parity_matrix_z(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(copy(H), ranks=true)
    k = n - r
    cx = size(Hx, 1)
    cz = size(Hz, 1)
    fm = faults_matrix(H)
    tabledecoderx = ClassicalTableDecoder(Matrix{Bool}(Hx); error_weight)
    tabledecoderz = ClassicalTableDecoder(Matrix{Bool}(Hz); error_weight)
    return CSSTableDecoder(H, fm, n, s, k, cx, cz, error_weight, tabledecoderx, tabledecoderz)
end

parity_checks(d::CSSTableDecoder) = d.H

function decode(d::CSSTableDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.cx]
    row_z = @view syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z = decode(d.tabledecoderx, row_x)
    guess_x = decode(d.tabledecoderz, row_z)
    return isnothing(guess_x) || isnothing(guess_z) ? nothing : vcat(guess_x, guess_z)
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
