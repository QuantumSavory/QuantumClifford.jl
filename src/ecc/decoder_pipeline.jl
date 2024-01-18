"""An abstract type for QECC syndrome decoding algorithms.

All `AbstractSyndromeDecoder` types are expected to:
- have a `parity_checks` method giving the parity checks for the code under study
- have a `decode` method that guesses error which caused the syndrome
- have an `evaluate_decoder` method which runs a full simulation but it supports only a small number of ECC protocols"""
abstract type AbstractSyndromeDecoder end

"""An abstract type mostly used by [`evaluate_decoder`](@ref) to specify in what context to evaluate an ECC."""
abstract type AbstractECCSetup end

"""A helper function that takes a parity check tableau and an `AbstractECCSetup` type and provides the circuit that needs to be simulated."""
function physical_ECC_circuit end # XXX Do not export! This might need to be refactored as we add more interesting setups!

"""Configuration for ECC evaluators that simulate the Shor-style syndrome measurement (without a flag qubit).

The simulated circuit includes:
- perfect noiseless encoding (encoding and its fault tolerance are not being studied here)
- one round of "memory noise" after the encoding but before the syndrome measurement
- perfect preparation of entangled ancillary qubits
- noisy Shor-style syndrome measurement (only two-qubit gate noise)
- noiseless "logical state measurement" (providing the comparison data when evaluating the decoder)
"""
struct ShorSyndromeECCSetup <: AbstractECCSetup
    mem_noise::Float64
    two_qubit_gate_noise::Float64
    function ShorSyndromeECCSetup(mem_noise, two_qubit_gate_noise)
        0<=mem_noise<=1 || throw(DomainError(mem_noise, "The memory noise in `ShorSyndromeECCSetup` should be between 0 and 1."))
        0<=two_qubit_gate_noise<=1 || throw(DomainError(two_qubit_gate_noise, "The two-qubit gate noise in `ShorSyndromeECCSetup` should be between 0 and 1."))
        new(mem_noise, two_qubit_gate_noise)
    end
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

    X_error = evaluate_decoder(
        d, nsamples,
        [encoding_circ..., physical_noisy_circ..., logZ_circ...],
        syndrome_bits, logZ_bits, O[length(logZ_bits)+1:end,:])
    Z_error = evaluate_decoder(
        d, nsamples,
        [preX..., encoding_circ..., physical_noisy_circ..., logX_circ...],
        syndrome_bits, logX_bits, O[1:length(logZ_bits),:])
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
    decoded = 0
    for i in 1:nsamples
        guess = decode(d, @view syndromes[i,:])
        isnothing(guess) && continue
        guess_faults = faults_submatrix * guess
        if guess_faults == @view measured_faults[i,:]
            decoded += 1
        end
    end

    return (nsamples - decoded) / nsamples
end

struct TableDecoder <: AbstractSyndromeDecoder
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
    """The lookup table corresponding to the code, slow to create"""
    lookup_table
end

function TableDecoder(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    lookup_table = create_lookup_table(H)
    fm = faults_matrix(H)
    return TableDecoder(H, n, s, k, fm, lookup_table)
end

parity_checks(d::TableDecoder) = d.H

function create_lookup_table(code::Stabilizer)
    lookup_table = Dict()
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
    return get(d.lookup_table, syndrome_sample, nothing)
end

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
end

function BeliefPropDecoder(Hx, Hz)
    c = CSS(Hx, Hz)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    println(typeof(H))
    faults_matrix = faults_matrix(H)
    log_probabs = zeros(n)
    channel_probs = fill(p_init, n)

    numchecks_X = size(Cx)[1]
    b2c_X = zeros(numchecks_X, n)
    c2b_X = zeros(numchecks_X, n)

    numchecks_Z = size(Cz)[1]
    b2c_Z = zeros(numchecks_Z, n)
    c2b_Z = zeros(numchecks_Z, n)
    err = zeros(n)

    sparse_Cx = sparse(Hx)
    sparse_CxT = sparse(Hx')
    sparse_Cz = sparse(Hz)
    sparse_CzT = sparse(Hz')
    return BeliefPropDecoder(H, faults_matrix, n, s, k, log_probabs, channel_probs, numchecks_X, b2c_X, c2b_X, numchecks_Z, b2c_Z, c2b_Z, err, sparse_Cx, sparse_CxT, sparse_Cz, sparse_CzT)
end

parity_checks(d::BeliefPropDecoder) = d.H

function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.numchecks_X]
    row_z = syndrome_sample[d.numchecks_X+1:d.numchecks_X+d.numchecks_Z]

    KguessX, success = syndrome_decode(d.sparse_Cx, d.sparse_CxT, d.row_x, d.max_iters, d.channel_probs, d.b2c_X, d.c2b_X, d.log_probabs, Base.copy(d.err))
    KguessZ, success = syndrome_decode(d.sparse_Cz, d.sparse_CzT, d.row_z, d.max_iters, d.channel_probs, d.b2c_Z, d.c2b_Z, d.log_probabs, Base.copy(d.err))
    guess = vcat(KguessZ, KguessX)
end
