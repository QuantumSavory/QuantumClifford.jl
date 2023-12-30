abstract type AbstractSyndromeDecoder end

function evaluate_decoder(d::AbstractSyndromeDecoder, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func)
    pre_X = [sHadamard(i) for i in n-k+1:n]
    X_error = evaluate_classical_decoder(d, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logicalxview, 1, d.k, pre_X)
    Z_error = evaluate_classical_decoder(d, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logicalzview, d.k + 1, 2 * d.k)
    return (X_error, Z_error)
end

function evaluate_classical_decoder(d::AbstractSyndromeDecoder, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logical_view_function, guess_start, guess_stop, pre_circuit = nothing)
    H = d.H
    O = d.faults_matrix
    syndrome_circuit = syndrome_circuit_func(H)

    n = d.n
    s = d.s
    k = d.k

    errors = [PauliError(i, init_error) for i in 1:n];

    md = MixedDestabilizer(H)

    full_circuit = []

    logview = logical_view_function(md)
    logcirc, _ = syndrome_circuit_func(logview)

    noisy_syndrome_circuit = add_two_qubit_gate_noise(syndrome_circuit, gate_error);

    for gate in logcirc
        type = typeof(gate)
        if type == sMRZ
            push!(syndrome_circuit, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(syndrome_circuit, type(gate.q1, gate.q2+s))
        end
    end

    ecirc = encoding_circuit_func(syndrome_circuit)
    if isnothing(pre_circuit)
        full_circuit = vcat(pre_circuits, ecirc, errors, noisy_syndrome_circuit)
    else
        full_circuit = vcat(ecirc, errors, noisy_syndrome_circuit)
    end

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, full_circuit)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logical_syndromes = pfmeasurements(frames)[:, s+1: s+k]

    for i in 1:nsamples
        guess = decode(d, syndromes[i])

        # result should be concatinated guess of the X and Z checks
        result = (O * (guess))[guess_start:guess_stop]

        if result == logical_syndromes[i]
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
    """The time taken to create the lookup table + decode the code a specified number of time"""
    time
end

function TableDecoder(Hx, Hz)
    c = CSS(Hx, Hz)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    lookup_table, time, _ = @timed create_lookup_table(H)
    faults_matrix = faults_matrix(H)
    return TableDecoder(H, n, s, k, faults_matrix, lookup_table, time)
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

function decode(d::TableDecoder, syndrome_sample)
    return get(d.lookup_table, syndrome_sample, nothing)
end

function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.numchecks_X]
    row_z = syndrome_sample[d.numchecks_X+1:d.numchecks_X+d.numchecks_Z]

    KguessX, success = syndrome_decode(d.sparse_Cx, d.sparse_CxT, d.row_x, d.max_iters, d.channel_probs, d.b2c_X, d.c2b_X, d.log_probabs, Base.copy(d.err))
    KguessZ, success = syndrome_decode(d.sparse_Cz, d.sparse_CzT, d.row_z, d.max_iters, d.channel_probs, d.b2c_Z, d.c2b_Z, d.log_probabs, Base.copy(d.err))
    guess = vcat(KguessZ, KguessX)
end


## NOT WORKING
function evaluate_classical_decoder(H, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logical_view_func, decoder_func, pre_circuit = nothing)
    decoded = 0

    H_stab = Stabilizer(fill(0x0, size(Hx, 2)), H, zeros(Bool, size(H)))

    O = faults_matrix(H_stab)
    syndrome_circuit = syndrome_circuit_func(H_stab)

    s, n = size(H)
    k = n - s

    errors = [PauliError(i, init_error) for i in 1:n];

    md = MixedDestabilizer(H_stab)

    full_circuit = []

    logview = logical_view_func(md)
    logcirc, _ = syndrome_circuit_func(logview)

    noisy_syndrome_circuit = add_two_qubit_gate_noise(syndrome_circuit, gate_error);

    for gate in logcirc
        type = typeof(gate)
        if type == sMRZ
            push!(circuit, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(circuit, type(gate.q1, gate.q2+s))
        end
    end

    ecirc = encoding_circuit_func(syndrome_circuit)
    if isnothing(pre_circuit)
        full_circuit = vcat(pre_circuits, ecirc, errors, noisy_syndrome_circuit)
    else
        full_circuit = vcat(ecirc, errors, noisy_syndrome_circuit)
    end

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, full_circuit)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logical_syndromes = pfmeasurements(frames)[:, s+1: s+k]

    for i in 1:nsamples
        guess = decode(decoder_obj, syndromes[i]) # TODO: replace 'decoder_obj' with proper object

        result = (O * (guess))

        if result == logical_syndromes[i]
            decoded += 1
        end
    end

    return (nsamples - decoded) / nsamples
end
