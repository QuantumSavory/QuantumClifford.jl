"""Generate a lookup table for decoding single qubit errors."""
function create_lookup_table(code::Stabilizer)
    lookup_table = Dict()
    constraints, qubits = size(code)
    # In the case of no errors
    lookup_table[ zeros(UInt8, constraints) ] = zero(PauliOperator, qubits)
    # In the case of single bit errors
    for bit_to_be_flipped in 1:qubits
        for error_type in [single_x, single_y, single_z]
            # Generate e⃗
            error = error_type(qubits, bit_to_be_flipped)
            # Calculate s⃗
            # (check which stabilizer rows do not commute with the Pauli error)
            syndrome = comm(error, code)
            # Store s⃗ → e⃗
            lookup_table[syndrome] = error
        end
    end
    lookup_table
end

"""Given a syndrome circuit, returns the fault tolerant encoding circuit. Basically just a copy of the syndrome circuit that throws away the measurement results."""
function fault_tolerant_encoding(scirc)
    ecirc = Vector{QuantumClifford.AbstractOperation}()
    for gate in scirc
        if isa(gate,sMRZ)
            push!(ecirc, sMZ(gate.qubit))
        elseif isa(gate,sMRX)
            push!(ecirc, sMX(gate.qubit))
        elseif isa(gate,sMRY)
            push!(ecirc, sMY(gate.qubit))
        else
            push!(ecirc, gate)
        end
    end
    return ecirc
end

""" Most simple simulation function
* Uses a lookup table decoder over single qubit errors ['create_lookup_table'](@ref)
* Uses ['naive_syndrome_circuit'](@ref)
Returns a logical X and Z error for the provided p_error - physical error rate per qubit
"""
function naive_error_correction_pipeline(checks::Stabilizer, p_error; nframes=10_000, ecirc=nothing, encoding_locs=nothing, scirc=nothing)
    if isnothing(scirc)
        scirc , _= naive_syndrome_circuit(checks)
    end
    lookup_table = create_lookup_table(checks)
    O = faults_matrix(checks)
    circuit_Z = Base.copy(scirc)
    circuit_X = Base.copy(scirc)

    s, n = size(checks)
    k = n-s

    if isnothing(encoding_locs)
        pre_X = [sHadamard(i) for i in n-k+1:n]
    else
        pre_X = [sHadamard(i) for i in encoding_locs]
    end
    
    md = MixedDestabilizer(checks)
    logview_Z = [ logicalzview(md);]
    logcirc_Z, _ = naive_syndrome_circuit(logview_Z) # numLogBits shoudl equal k

    logview_X = [ logicalxview(md);]
    logcirc_X, _ = naive_syndrome_circuit(logview_X)

    # Z logic circuit
    for gate in logcirc_Z
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_Z, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(circuit_Z, type(gate.q1, gate.q2+s))
        end
    end

   # X logic circuit
    for gate in logcirc_X
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_X, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(circuit_X, type(gate.q1, gate.q2+s))
        end
    end

    # Z simulation
    errors = [PauliError(i,p_error) for i in 1:n]
    if isnothing(ecirc)
        ecirc_z = fault_tolerant_encoding(circuit_Z) # double syndrome encoding
    else
        ecirc_z = ecirc
    end
    fullcircuit_Z = vcat(ecirc_z, errors, circuit_Z)

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, fullcircuit_Z)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logicalSyndromes = pfmeasurements(frames)[:, s+1: s+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        guess = get(lookup_table,row,nothing)
        if isnothing(guess)
            continue
        else
            result_Z = (O * stab_to_gf2(guess))[k+1:2k]
            if result_Z == logicalSyndromes[i,:]
                decoded += 1
            end
        end
    end
    z_error = 1 - decoded / nframes

    # X simulation
    if isnothing(ecirc)
        ecirc_x = fault_tolerant_encoding(circuit_X)  # double syndrome encoding
    else
        ecirc_x = ecirc
    end
    fullcircuit_X = vcat(pre_X, ecirc_x, errors, circuit_X)
    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, fullcircuit_X)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logicalSyndromes = pfmeasurements(frames)[:, s+1: s+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        guess = get(lookup_table,row,nothing)
        if isnothing(guess)
            continue
        else
            result_X = (O * stab_to_gf2(guess))[1:k]
            if result_X == logicalSyndromes[i, :]
                decoded += 1
            end
        end
    end
    x_error = 1 - decoded / nframes

    return x_error, z_error
end

""" Similiar to ['naive_error_correction_pipeline'](@ref) but now uses shor style fault tolerant syndrome measurement
* Uses a lookup table decoder over single qubit errors ['create_lookup_table'](@ref)
* Uses ['shor_syndrome_circuit'](@ref)
Returns a logical X and Z error for the provided p_error - physical error rate per qubit
"""
function shor_error_correction_pipeline(checks::Stabilizer, p_init; nframes=10_000, cat=nothing, scirc=nothing, ecirc=nothing, encoding_locs=nothing)
    if isnothing(scirc) || isnothing(cat)
        cat, scirc, _ = shor_syndrome_circuit(checks)
    end

    lookup_table = create_lookup_table(checks)
    O = faults_matrix(checks)
    circuit_Z = Base.copy(scirc)
    circuit_X = Base.copy(scirc)

    s, n = size(checks)
    k = n-s

    if isnothing(encoding_locs)
        pre_X = [sHadamard(i) for i in n-k+1:n]
    else
        pre_X = [sHadamard(i) for i in encoding_locs]
    end

    anc_qubits = 0
    for pauli in checks
        anc_qubits += mapreduce(count_ones,+, xview(pauli) .| zview(pauli))
    end
    regbits = anc_qubits + s

    md = MixedDestabilizer(checks)
    logview_Z = logicalzview(md)
    logcirc_Z, _ = naive_syndrome_circuit(logview_Z)

    logview_X = logicalxview(md)
    logcirc_X, _ = naive_syndrome_circuit(logview_X)

    # Z logic circuit
    for gate in logcirc_Z
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_Z, sMRZ(gate.qubit+anc_qubits, gate.bit+regbits))
        else
            push!(circuit_Z, type(gate.q1, gate.q2+anc_qubits))
        end
    end

   # X logic circuit
    for gate in logcirc_X
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_X, sMRZ(gate.qubit+anc_qubits, gate.bit+regbits))
        else
            push!(circuit_X, type(gate.q1, gate.q2+anc_qubits))
        end
    end

    # Z simulation
    errors = [PauliError(i,p_init) for i in 1:n]
    if isnothing(ecirc)
        ecirc_z = fault_tolerant_encoding(vcat(cat,circuit_Z)) # double syndrome encoding
        fullcircuit_Z = vcat(ecirc_z, errors, circuit_Z)
    else
        ecirc_z = ecirc
        fullcircuit_Z = vcat(ecirc_z, errors, cat, circuit_Z)
    end
    
    frames = PauliFrame(nframes, n+anc_qubits+k, regbits+k)
    pftrajectories(frames, fullcircuit_Z)
    syndromes = pfmeasurements(frames)[:, anc_qubits+1:regbits]
    logicalSyndromes = pfmeasurements(frames)[:, regbits+1:regbits+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        guess = get(lookup_table,row,nothing)
        if isnothing(guess)
            continue
        else
            result_Z = (O * stab_to_gf2(guess))[k+1:2k]
            if result_Z == logicalSyndromes[i,:]
                decoded += 1
            end
        end
    end
    z_error = 1 - decoded / nframes

    # X simulation
    if isnothing(ecirc)
        ecirc_x = fault_tolerant_encoding(vcat(cat,circuit_X)) # double syndrome encoding
        fullcircuit_X = vcat(pre_X, ecirc_x, errors, circuit_X)
    else
        ecirc_x = ecirc
        fullcircuit_X = vcat(pre_X, ecirc_x, errors, cat, circuit_X)
    end
    
    frames = PauliFrame(nframes, n+anc_qubits+k, regbits+k)
    pftrajectories(frames, fullcircuit_X)
    syndromes = pfmeasurements(frames)[:, anc_qubits+1:regbits]
    logicalSyndromes = pfmeasurements(frames)[:, regbits+1:regbits+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        guess = get(lookup_table,row,nothing)
        if isnothing(guess)
            continue
        else
            result_X = (O * stab_to_gf2(guess))[1:k]
            if result_X == logicalSyndromes[i, :]
                decoded += 1
            end
        end
    end
    x_error = 1 - decoded / nframes

    return x_error, z_error
end

struct CSS; tableau; cx; cz; end

"""Naive syndrome measurement on a CSS ECC with a Cx and Cz matrix
- Only wroks with fault tolerant encoding. 
"""
function CSS_naive_error_correction_pipeline(code::CSS, p_init; nframes=1_000, scirc=nothing, max_iters = 25)
    if isnothing(scirc)
        scirc , _= naive_syndrome_circuit(code.tableau)
    end

    O = faults_matrix(code.tableau)
    circuit_Z = Base.copy(scirc)
    circuit_X = Base.copy(scirc)

    @assert size(code.cx, 2) == size(code.cz, 2) == nqubits(code.tableau)
    @assert size(code.cx, 1) + size(code.cz, 1) == length(code.tableau)

    s, n = size(code.tableau)
    _, _, r = canonicalize!(Base.copy(code.tableau), ranks=true)
    k = n - r

    # Krishna decoder
    log_probabs = zeros(n)
    channel_probs = fill(p_init, n)

    numchecks_X = size(code.cx)[1]
    b2c_X = zeros(numchecks_X, n)
    c2b_X = zeros(numchecks_X, n)

    numchecks_Z = size(code.cz)[1]
    b2c_Z = zeros(numchecks_Z, n)
    c2b_Z = zeros(numchecks_Z, n)
    err = zeros(n)

    pre_X = [sHadamard(i) for i in n-k+1:n]

    md = MixedDestabilizer(code.tableau)
    logview_Z = logicalzview(md)
    logcirc_Z, numLogBits_Z, _ = naive_syndrome_circuit(logview_Z)
    @assert numLogBits_Z == k

    logview_X = logicalxview(md)
    logcirc_X, numLogBits_X, _ = naive_syndrome_circuit(logview_X)
    @assert numLogBits_X == k

    # Z logic circuit
    for gate in logcirc_Z
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_Z, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(circuit_Z, type(gate.q1, gate.q2+s))
        end
    end

   # X logic circuit
    for gate in logcirc_X
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_X, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(circuit_X, type(gate.q1, gate.q2+s))
        end
    end

    errors = [PauliError(i,p_init) for i in 1:n]

    # Z simulation
    ecirc = fault_tolerant_encoding(circuit_Z)
    fullcircuit_Z = vcat(ecirc, errors, circuit_Z)

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, fullcircuit_Z)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logicalSyndromes = pfmeasurements(frames)[:, s+1: s+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        row_x = row[1:numchecks_X]
        row_z = row[numchecks_X+1:numchecks_X+numchecks_Z]

        KguessX, success = syndrome_decode(sparse(code.cx), sparse(code.cx'), row_x, max_iters, channel_probs, b2c_X, c2b_X, log_probabs, Base.copy(err))
        KguessZ, success = syndrome_decode(sparse(code.cz), sparse(code.cz'), row_z, max_iters, channel_probs, b2c_Z, c2b_Z, log_probabs, Base.copy(err))
        guess = vcat(KguessZ, KguessX)
        
        result_Z = (O * (guess))[k+1:2k]
        if result_Z == logicalSyndromes[i,:]
            decoded += 1
        end
    end
    z_error = 1 - decoded / nframes

    # X simulation
    ecirc = fault_tolerant_encoding(circuit_X)
    fullcircuit_X = vcat(pre_X, ecirc, errors, circuit_X)

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, fullcircuit_X)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logicalSyndromes = pfmeasurements(frames)[:, s+1: s+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        row_x = row[1:numchecks_X]
        row_z = row[numchecks_X+1:numchecks_X+numchecks_Z]

        KguessX, success = syndrome_decode(sparse(code.cx), sparse(code.cx'), row_x, max_iters, channel_probs, b2c_X, c2b_X, log_probabs, Base.copy(err))
        KguessZ, success = syndrome_decode(sparse(code.cz), sparse(code.cz'), row_z, max_iters, channel_probs, b2c_Z, c2b_Z, log_probabs, Base.copy(err))
        guess = vcat(KguessZ, KguessX)
        
        result_X = (O * (guess))[1:k]
        if result_X == logicalSyndromes[i, :]
            decoded += 1
        end
    end
    x_error = 1 - decoded / nframes

    return x_error, z_error
end

"""Shor syndrome measurement on a CSS ECC with a Cx and Cz matrix
- Only wroks with fault tolerant encoding. 
"""
function CSS_shor_error_correction_pipeline(code::CSS, p_init;  nframes=10_000, cat=nothing, scirc=nothing, max_iters = 25)
    if isnothing(scirc) || isnothing(cat)
        cat, scirc, _ = shor_syndrome_circuit(code.tableau)
    end

    O = faults_matrix(code.tableau)
    circuit_Z = Base.copy(scirc)
    circuit_X = Base.copy(scirc)

    @assert size(code.cx, 2) == size(code.cz, 2) == nqubits(code.tableau)
    @assert size(code.cx, 1) + size(code.cz, 1) == length(code.tableau)

    s, n = size(code.tableau)
    _, _, r = canonicalize!(Base.copy(code.tableau), ranks=true)
    k = n - r

    # Krishna decoder
    log_probabs = zeros(n)
    channel_probs = fill(p_init, n)

    numchecks_X = size(code.cx)[1]
    b2c_X = zeros(numchecks_X, n)
    c2b_X = zeros(numchecks_X, n)

    numchecks_Z = size(code.cz)[1]
    b2c_Z = zeros(numchecks_Z, n)
    c2b_Z = zeros(numchecks_Z, n)
    err = zeros(n)

    pre_X = [sHadamard(i) for i in n-k+1:n]

    anc_qubits = 0
    for pauli in code.tableau
        anc_qubits += mapreduce(count_ones,+, xview(pauli) .| zview(pauli))
    end

    regbits = anc_qubits + s

    md = MixedDestabilizer(code.tableau)
    logview_Z = logicalzview(md)
    logcirc_Z, _ = naive_syndrome_circuit(logview_Z)

    logview_X = logicalxview(md)
    logcirc_X, _ = naive_syndrome_circuit(logview_X)

    # Z logic circuit
    for gate in logcirc_Z
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_Z, sMRZ(gate.qubit+anc_qubits, gate.bit+regbits))
        else
            push!(circuit_Z, type(gate.q1, gate.q2+anc_qubits))
        end
    end

   # X logic circuit
    for gate in logcirc_X
        type = typeof(gate)
        if type == sMRZ
            push!(circuit_X, sMRZ(gate.qubit+anc_qubits, gate.bit+regbits))
        else
            push!(circuit_X, type(gate.q1, gate.q2+anc_qubits))
        end
    end

    errors = [PauliError(i,p_init) for i in 1:n]

    # Z simulation
    ecirc = fault_tolerant_encoding(vcat(cat, circuit_Z))# notice that the ecirc now contains the cat state
    fullcircuit_Z = vcat(ecirc, errors, circuit_Z)

    frames = PauliFrame(nframes, n+anc_qubits+k, regbits+k)
    pftrajectories(frames, fullcircuit_Z)
    syndromes = pfmeasurements(frames)[:, anc_qubits+1:regbits]
    logicalSyndromes = pfmeasurements(frames)[:, regbits+1:regbits+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        row_x = row[1:numchecks_X]
        row_z = row[numchecks_X+1:numchecks_X+numchecks_Z]

        KguessX, success = syndrome_decode(sparse(code.cx), sparse(code.cx'), row_x, max_iters, channel_probs, b2c_X, c2b_X, log_probabs, Base.copy(err))
        KguessZ, success = syndrome_decode(sparse(code.cz), sparse(code.cz'), row_z, max_iters, channel_probs, b2c_Z, c2b_Z, log_probabs, Base.copy(err))
        guess = vcat(KguessZ, KguessX)
        
        result_Z = (O * (guess))[k+1:2k]
        if result_Z == logicalSyndromes[i,:]
            decoded += 1
        end
    end
    z_error = 1 - decoded / nframes

    # X simulation
    ecirc = fault_tolerant_encoding(vcat(cat, circuit_X))
    fullcircuit_X = vcat(pre_X, ecirc, errors, circuit_X) # notice that the ecirc now contains the cat state

    frames = PauliFrame(nframes, n+anc_qubits+k, regbits+k)
    pftrajectories(frames, fullcircuit_X)
    syndromes = pfmeasurements(frames)[:, anc_qubits+1:regbits]
    logicalSyndromes = pfmeasurements(frames)[:, regbits+1:regbits+k]

    decoded = 0
    for i in 1:nframes
        row = syndromes[i,:]
        row_x = row[1:numchecks_X]
        row_z = row[numchecks_X+1:numchecks_X+numchecks_Z]

        KguessX, success = syndrome_decode(sparse(code.cx), sparse(code.cx'), row_x, max_iters, channel_probs, b2c_X, c2b_X, log_probabs, Base.copy(err))
        KguessZ, success = syndrome_decode(sparse(code.cz), sparse(code.cz'), row_z, max_iters, channel_probs, b2c_Z, c2b_Z, log_probabs, Base.copy(err))
        guess = vcat(KguessZ, KguessX)
        
        result_X = (O * (guess))[1:k]
        if result_X == logicalSyndromes[i, :]
            decoded += 1
        end
    end
    x_error = 1 - decoded / nframes

    return x_error, z_error
end