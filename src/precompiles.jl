function _precompile_()
    rng = Random.default_rng()
    saved_rng = copy(rng)
    try
        Random.seed!(rng, 0x5143)
        ds = random_destabilizer(3)
        s = random_stabilizer(3)
        canonicalize!(s)
        canonicalize_rref!(s)
        canonicalize_gott!(s)
        canonicalize!(s; phases=false)
        canonicalize_gott!(s; phases=false)
        canonicalize_rref!(s; phases=false)
        c = random_clifford(3)
        apply!(s, c)
        apply!(s, [1, 2], tCNOT)
        apply!(s, sCNOT(1, 2))
        project!(s, P"XXX")
        s = S"XX
              ZZ"
        md = MixedDestabilizer(s)
        p = P"XX"
        c = C"X_
              _X
              Z_
              _Z"
        p * md
        c * md
        apply!(md, p; phases=false)
        apply!(md, c; phases=false)
        random_clifford(3)
        random_pauli(3)
        for op in [sHadamard, sId1, sInvPhase, sPhase, sX, sY, sZ]
            op(2) * s
        end
        for op in [sCNOT, sSWAP, sCPHASE]
            op(2, 1) * s
        end
        project!(s, p)
        project!(md, p)
    finally
        copy!(rng, saved_rng)
    end
end

import Random
using PrecompileTools

@setup_workload begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        _precompile_()
    end
end

@setup_workload let
    @compile_workload begin
        @assert P"X" * P"Z" == P"-iY"
        @assert P"X" ⊗ P"Z" == P"XZ"
        bell_state = S"-XX
                       +ZZ"
        expected = S"-X_
                     +_Z"
        @assert tCNOT * bell_state == expected
        @assert apply!(copy(bell_state), tCNOT) == expected

        state = S"-XXX
                  +ZZI
                  -IZZ"
        canonicalize!(state)
        mixed = MixedDestabilizer(state)
        projected, anticommuting_index, outcome = project!(mixed, P"ZII")
        @assert anticommuting_index != 0
        @assert isnothing(outcome)
        @assert nqubits(projected) == 3
        _, repeated_index, repeated_outcome = project!(copy(projected), P"ZII")
        @assert repeated_index == 0
        @assert repeated_outcome == 0x00
    end
end

@setup_workload let
    rng = Random.default_rng()
    saved_rng = copy(rng)
    try
        Random.seed!(rng, 0x5143)
        @compile_workload begin
            code = ECC.Steane7()
            encoding_circuit = ECC.naive_encoding_circuit(code)
            syndrome_circuit, ancillaries, syndrome_bits = ECC.naive_syndrome_circuit(code)
            circuit = [encoding_circuit...; syndrome_circuit...]
            frames = pftrajectories(circuit; trajectories=4, threads=false)
            syndromes = measurements(frames)
            @assert ancillaries == ECC.code_s(code) == 6
            @assert syndrome_bits == 1:6
            @assert size(syndromes) == (4, 6)
            @assert all(iszero, syndromes)
        end
    finally
        copy!(rng, saved_rng)
    end
end
