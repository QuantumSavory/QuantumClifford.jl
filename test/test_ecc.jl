using Test
using QuantumClifford
using QuantumClifford.ECC: AbstractECC, Cleve8, Steane7, Shor9, Bitflip3, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate, rank, standard_tab_gott

codes = [
    Bitflip3(),
    Steane7(),
    Shor9(),
    # Cleve8(),
]

##

# These are the old, manually written encoding circuits. They are used to test whether the new algorithmically constructed encoding circuit function works.

function manual_encoding_circuit(c::Bitflip3)
    c1 = sCNOT(1,2)
    c2 = sCNOT(1,3)
    return [c1,c2]
end

function manual_encoding_circuit(c::Shor9)
    c1 = sCNOT(1,4)
    c2 = sCNOT(1,7)

    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    c3 = sCNOT(1,2)
    c4 = sCNOT(4,5)
    c5 = sCNOT(7,8)

    c6 = sCNOT(1,3)
    c7 = sCNOT(4,6)
    c8 = sCNOT(7,9)

    # XXX: The extra sHadamard(1) at the start is due to a popular mismatch in
    # conventions for which logical operator is the X one and which is the Z one
    return [sHadamard(1),c1,c2,h1,h2,h3,c3,c4,c5,c6,c7,c8]
end

function manual_encoding_circuit(c::Steane7)
    sc1 = sCNOT(1,2)
    sc2 = sCNOT(1,3)

    sh1 = sHadamard(5)
    sh2 = sHadamard(6)
    sh3 = sHadamard(7)

    sc3 = sCNOT(7,4)
    sc4 = sCNOT(7,2)
    sc5 = sCNOT(7,1)
    sc6 = sCNOT(6,4)
    sc7 = sCNOT(6,3)
    sc8 = sCNOT(6,1)
    sc9 = sCNOT(5,4)
    sc10 = sCNOT(5,3)
    sc11 = sCNOT(5,2)

    return [sc1,sc2,sh1,sh2,sh3,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10,sc11]
end

@testset "encoding circuits - manual vs algorithmic" begin
    for c in codes
        manual = manual_encoding_circuit(c)
        algorithmic, perm = encoding_circuit(c)
        # init = random_stabilizer(code_k(c))⊗one(Stabilizer,code_s(c))
        init = one(Stabilizer, code_k(c))⊗one(Stabilizer,code_s(c))
        fin_m = mctrajectory!(copy(init), manual)[1] |> canonicalize!
        fin_a = mctrajectory!(copy(init), algorithmic)[1] |> canonicalize!
        @test fin_m == fin_a
    end
end

##

function test_phys_log_op_encoding(c::AbstractECC)
    # encode k physical qubits into n physical qubits (turning them into k logical qubits)
    # apply operations in an encoded or unencoded fashion
    # compare results
    physicalqubit = random_stabilizer(code_k(c))

    gate = rand((:x,:z))
    target = rand(1:code_k(c))

    physicalgate, logicalgate =
    if gate==:x
        P"X",logx_ops(c)[target]
    elseif gate == :z
        P"Z",logz_ops(c)[target]
    end

    #run 1
    #start physical state
    physicalqubit1 = copy(physicalqubit)
    #apply physical gate
    apply!(physicalqubit1,physicalgate)
    #encode into logical state
    bufferqubits1 = one(Stabilizer,code_s(c))
    logicalqubit1 = physicalqubit1⊗bufferqubits1 # pad up the k physical qubits into a state of n physical qubits
    mctrajectory!(logicalqubit1, encoding_circuit(c))

    #run 2
    #start same physical state
    physicalqubit2 = copy(physicalqubit)
    #encode logical state
    bufferqubits2 = one(Stabilizer,code_s(c))
    logicalqubit2 = physicalqubit2⊗bufferqubits2
    mctrajectory!(logicalqubit2, encoding_circuit(c))
    #apply logical gate
    apply!(logicalqubit2,logicalgate)

    @test canonicalize!(logicalqubit1) == canonicalize!(logicalqubit2)
end

@testset "physical vs optical operators - check of encoding circuit" begin
    for c in codes, _ in 1:2
        test_phys_log_op_encoding(c)
    end
end

##

function test_naive_syndrome(c::AbstractECC, e::Bool=false)
    # create a random logical state
    unencoded_qubits = random_stabilizer(code_k(c))
    bufferqubits = one(Stabilizer,code_s(c))
    logicalqubits = unencoded_qubits⊗bufferqubits
    mctrajectory!(logicalqubits, encoding_circuit(c))
    if e
        #add some noise to logicalqubits
        apply!(logicalqubits, P"X", rand(1:code_n(c)))
        apply!(logicalqubits, P"Z", rand(1:code_n(c)))
    end
    # measure using `project!`
    s1 = copy(logicalqubits)
    syndrome1 = [project!(s1, check)[3] for check in parity_checks(c)]
    # measure using `naive_syndrome_circuit`
    naive_circuit = naive_syndrome_circuit(c)
    ancillaryqubits = one(Stabilizer,code_s(c))
    s2 = copy(logicalqubits)
    syndrome2 = Register(s2⊗ancillaryqubits, falses(code_s(c)))
    mctrajectory!(syndrome2, naive_circuit)
    if !e
        @test all(syndrome1 .== 0)
        @test all(bitview(syndrome2) .== 0)
    end
    @test bitview(syndrome2) == syndrome1.÷2

    # TODO test when there is potential for errors / non-commuting operators
end

@testset "naive syndrome circuits - zero syndrome for logical states" begin
    for c in codes, _ in 1:2
        test_naive_syndrome(c)
        test_naive_syndrome(c, true)
    end
end

##

function test_with_pframes(code)
    ecirc = encoding_circuit(code)
    scirc = naive_syndrome_circuit(code)
    nframes = 10
    dataqubits = code_n(code)
    ancqubits = code_s(code)
    regbits = ancqubits
    frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
    circuit = [ecirc..., scirc...]
    pftrajectories(frames, circuit)
    @test sum(pfmeasurements(frames)) == 0
end

@testset "naive syndrome circuits - zero syndrome for logical states" begin
    for c in codes, _ in 1:2
        test_with_pframes(c)
    end
end

##

@testset "is degenerate function - test on popular codes" begin
    @test isdegenerate(Shor9()) == true
    @test isdegenerate(Steane7()) == false
    @test isdegenerate(Steane7(), 2) == true
    @test isdegenerate(Bitflip3()) == true
end
