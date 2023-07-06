using Test
using QuantumClifford
using QuantumClifford.ECC: AbstractECC, Steane7, Shor9, Bitflip3, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops
include("../src/ecc/ECC.jl")

codes = [
    Bitflip3(),
    Steane7(),
    Shor9(),
]

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

function test_naive_syndrome(c::AbstractECC)
    # create a random logical state
    unencoded_qubits = random_stabilizer(code_k(c))
    bufferqubits = one(Stabilizer,code_s(c))
    logicalqubits = unencoded_qubits⊗bufferqubits
    mctrajectory!(logicalqubits, encoding_circuit(c))
    # measure using `project!`
    s1 = copy(logicalqubits)
    syndrome1 = [project!(s1, check)[3] for check in parity_checks(c)]
    # measure using `naive_syndrome_circuit`
    naive_circuit = naive_syndrome_circuit(c)
    ancillaryqubits = one(Stabilizer,code_s(c))
    s2 = copy(logicalqubits)
    syndrome2 = Register(s2⊗ancillaryqubits, falses(code_s(c)))
    mctrajectory!(syndrome2, naive_circuit)
    @test all(syndrome1 .== 0)
    @test all(bitview(syndrome2) .== 0)
    @test bitview(syndrome2) == syndrome1.÷2

    # TODO test when there is potential for errors / non-commuting operators
end

@testset "naive syndrome circuits - zero syndrome for logical states" begin
    for c in codes, _ in 1:2
        test_naive_syndrome(c)
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


function test_is_degenerate(c::AbstractECC)
    if c == Shor9()
        @test is_degenerate(c) == true
    elseif c == Steane7()
        @test is_degenerate(c) == false
    elseif c== Bitflip3()
        @test is_degenerate(c) == true
    end
end

@testset "is degenerate function - test on popular codes" begin
    for c in codes
        test_is_degenerate(c)
    end
end