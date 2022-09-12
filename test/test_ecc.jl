using QuantumClifford
using QuantumClifford.ECC: AbstractECC, Steane5, Steane7, Shor9, Bitflip3, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate
#using Tests

function test_op(c::AbstractECC)
    @testset "Physical and Logical qubit" begin
        physicalqubit = random_stabilizer(code_k(c))

        gate = rand((:x,:z))
        physicalgate, logicalgate = if gate==:x
            P"X",logx_ops(c)[1]
        elseif gate == :z
            P"Z",logz_ops(c)[1]
        end

        #run 1
        #start physical state
        physicalqubit1 = copy(physicalqubit)

        #perform physical gate
        apply!(physicalqubit1,physicalgate)
        #encode into logical state
        bufferqubits1 = one(Stabilizer,code_s(c))
        logicalqubit1 = physicalqubit1⊗bufferqubits1 
        for gate in encoding_circuit(c)
            apply!(logicalqubit1,gate)
        end
        

        #run 2
        #start same physical state
        physicalqubit2 = copy(physicalqubit)
        #encode logical state
        bufferqubits2 = one(Stabilizer,code_s(c))
        logicalqubit2 = physicalqubit2⊗bufferqubits2 
        for gate in encoding_circuit(c)
            apply!(logicalqubit2,gate)
        end
        #apply logical gate
        apply!(logicalqubit2,logicalgate)

        @test canonicalize!(logicalqubit1) == canonicalize!(logicalqubit2)
       
    end
end


codes = [Steane5(),Steane7(),Shor9(),Bitflip3()] #fix other encoding circuits

for c in codes
    test_op(c)
end
