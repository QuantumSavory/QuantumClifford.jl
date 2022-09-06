using QuantumClifford
using QuantumClifford.ECC: AbstractECC, Steane5, Steane7, Shor9, Bitflip3, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate
#using Tests

function test_op(c::AbstractECC)
    @testset "Physical and Logical qubit" begin
        #physicalqubit
        encoding_circuit_physical = encoding_circuit(c)
        physicalqubit = S"X"
        apply!(physicalqubit,P"X")

        #logicalqubit
        encoding_circuit_logical = encoding_circuit(c)

        if c == Steane5()
            ancillary_qubit_count = 3
        elseif c == Steane7()
            ancillary_qubit_count = 4
        elseif c == Shor9()
            ancillary_qubit_count = 8
        elseif c == Bitflip3()
            ancillary_qubit_count = 2
        end

        

        bufferqubits = one(Stabilizer,ancillary_qubit_count)
        logicalqubit = physicalqubit⊗bufferqubits 
        for gate in encoding_circuit_logical
            apply!(logicalqubit,gate)
        end
        #=
        for gate in logx_ops(c)
            apply!(logicalqubit,gate) #logical gate
        end
        =#
        canonicalize!(logicalqubit)

        for gate in encoding_circuit(c)
            @test encoding_circuit_physical == encoding_circuit_logical
        end
    end
end


codes = [Steane5(),Steane7(),Shor9(),Bitflip3()] #requires code generators

for c in codes
    test_op(c)
end
