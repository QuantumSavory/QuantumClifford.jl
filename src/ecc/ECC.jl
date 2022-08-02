
module ECC

    using QuantumClifford

    abstract type AbstractECC end

    function sCNOT_gatechange(gate::sCNOT) sCNOT(gate.q2,gate.q1) end

    #syndrome extraction circuit 
    function naive_syndrome_circuit(code::AbstractECC) #TODO: add the other 3 types of syndrome circuit
        #AbstractECC type
        encondingc = encoding_circuit(code)
        naive_sc = []
        dim_encondingc = length(encondingc)
        ancilla_qubit = dim_encondingc+1
        tracking1 = (dim_encondingc - 1)
        tracking2 = 2
        #iterating through all the steps of the encoding circuit
        for qubit in 1:tracking1
            append!(naive_sc, sCNOT(1,ancilla_qubit))
            append!(naive_sc, sCNOT(tracking2,ancilla_qubit))
            ancilla_qubit + 1
            tracking2 +1
        end
        return naive_sc
        #fault tolerant (3 types) -Neil, Steane, Shor
    end 

    #=
    function shor_syndrome_circuit(code::AbstractECC) #TODO: add the other 3 types of syndrome circuit
        #p20 of https://young.physics.ucsc.edu/150/error_corr.pdf
        #AbstractECC type
        encondingc = encoding_circuit(code)
        naive_sc = []
        ancilla_qubit = 1
        #iterating through all the steps of the encoding circuit
        dim_encondingc = length(encondingc)
        for gate in parity_checks(code) 
            for qubit in 1:dim_encondingc
                if non identy gate for qubit in gate
                    append!(naive_sc, sCNOT(ancilla_qubit,qubit))
                end
                
            end
            ancilla_qubit + 1
            return naive_sc
        end
        #fault tolerant (3 types) -Neil, Steane, Shor
    end 
    =#

    """The encoding circuit of a given code."""
    function encoding_circuit end 

    """The number of physical qubits in a code."""
    function code_n end

    """The number of stabilizer checks in a code."""
    function code_s(code::AbstractECC)
        #s = sizeof(parity_checks(code)) / code_n(code)
        #s = nrow(parity_checks(code))
        size = size(parity_checks(code))[1]
        s = size / code_n(code)
        return s
    end

    """The number of logical qubits in a code."""
    function code_k(code::AbstractECC)
        k = css_n(code) - code_s(code)
        return k
    end


    """The rate of a code."""
    function rate(code::AbstractECC)
        rate = code_k(code)/code_s(code)
        return rate
    end

    """The distance of a code."""
    function distance end

    """Parity check tableau of a code."""
    function parity_checks end

    """Parity matrix of a code."""
    function parity_matrix end

    """Logical X operations of a code."""
    function logx_ops end # can be computed from the parity checks

    """Logical Z operations of a code."""
    function logz_ops end # can be computed from the parity checks

    """Logical Y operations of a code."""
    function logy_ops end # can be computed from the parity checks

    """Is the code degenerate"""
    function isdegenerate end

    #------IN PROGRESS-------------------------
    # 3 qubit bit flip code 
    #include("./bitflipcode.jl")

    # Shor code
    include("./shorcode.jl")

    # Steane, 7 qubit, 5 qubit
    include("./steanecode.jl")

    # CSS codes
    #include("./csscodes.jl")

    #------TODO--------------------------------

    # Toric codes
    #include("./toriccode.jl") #has to b commented for tsting

    # Surface codes
    #include("./surface_codes.jl")
    #include...
    # repetition codes
    # concatenated codes
    # color codes
    # LDPC expander and hypergraph codes
    # reed muller and reed solomon codes

    # Add others ------------------------

end
