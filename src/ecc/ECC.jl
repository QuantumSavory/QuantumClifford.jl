module Check end
module ECC

    using QuantumClifford
    abstract type AbstractECC end
    #export AbstractECC

    function sCNOT_gatechange(gate::sCNOT) sCNOT(gate.q2,gate.q1) end

    """The encoding circuit of a given code."""
    function encoding_circuit end 

    """Parity check tableau of a code."""
    function parity_checks end

    """The number of physical qubits in a code."""
    function code_n end

    """The number of stabilizer checks in a code."""
    function code_s end

    """The number of logical qubits in a code."""
    function code_k end

    """The rate of a code."""
    function rate end

    """Naive syndrome circuit""" #TODO: add the other 3 types of syndrome circuit
    function naive_syndrome_circuit(c::AbstractECC)
        naive_sc = []
        dim_encondingc = code_n(c) -1
    
        ancilla_qubit = dim_encondingc
        tracking1 = dim_encondingc 
    
        #iterating through all the steps of the encoding circuit
        for qubit in 0:tracking1
            tracking2 = qubit+1
            if qubit < dim_encondingc
                push!(naive_sc, sCNOT(qubit,ancilla_qubit)) 
                push!(naive_sc, sCNOT(tracking2,ancilla_qubit)) 
                ancilla_qubit = ancilla_qubit + 1
                tracking2 = tracking2 + 1
            end
        end

        return naive_sc
    end

    """The distance of a code."""
    function distance end

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
    include("./bitflipcode.jl")

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
