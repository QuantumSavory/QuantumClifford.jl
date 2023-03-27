module ECC

    using QuantumClifford
    abstract type AbstractECC end
    #export AbstractECC

    function sCNOT_gatechange(gate::sCNOT) 
        sCNOT(gate.q2,gate.q1) 
    end

    """The encoding circuit of a given code."""
    function encoding_circuit end 

    """Parity check tableau of a code."""
    function parity_checks end

    """The number of physical qubits in a code."""
    function code_n(c::AbstractECC)
        return length(c[0]) #do we even need to return in this function?
    end
    
    """The number of stabilizer checks in a code."""
    function code_s(c::AbstractECC) 
        s = length(parity_checks(c))
        return s
    end
    
    """The number of logical qubits in a code."""
    function code_k(c::AbstractECC)
        k = code_n(c) - code_s(c)
        return k
    end
    
    """The rate of a code."""
    function rate(c::AbstractECC)
        rate = code_k(c)//code_s(c)
        return rate
    end
    
    """Naive syndrome circuit""" 
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

        convert(Vector{AbstractSymbolicOperator},naive_sc)

        return naive_sc
    end

    """The distance of a code."""
    function distance end
    
    """Parity matrix of a code."""
    function parity_matrix(c::AbstractECC) 
        paritym = stab_to_gf2(parity_checks(c::AbstractECC)) 
        return paritym
    end
    
    """Logical X operations of a code."""
    function logx_ops(c::AbstractECC)
        MixedDest = MixedDestabilizer(parity_checks(c))
        logicalxview(MixedDest)
    end 

    """Logical Z operations of a code."""
    function logz_ops(c::AbstractECC)
        MixedDest = MixedDestabilizer(parity_checks(c))
        logicalzview(MixedDest)
    end 


    #NEW!
    """Is the code degenerate"""
    function isdegenerate(c::AbstractECC, errors)
        syndromes = Set()
        for e in errors
            s = c ⊻ e
            if s in syndromes
                return true
            end
            syndromes.add(s)
        end
        return false
    end

    #------IN PROGRESS-------------------------
    # 3 qubit bit flip code 
    include("./bitflipcode.jl")

    # Shor code
    include("./shorcode.jl")

    # Steane, 7 qubit, 5 qubit
    include("./steanecode.jl")

    # CSS codes
    include("./csscodes.jl")

    #Hamming code generator (for CSS codes)
    include("./hammingcodegenerator.jl")

    # Toric codes
    include("./toriccode.jl") 

    # Surface codes
    include("./surfacecodes.jl")

end #module
