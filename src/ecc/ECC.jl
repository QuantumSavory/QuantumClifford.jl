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
    function code_n end
    
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
    function naive_syndrome_circuit(c:: AbstractECC)
        naive_sc = []        
         
        ancilla_bit = 1  
        # ancilla_qubit = code_n(c) + ancilla_bit    
        pc = parity_checks(c)
        for check in pc
            ancilla_qubit = code_n(c) + ancilla_bit
            for qubit in 1: code_n(c)
                if check[qubit] == (1,0)
                    h1 = sHadamard(qubit, ancilla_qubit)
                    c1 = sCNOT(qubit, ancilla_qubit)
                    h2 = sHadamard(qubit, ancilla_qubit)
                    push!(naive_sc, h1)
                    push!(naive_sc, c1)
                    push!(naive_sc, h2)
                else if check[qubit] == (0,1)
                    c1 = sCNOT(qubit, ancilla_qubit)
                    push!(naive_sc, c1)
                end
                
            end
            mz = sMZ(ancilla_qubit, ancilla_bit)
            push!(naive_sc, mz)
            ancilla_bit +=1             
        end

        return naive_sc
    end

    # function naive_syndrome_circuit(c::AbstractECC)
    #     naive_sc = []
    #     dim_encondingc = code_s(c) #previously code_n(c)-1?
    
    #     ancilla_qubit = dim_encondingc +1
    #     tracking1 = dim_encondingc 
    
    #     #iterating through all the steps of the encoding circuit
    #     for qubit in 1:tracking1+1
    #         tracking2 = qubit+1
    #         if qubit < dim_encondingc
    #             push!(naive_sc, sCNOT(qubit,ancilla_qubit)) 
    #             push!(naive_sc, sCNOT(tracking2,ancilla_qubit)) 
    #             push!(naive_sc, sMX(ancilla_qubit))
    #             ancilla_qubit = ancilla_qubit + 1
    #             tracking2 = tracking2 + 1
    #         end
    #     end

    #     convert(Vector{AbstractSymbolicOperator},naive_sc)

    #     return naive_sc
    # end

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
            s = parity_matrix(c) * e 
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
end #module
