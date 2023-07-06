module ECC

using QuantumClifford
using QuantumClifford: AbstractOperation

abstract type AbstractECC end

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

"""Number of physical qubits for a given parity check tableau"""
function code_n(parity_check_tableau)
    return size(parity_check_tableau)[2]
end

"""Wrapper function for codes of type AbstractECC"""
function naive_syndrome_circuit(code_type::AbstractECC)
    naive_syndrome_circuit(parity_checks(code_type))
end

"""Naive syndrome circuit"""
function naive_syndrome_circuit(parity_check_tableau)
    naive_sc = AbstractOperation[]

    ancilla_bit = 1
    # ancilla_qubit = code_n(c) + ancilla_bit
    for check in parity_check_tableau
        ancilla_qubit = code_n(parity_check_tableau) + ancilla_bit
        for qubit in 1: code_n(parity_check_tableau)
            if check[qubit] == (1,0) # TODO simplify this branch using sXCX and similar gates
                push!(naive_sc, sXCX(qubit, ancilla_qubit))
            elseif check[qubit] == (0,1)
                push!(naive_sc, sCNOT(qubit, ancilla_qubit))
            elseif check[qubit] == (1,1)
                push!(naive_sc, sYCX(qubit, ancilla_qubit))
            end
        end
        mz = sMZ(ancilla_qubit, ancilla_bit)
        push!(naive_sc, mz)
        ancilla_bit +=1
    end

    return naive_sc
end

"""The distance of a code."""
function distance end

"""Parity matrix of a code, given as a stabilizer tableau."""
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

""" The bimatrix representation of a code"""
function tableau_representation(c::AbstractECC)
    n, s = code_n(c), code_s(c)
    tableau = zeros(s, 2 * n)
    parity_check_tableau = parity_checks(c)
    for check in 1:s
        for qubit in 1: n
            tableau[check, qubit], tableau[check, qubit + n] = parity_check_tableau[check, qubit]
        end
    end

    return tableau
end

""" Check if the code is degenerate or not """
function is_degenerate(c::AbstractECC)
    tableau = tableau_representation(c)
    n = code_n(c)
    dictionary = Set()
    for column in 1:2*n 
        temp = tableau[:, column]
        if temp in dictionary
              return true
        else
              push!(dictionary, temp)
        end
    end
    return false
end

""" Canonicalize the logicals x operators of a code by @gottesman1997 arxiv:quant-ph/9705052 """ # TODO: Implement the same code for the Z operators
function canonicalize_logicals(c::AbstractECC)
    n, s, k = code_n(c), code_s(c), code_k(c)
    logx = logx_ops(c)
    tabx = tableau_representation(logx)
    tab = tableau_representation(canonicalize_gott!(parity_checks(c)))
    # Finding the rank r of the logical  X
    r =0
    for i in 1: n 
        pivot =findfirst(tab[:, i])
        if pivot !== nothing
            r += 1
        end
    end
    # standardize u1 and v2 for each element of logX (u1u2u3|v1v2v3)    
    for i in 1:k
        op = tabx[i, :]
        # standardize the first n-k qubits (the u1 and v2 component)
        for j in 1:n-k
            if (j <= r && op[j] == 1) || (j >= r+1 && op[j+n] ==1)
                tabx[i] += tab[j] # TODO: fix this xor plus                 
            end
        end
    end
    # setting u3 = I and v3 = 0
    for i in 1:k
        op = tabx[i, :]
        for j in n-k+1:n
            if j - (n-k) == i
                op[j]=1
            else
                op[j]=0
            end
            op[j+n] = 0            
        end
    end

    return tabx
end

""" The naive implementation of the encoding circuit by arXiv:quant-ph/9607030 """
function naive_encoding_circuit(c::AbstractECC)
    n, k, s = code_n(c), code_k(c), code_s(c)
    naive_ec = AbstractOperation[]
    # Applying the hadamard gate to the last r qubits
    for i in n: -1: n-r+1
        push!(naive_ec, sHadamard(i))
    end
    # The standard form is 
    # I A1 A2 | B C1 C2
    # 0  0 0  | D  I  E  
    # and we augment the following third line (for logical qubits)
    # 0 E^T I | 0  0  0
    # Then we apply the gates line by line bottom up in accordance with the formalisms here: arXiv:quant-ph/9607030
    standard_tab = canonicalize_gott!(parity_checks(c))
    for i in 1: k
        # can we use canonicalize_gott for the entire mixedDestabilizer? (i.e. the stab + log parts = n rows)
        augment = zeros(2*n)
        for j in 1:n
            if j > r && j <= n - k
                augment[j] = standard_tab[r+j, 2*n-k+i] # the corresponding column of E in E^T
            elseif j == n-k+i 
                augment[j] = 1
            end
        end
        push!(standard_tab, augment)
    end

    for i in n: -1: 1 # implement the decoder from the augmented bimatrix bottom up and from right to left
        if standard_tab[i, i] ==1
            for j in n: -1: 1
                if j == i continue end
                gate = (standard_tab[i,j], standard_tab[i,j+n])
                if gate == (1,0)
                    push!(naive_ec, sXCX(j, i))
                elseif gate == (0,1)
                    push!(naive_ec, sCNOT(j, i))
                elseif gate == (1,1)
                    push!(naive_ec, sYCX(j, i))
                end
            end
        end        
    end
    return naive_ec
end 


include("./bitflipcode.jl")
include("./shorcode.jl")
include("./steanecode.jl")

end #module
