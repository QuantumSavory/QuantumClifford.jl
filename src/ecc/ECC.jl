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

"""Check if the code is degenerate with respect to single-qubit physical errors."""
function is_degenerate(c::AbstractECC)
    tableau = stab_to_gf2(parity_checks(c))
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

"""The rank of the bimatrix of a code."""
function rank(c::AbstractECC)
    destab_gott = MixedDestabilizer(parity_checks(c), undoperm=false)
    bimat = stab_to_gf2(stabilizerview(destab_gott))
    rank = 0
    for i in 1:code_s(c)
        if bimat[i, i] == 1
            rank +=1
        end
    end
    return rank
end

"""The standardized logical tableau of a code by [PhysRevA.56.76](@cite)"""
function standard_tab_gott(c::AbstractECC)
    n, s, k, r = code_n(c), code_s(c), code_k(c), rank(c)
    # The standard form is 
    # I A1 A2 | B C1 C2
    # 0  0 0  | D  I  E  
    # and we augment the following third line (for logical qubits)
    # 0 E^T I | 0  0  0
    # Then we apply the gates line by line bottom up in accordance with the formalisms here: arXiv:quant-ph/9607030
    standard_tab = stab_to_gf2(stabilizerview(MixedDestabilizer(parity_checks(c), undoperm=false)))
    for i in 1: k
        # can we use canonicalize_gott for the entire mixedDestabilizer? (i.e. the stab + log parts = n rows)
        augment = zeros(Int8, (1, 2*n))
        for j in 1:n
            if j > r && j <= n - k
                augment[j] = standard_tab[j, 2*n-k+i] # the corresponding column of E in E^T 
            elseif j == n-k+i 
                augment[j] = 1
            end
        end
        standard_tab = vcat(standard_tab, augment)
    end
    # Flipping the table so it has the same format as the papercode
    res = zeros(Int8, (n, 2n))
    for i in 1:n
        for j in 1:n
            res[i, j] = standard_tab[n+1-j, n+1-i]
        end
        for j in n+1:2*n
            res[i,j] = standard_tab[2*n+1-j, n+1-i]
        end
    end
    return res
end

function standard_code_tab(c::AbstractECC)
    n, s, k = code_n(c), code_s(c), code_k(c)
    standard_tab = stab_to_gf2(stabilizerview(MixedDestabilizer(parity_checks(c), undoperm=false)))
    md = MixedDestabilizer(parity_checks(c), undoperm=false)
    XZᵗ = stab_to_gf2(stabilizerview(md))
    X₂ = reverse(XZᵗ[1:s,1:n]', dims=(1,2))
    Z₂ = reverse(XZᵗ[1:s,n+1:2n]', dims=(1,2))
    X = falses(n, n)
    Z = falses(n, n)
    X[:, k+1:end] .= X₂
    Z[:, k+1:end] .= Z₂
    X[:, 1:k] .= Xˡ
    return X, Z
    # TODO The permutations need to be reverted at some point, otherwise the generated circuit will have permuted qubits. It is not clear to me whether it is more convenient to revert the permutation here or in naive_encoding_circuit
    # TODO This function does not seem to actually be necessary. It seems like iterating over the `MixedDestabilizer` object is enough. Out of convenience, let's keep if for the moment, but also keep this TODO
end

""" The naive implementation of the encoding circuit by arXiv:quant-ph/9607030 """
function naive_encoding_circuit(c::AbstractECC)
    n, k, s, r = code_n(c), code_k(c), code_s(c), rank(c)
    naive_ec = AbstractOperation[]
    # Applying the hadamard gate to the last r qubits
    for i in n: -1: n-r+1
        push!(naive_ec, sHadamard(i))
    end

    standard_tab = standard_tab_gott(c)
    
    for i in 1 : n
        if standard_tab[i, i] ==1
            for j in 1:n
                if j == i continue end
                gate = (standard_tab[j, i], standard_tab[j, i+n])
                if gate == (1,0)
                    push!(naive_ec, sCNOT(j, i))
                elseif gate == (0,1)
                    push!(naive_ec, sXCZ(j, i))
                elseif gate == (1,1)
                    push!(naive_ec, sXCY(j, i))
                end
            end
        end        
    end
    return naive_ec
end 


include("./bitflipcode.jl")
include("./shorcode.jl")
include("./steanecode.jl")
include("./papercode.jl")

end #module
