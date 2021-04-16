"""
A module for sampling random n-qubit Clifford gates.
Implements the algorithm in https://arxiv.org/abs/2003.09412
"""

#= from Bravyi and Maslov Algorithm 1
sample (h, S) from the distribution P_n(h, S) =#
function quantum_mallows(n)
    arr = [1:n;]
    hadamard = falses(n)
    perm = zeros(Int64, n)
    for idx in [1:n;]
        m = length(arr)
        # sample h_i from given prob distribution
        r = rand()
        weight = Int64(2 * m - ceil(log2(r*(BigInt(4)^m-1) + 1)))
        hadamard[idx] = (weight < m)
        k = weight < m ? weight : 2*m - weight - 1
        perm[idx] = popat!(arr, k + 1) # beware of indexing in julia
    end
    return hadamard, perm
end

#= replacement of QuantumClifford's random_clifford function
from Algorithm 2 of Bravyi and Maslov, with code idioms following
the Python implementation in Qiskit =#
function rand_clifford(n)

    @assert n < 200 # otherwise matrix operations could fail

    hadamard, perm = quantum_mallows(n)
    had_idxs = findall(i -> hadamard[i], [1:n;])
    
    # delta, delta', gamma, gamma' appear in the canonical form
    # of a Clifford operator (Eq. 3/Theorem 1)
    # delta is unit lower triangular, gamma is symmetric
    delta = Array{Int8}(LinearAlgebra.I, n, n)
    delta_p = Array{Int8}(LinearAlgebra.I, n, n)
    gamma = zeros(Int8, n, n)
    gamma_p = Array{Int8}(LinearAlgebra.Diagonal(rand(0:1, n)))
    
    # Gamma_ii is zero if h[i] = 0
    for idx in had_idxs
        gamma[idx, idx] = rand(0:1)
    end

    # gamma' and delta' are unconstrained on the lower triangular
    fill_tril(gamma_p, n, symmetric = true)
    fill_tril(delta_p, n)

    # off diagonal: Gamma, Delta must obey conditions C1-C5
    for row in [1:n;], col in [1:row-1;]
        if hadamard[row] && hadamard[col]
            b = rand(0:1)
            gamma[row, col] = b
            gamma[col, row] = b
            # otherwise delta[row,col] must be zero by C4
            if perm[row] > perm[col]
                 delta[row, col] = rand(0:1)
            end

        elseif hadamard[row] && (!hadamard[col]) && perm[row] < perm[col]
            # C5 imposes delta[row, col] = 0 for h[row]=1, h[col]=0
            # if perm[row] > perm[col] then C2 imposes gamma[row,col] = 0
            b = rand(0:1)
            gamma[row, col] = b
            gamma[col, row] = b
            
        elseif (!hadamard[row]) && hadamard[col]
            delta[row, col] = rand(0:1)

            # not sure what condition imposes this
            if perm[row] > perm[col]
                 b = rand(0:1)
                 gamma[row, col] = b
                 gamma[col, row] = b
            end

        elseif (!hadamard[row]) && (!hadamard[col]) && perm[row] < perm[col]
            # C1 imposes gamma[row, col] = 0 for h[row]=h[col] = 0
            # if perm[row] > perm[col] then C3 imposes delta[row,col] = 0
            delta[row, col] = rand(0:1)
        end
    end

    # now construct the tableau representation for F(I, Gamma, Delta)
    prod = gamma * delta
    prod_p = gamma_p * delta_p
    inv_delta = Array(inv(transpose(delta)))
    inv_delta_p = Array(inv(transpose(delta_p)))
    
    # block matrix form
    F1 = Array{Int8}(mod.([delta zeros(Int8, n, n); prod inv_delta], 2))
    F2 = Array{Int8}(mod.([delta_p zeros(Int8, n, n); prod_p inv_delta_p],2))

    # apply qubit permutation S to F2
    perm_inds = vcat(perm, perm .+ n)
    U = F2[perm_inds,:]
    
    # apply layer of hadamards
    lhs_inds = vcat(had_idxs, had_idxs .+ n)
    rhs_inds = vcat(had_idxs .+ n, had_idxs)
    U[lhs_inds, :] = U[rhs_inds, :]
 
    # apply F1
    xzs = Array{Bool}(mod.(F1 * U,2))
 
    # random Pauli matrix just amounts to phases on the stabilizer tableau
    phases = Array{UInt8}(rand([0x0,0x2], 2 * n))
    return CliffordOperator(Stabilizer(phases, xzs))
end

#= simplified version of Algorithm 2 of Bravyi and Maslov
(closely follows the Python code in Qiskit) =#
function rand_clifford_qiskit(n)

    @assert n < 200 

    hadamard, perm = quantum_mallows(n)

    # delta, delta', gamma, gamma' appear in the canonical form
    # of a Clifford operator (Eq. 3/Theorem 1)
    # delta is unit lower triangular, gamma is symmetric
    delta = Array{Int8}(LinearAlgebra.I, n, n)
    delta_p = Array{Int8}(LinearAlgebra.I, n, n)
    gamma = Array{Int8}(LinearAlgebra.Diagonal(rand(0:1, n)))
    gamma_p = Array{Int8}(LinearAlgebra.Diagonal(rand(0:1, n)))

    fill_tril(gamma, n, symmetric = true)
    fill_tril(gamma_p, n, symmetric = true)
    fill_tril(delta, n)
    fill_tril(delta_p, n)

    # now construct the tableau representation for F(I, Gamma, Delta)
    prod = gamma * delta
    prod_p = gamma_p * delta_p
    inv_delta = inv(transpose(delta))
    inv_delta_p = inv(transpose(delta_p))
    
    # block matrix form
    F1 = Array{Int8}(mod.([delta zeros(Int8, n, n); prod inv_delta], 2))
    F2 = Array{Int8}(mod.([delta_p zeros(Int8, n, n); prod_p inv_delta_p],2))

    # apply qubit permutation S to F2
    perm_inds = vcat(perm, perm .+ n)
    U = F2[perm_inds,:]
    
    # apply layer of hadamards
    had_idxs = findall(i -> hadamard[i], [1:n;])
    lhs_inds = vcat(had_idxs, had_idxs .+ n)
    rhs_inds = vcat(had_idxs .+ n, had_idxs)
    U[lhs_inds, :] = U[rhs_inds, :]
    
    # apply F1
    xzs = Array{Bool}(mod.(F1 * U,2))
    #println(xzs)
    # random Pauli matrix
    phases = Array{UInt8}(2 .* rand(0:1, 2 * n))
    return CliffordOperator(Stabilizer(phases, xzs))
end

# assign (symmetric) random ints to off diagonals of matrix
# from Qiskit
function fill_tril(matrix, n; symmetric::Bool=false)
    # Add (symmetric) random ints to off diagonals
    for row in [1:n;], col in [1:row-1;]
        b = rand(0:1)
        matrix[row, col] = b
        if symmetric
            matrix[col, row] = b
        end
    end
    matrix
end
