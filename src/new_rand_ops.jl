using LinearAlgebra
using Nemo: ResidueRing, MatrixSpace

const binaryring = ResidueRing(ZZ, 2)

"""Inverting a binary matrix: uses floating point for small matrices and Nemo for large matrices."""
function precise_inv(a)
    n = size(a,1)
    if n<200
        return inv(a)
    else
        M = MatrixSpace(binaryring, n, n)
        inverted = inv(M(Matrix{Int}(a))) # Nemo is very picky about input data types
        return (x->x.data).(inverted)
    end
end

"""Sample (h, S) from the distribution P_n(h, S) from Bravyi and Maslov Algorithm 1."""
function quantum_mallows(n)
    if n<500 # TODO Do in a prettier way without repetition.
        quantum_mallows_float(n)
    else
        quantum_mallows_bigint(n)
    end
end

function quantum_mallows_float(n)
    arr = collect(1:n)
    hadamard = falses(n)
    perm = zeros(Int64, n)
    for idx in 1:n
        m = length(arr)
        # sample h_i from given prob distribution
        r = rand()
        weight = Int64(2 * m - ceil(log2(r*(4.0^m-1) + 1)))
        hadamard[idx] = (weight < m)
        k = weight < m ? weight : 2*m - weight - 1
        perm[idx] = popat!(arr, k + 1) # beware of indexing in julia
    end
    return hadamard, perm
end

function quantum_mallows_bigint(n)
    arr = collect(1:n)
    hadamard = falses(n)
    perm = zeros(Int64, n)
    for idx in 1:n
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

"""Assign (symmetric) random ints to off diagonals of matrix."""
function fill_tril(matrix, n; symmetric::Bool=false)
    # Add (symmetric) random ints to off diagonals
    for row in 1:n, col in 1:row-1
        b = rand(0:1)
        matrix[row, col] = b
        if symmetric
            matrix[col, row] = b
        end
    end
    matrix
end

"""from Algorithm 2 of Bravyi and Maslov"""
function rand_clifford(n)
    hadamard, perm = quantum_mallows(n)
    had_idxs = findall(i -> hadamard[i], 1:n)
    
    # delta, delta', gamma, gamma' appear in the canonical form
    # of a Clifford operator (Eq. 3/Theorem 1)
    # delta is unit lower triangular, gamma is symmetric
    F1 = zeros(Int8, 2n, 2n)
    F2 = zeros(Int8, 2n, 2n)
    delta   = @view F1[1:n, 1:n]
    delta_p = @view F2[1:n, 1:n]
    prod   = @view F1[n+1:2n, 1:n]
    prod_p = @view F2[n+1:2n, 1:n]
    gamma   = @view F1[1:n, n+1:2n]
    gamma_p = @view F2[1:n, n+1:2n]
    inv_delta   = @view F1[n+1:2n, n+1:2n]
    inv_delta_p = @view F2[n+1:2n, n+1:2n]
    for i in 1:n
        delta[i,i] = 1
        delta_p[i,i] = 1
        gamma_p[i,i] = rand(0:1)
    end
    
    # gamma_ii is zero if h[i] = 0
    for idx in had_idxs
        gamma[idx, idx] = rand(0:1)
    end

    # gamma' and delta' are unconstrained on the lower triangular
    fill_tril(gamma_p, n, symmetric = true)
    fill_tril(delta_p, n)

    # off diagonal: gamma, delta must obey conditions C1-C5
    for row in 1:n, col in 1:row-1
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
    mul!(prod, gamma, delta)
    mul!(prod_p, gamma_p, delta_p)
    inv_delta .= mod.(precise_inv(delta'), 2)
    inv_delta_p .= mod.(precise_inv(delta_p'), 2)
 
    # block matrix form
    F1 .= mod.(F1, 2)
    F2 .= mod.(F2, 2)
    gamma .= 0
    gamma_p .= 0

    # apply qubit permutation S to F2
    perm_inds = vcat(perm, perm .+ n)
    U = F2[perm_inds,:]
    
    # apply layer of hadamards
    lhs_inds = vcat(had_idxs, had_idxs .+ n)
    rhs_inds = vcat(had_idxs .+ n, had_idxs)
    U[lhs_inds, :] = U[rhs_inds, :]
 
    # apply F1
    xzs = mod.(F1 * U,2) .== 1
 
    # random Pauli matrix just amounts to phases on the stabilizer tableau
    phases = rand([0x0,0x2], 2 * n)
    return CliffordOperator(Stabilizer(phases, xzs))
end