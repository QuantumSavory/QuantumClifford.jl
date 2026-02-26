"""
    SubsystemHypergraphProduct(H1::AbstractMatrix, H2::AbstractMatrix)

Constructs a Subsystem Hypergraph Product (SHP) quantum code using two classical parity check matrices `H1` and `H2`.
The resulting code has `n_1 * n_2` physical qubits.
"""
struct SubsystemHypergraphProduct <: AbstractCSSCode
    H1::Matrix{Int}
    H2::Matrix{Int}
    gauge_generators::Tableau
    stabilizer::Stabilizer
end

function SubsystemHypergraphProduct(H1::AbstractMatrix, H2::AbstractMatrix)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    N = n1 * n2

    gauge_ops = PauliOperator[]
    
    # G_X = H1 \otimes I_{n2}
    I_n2 = Matrix{Int}(I, n2, n2)
    GX = kron(H1, I_n2) .% 2
    for r in 1:size(GX, 1)
        p = PauliOperator(0x0, falses(N), falses(N))
        empty = true
        for c in 1:N
            if GX[r, c] == 1
                p[c] = (true, false)
                empty = false
            end
        end
        if !empty push!(gauge_ops, p) end
    end
    
    # G_Z = I_{n1} \otimes H2
    I_n1 = Matrix{Int}(I, n1, n1)
    GZ = kron(I_n1, H2) .% 2
    for r in 1:size(GZ, 1)
        p = PauliOperator(0x0, falses(N), falses(N))
        empty = true
        for c in 1:N
            if GZ[r, c] == 1
                p[c] = (false, true)
                empty = false
            end
        end
        if !empty push!(gauge_ops, p) end
    end

    gauge_generators = isempty(gauge_ops) ? Tableau(zeros(PauliOperator, 0, N)) : Tableau(gauge_ops)

    F = GF(2)
    elem_type = typeof(F(1))

    # G1 = nullspace(H1)
    _, K1 = Nemo.kernel(matrix(F, H1))
    G1 = zeros(Int, Nemo.ncols(K1), n1)
    for i in 1:Nemo.ncols(K1), j in 1:n1
        G1[i,j] = K1[j,i] == F(1) ? 1 : 0
    end

    # G2 = nullspace(H2)
    _, K2 = Nemo.kernel(matrix(F, H2))
    G2 = zeros(Int, Nemo.ncols(K2), n2)
    for i in 1:Nemo.ncols(K2), j in 1:n2
        G2[i,j] = K2[j,i] == F(1) ? 1 : 0
    end

    stabs = PauliOperator[]

    # S_X = H1 \otimes G2
    # But wait, G2 has size k2 x n2, and H1 has size m1 x n1. 
    # The paper says S_X = X(h^T \otimes g) where h in row(H1) and g in row(G2).
    # This means kron(H1, G2)
    SX = kron(H1, G2) .% 2
    for r in 1:size(SX, 1)
        p = PauliOperator(0x0, falses(N), falses(N))
        empty = true
        for c in 1:N
            if SX[r, c] == 1
                p[c] = (true, false)
                empty = false
            end
        end
        if !empty push!(stabs, p) end
    end

    # S_Z = G1 \otimes H2
    SZ = kron(G1, H2) .% 2
    for r in 1:size(SZ, 1)
        p = PauliOperator(0x0, falses(N), falses(N))
        empty = true
        for c in 1:N
            if SZ[r, c] == 1
                p[c] = (false, true)
                empty = false
            end
        end
        if !empty push!(stabs, p) end
    end

    stabilizer = isempty(stabs) ? Stabilizer(zeros(PauliOperator, 0, N)) : Stabilizer(stabs)

    return SubsystemHypergraphProduct(Matrix{Int}(H1), Matrix{Int}(H2), gauge_generators, stabilizer)
end

parity_checks(c::SubsystemHypergraphProduct) = c.stabilizer

function parity_matrix_x(c::SubsystemHypergraphProduct)
    s = c.stabilizer
    # Return binary matrix of X stabilizers
    xs = falses(0, nqubits(s))
    for p in s
        if any(p.xzs[1:nqubits(s)]) && !any(p.xzs[nqubits(s)+1:end])
            xs = vcat(xs, transpose(p.xzs[1:nqubits(s)]))
        end
    end
    return xs
end

function parity_matrix_z(c::SubsystemHypergraphProduct)
    s = c.stabilizer
    # Return binary matrix of Z stabilizers
    zs = falses(0, nqubits(s))
    for p in s
        if !any(p.xzs[1:nqubits(s)]) && any(p.xzs[nqubits(s)+1:end])
            zs = vcat(zs, transpose(p.xzs[nqubits(s)+1:end]))
        end
    end
    return zs
end
