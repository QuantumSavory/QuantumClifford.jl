"""
    BravyiBaconShor(A::AbstractMatrix)

Constructs a Bravyi-Bacon-Shor subsystem quantum code using a binary matrix `A`.
The code qubits are located at the non-zero entries of `A`.
"""
struct BravyiBaconShor <: AbstractCSSCode
    A::Matrix{Int}
    gauge_generators::Tableau
    stabilizer::Stabilizer
end

function BravyiBaconShor(A::AbstractMatrix)
    nr, nc = size(A)
    # Find active qubits
    qubits = Tuple{Int,Int}[]
    for i in 1:nr, j in 1:nc
        if A[i,j] != 0
            push!(qubits, (i,j))
        end
    end
    N = length(qubits)
    idx_map = Dict(q => k for (k,q) in enumerate(qubits))

    gauge_ops = PauliOperator[]
    
    # X gauges: column
    for j in 1:nc
        cols_i = [i for i in 1:nr if A[i,j] != 0]
        for k in 1:length(cols_i)-1
            p = PauliOperator(0x0, falses(N), falses(N))
            p[idx_map[(cols_i[k], j)]] = (true, false)
            p[idx_map[(cols_i[k+1], j)]] = (true, false)
            push!(gauge_ops, p)
        end
    end
    
    # Z gauges: row
    for i in 1:nr
        rows_j = [j for j in 1:nc if A[i,j] != 0]
        for k in 1:length(rows_j)-1
            p = PauliOperator(0x0, falses(N), falses(N))
            p[idx_map[(i, rows_j[k])]] = (false, true)
            p[idx_map[(i, rows_j[k+1])]] = (false, true)
            push!(gauge_ops, p)
        end
    end
    
    gauge_generators = isempty(gauge_ops) ? Tableau(zeros(PauliOperator, 0, N)) : Tableau(gauge_ops)

    # Compute Stabilizers using Nemo.kernel
    # C2 = row(A), parity check H2 is right kernel of A
    # C1 = col(A), parity check H1 is right kernel of A^T (left kernel of A)
    F = GF(2)
    elem_type = typeof(F(1))
    A_F = matrix(F, A)
    _, K2 = Nemo.kernel(A_F)
    H2 = zeros(Int, Nemo.ncols(K2), nc)
    for i in 1:Nemo.ncols(K2), j in 1:nc
        H2[i,j] = K2[j,i] == F(1) ? 1 : 0
    end

    _, K1 = Nemo.kernel(transpose(A_F))
    H1 = zeros(Int, Nemo.ncols(K1), nr)
    for i in 1:Nemo.ncols(K1), j in 1:nr
        H1[i,j] = K1[j,i] == F(1) ? 1 : 0
    end

    stabs = PauliOperator[]

    # S_X = {X(diag(r) * A) : r in row(H1)}
    for r in 1:size(H1, 1)
        p = PauliOperator(0x0, falses(N), falses(N))
        empty = true
        for i in 1:nr
            if H1[r, i] == 1
                for j in 1:nc
                    if A[i,j] != 0
                        p[idx_map[(i,j)]] = (true, false)
                        empty = false
                    end
                end
            end
        end
        if !empty push!(stabs, p) end
    end

    # S_Z = {Z(A * diag(c)) : c in row(H2)}
    for c in 1:size(H2, 1)
        p = PauliOperator(0x0, falses(N), falses(N))
        empty = true
        for j in 1:nc
            if H2[c, j] == 1
                for i in 1:nr
                    if A[i,j] != 0
                        p[idx_map[(i,j)]] = (false, true)
                        empty = false
                    end
                end
            end
        end
        if !empty push!(stabs, p) end
    end

    stabilizer = isempty(stabs) ? Stabilizer(zeros(PauliOperator, 0, N)) : Stabilizer(stabs)

    return BravyiBaconShor(A, gauge_generators, stabilizer)
end

parity_checks(c::BravyiBaconShor) = c.stabilizer

function parity_matrix_x(c::BravyiBaconShor)
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

function parity_matrix_z(c::BravyiBaconShor)
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
