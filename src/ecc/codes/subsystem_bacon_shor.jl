"""
    $TYPEDEF

A Bravyi-Bacon-Shor (BBS) subsystem code defined by a binary matrix `A`.
Physical qubits sit at the nonzero entries of `A`, giving `n = nnz(A)`.
The number of logical qubits equals the GF(2)-rank of `A`.

X-type gauge generators are weight-2 operators pairing consecutive nonzero
entries in each column of `A`; Z-type gauge generators pair entries in each row.
Stabilizers are obtained from the left and right GF(2)-kernels of `A`.

Based on the construction from [bravyi2011subsystem](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/bacon_shor).

See also: [`SubsystemHypergraphProduct`](@ref)

### Fields
    $TYPEDFIELDS
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
    
    gauge_generators = Tableau(gauge_ops)

    # Compute Stabilizers using Nemo.kernel
    # Nemo.kernel returns left kernel (K*M=0), so transpose to get right kernel
    # H2 = right kernel of A (vectors v with A*v=0)
    # H1 = left kernel of A (vectors v with v*A=0)
    F = GF(2)
    A_F = matrix(F, A)
    K2 = Nemo.kernel(transpose(A_F))
    H2 = zeros(Int, Nemo.nrows(K2), nc)
    for i in 1:Nemo.nrows(K2), j in 1:nc
        H2[i,j] = K2[i,j] == F(1) ? 1 : 0
    end

    K1 = Nemo.kernel(A_F)
    H1 = zeros(Int, Nemo.nrows(K1), nr)
    for i in 1:Nemo.nrows(K1), j in 1:nr
        H1[i,j] = K1[i,j] == F(1) ? 1 : 0
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

    stabilizer = Stabilizer(stabs)

    return BravyiBaconShor(A, gauge_generators, stabilizer)
end

parity_checks(c::BravyiBaconShor) = c.stabilizer

gauge_generators(c::BravyiBaconShor) = c.gauge_generators

function code_k(c::BravyiBaconShor)
    F = GF(2)
    return Nemo.rank(matrix(F, c.A))
end

function code_g(c::BravyiBaconShor)
    return code_n(c) - code_k(c) - _stabilizer_rank(c)
end

parity_matrix_x(c::BravyiBaconShor) = _stab_to_parity_x(c.stabilizer)
parity_matrix_z(c::BravyiBaconShor) = _stab_to_parity_z(c.stabilizer)
