"""
    $TYPEDEF

A Bravyi-Bacon-Shor (BBS) subsystem code built from a binary matrix `A`.
Each nonzero entry of `A` is a physical qubit, so `n = nnz(A)`.
The number of logical qubits equals the GF(2) rank of `A`.

X-type gauge generators connect consecutive qubit pairs down each column of `A`;
Z-type do the same across each row. Stabilizers come from the left and right
null spaces of `A` over GF(2).

Based on [li2020numerical](@cite).
ECC Zoo: [Bacon-Shor code family](https://errorcorrectionzoo.org/c/bacon_shor).

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
    # find which positions in A actually have qubits (nonzero entries only)
    qubits = Tuple{Int,Int}[]
    for i in 1:nr, j in 1:nc
        if A[i,j] != 0
            push!(qubits, (i,j))
        end
    end
    N = length(qubits)
    idx_map = Dict(q => k for (k,q) in enumerate(qubits))

    gauge_ops = PauliOperator[]
    
    # X-type gauge generators -- pair up consecutive qubits going down each column
    for j in 1:nc
        cols_i = [i for i in 1:nr if A[i,j] != 0]
        for k in 1:length(cols_i)-1
            p = PauliOperator(0x0, falses(N), falses(N))
            p[idx_map[(cols_i[k], j)]] = (true, false)
            p[idx_map[(cols_i[k+1], j)]] = (true, false)
            push!(gauge_ops, p)
        end
    end
    
    # Z-type gauge generators -- same idea, but pair consecutive qubits across each row
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

    # now compute the stabilizers -- we need both kernels of A over GF(2)
    # Nemo.kernel gives the LEFT kernel (K*M=0), so we transpose A to get the right kernel
    # H2 = right kernel of A, i.e. vectors v where A*v = 0
    # H1 = left kernel of A, i.e. vectors v where v*A = 0
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

    # build X stabilizers -- one per row of H1 (the left kernel of A)
    # each row picks a subset of rows of A; we put X on all qubits in those rows
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

    # build Z stabilizers -- same idea but using H2 (right kernel of A)
    # each row of H2 picks a subset of columns; we put Z on all qubits in those columns
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
