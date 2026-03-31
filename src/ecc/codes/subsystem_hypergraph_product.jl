"""
    $TYPEDEF

A Subsystem Hypergraph Product (SHP) code built from two classical parity
check matrices `H1` (m₁ × n₁) and `H2` (m₂ × n₂).

The code has `n = n₁ · n₂` physical qubits and `k = nullity(H1) · nullity(H2)`
logical qubits. X-type gauge generators are `H1 ⊗ I`, Z-type are `I ⊗ H2`.
Stabilizers are `S_X = H1 ⊗ ker(H2)` and `S_Z = ker(H1) ⊗ H2`.

Based on the construction from [quintavalle2021subsystem](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/subsystem_product).

See also: [`BravyiBaconShor`](@ref), [`SubsystemHypergraphProductSimplex`](@ref)

### Fields
    $TYPEDFIELDS
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
    
    # G_X = H1 ⊗ I_{n2}
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
    
    # G_Z = I_{n1} ⊗ H2
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

    gauge_generators = Tableau(gauge_ops)

    F = GF(2)

    # G1 = right nullspace of H1 (vectors v with H1*v=0)
    # Nemo.kernel returns left kernel (K*M=0), so transpose to get right kernel
    K1 = Nemo.kernel(transpose(matrix(F, H1)))
    G1 = zeros(Int, Nemo.nrows(K1), n1)
    for i in 1:Nemo.nrows(K1), j in 1:n1
        G1[i,j] = K1[i,j] == F(1) ? 1 : 0
    end

    # G2 = right nullspace of H2
    K2 = Nemo.kernel(transpose(matrix(F, H2)))
    G2 = zeros(Int, Nemo.nrows(K2), n2)
    for i in 1:Nemo.nrows(K2), j in 1:n2
        G2[i,j] = K2[i,j] == F(1) ? 1 : 0
    end

    stabs = PauliOperator[]

    # S_X = H1 ⊗ G2 (X-type stabilizers from rows of H1 tensored with rows of G2)
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

    # S_Z = G1 ⊗ H2 (Z-type stabilizers from rows of G1 tensored with rows of H2)
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

    stabilizer = Stabilizer(stabs)

    return SubsystemHypergraphProduct(Matrix{Int}(H1), Matrix{Int}(H2), gauge_generators, stabilizer)
end

parity_checks(c::SubsystemHypergraphProduct) = c.stabilizer

gauge_generators(c::SubsystemHypergraphProduct) = c.gauge_generators

function code_k(c::SubsystemHypergraphProduct)
    F = GF(2)
    n1 = size(c.H1, 2)
    n2 = size(c.H2, 2)
    r1 = Nemo.rank(matrix(F, c.H1))
    r2 = Nemo.rank(matrix(F, c.H2))
    return (n1 - r1) * (n2 - r2)
end

function code_g(c::SubsystemHypergraphProduct)
    return code_n(c) - code_k(c) - _stabilizer_rank(c)
end

parity_matrix_x(c::SubsystemHypergraphProduct) = _stab_to_parity_x(c.stabilizer)
parity_matrix_z(c::SubsystemHypergraphProduct) = _stab_to_parity_z(c.stabilizer)
