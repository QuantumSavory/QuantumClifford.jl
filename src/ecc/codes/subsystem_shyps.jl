"""
    SubsystemHypergraphProductSimplex(r::Int)

Constructs the SHYPS(r) code — a Subsystem Hypergraph Product code instantiated
with the parity check matrix of the classical simplex code `C(r)`.

The SHYPS(r) code has parameters `[(2ʳ - 1)², r², 2ʳ⁻¹]`, where:
- Physical qubits: `n = (2ʳ - 1)²`
- Logical qubits: `k = r²`
- Distance: `d = 2ʳ⁻¹`

The gauge generators have weight at most 3 (inherited from the simplex code's
constant row weight), making SHYPS codes quantum LDPC (QLDPC).

SHYPS codes support a universal set of transversal logical gates on individual
logical qubits due to the large automorphism group `GL_r(2)` of the simplex code.

Based on the construction from [malcolm2025computing](@cite).

See also: [`SubsystemHypergraphProduct`](@ref), [`Simplex`](@ref)
"""
struct SubsystemHypergraphProductSimplex <: AbstractCSSCode
    r::Int
    shp::SubsystemHypergraphProduct
end

function SubsystemHypergraphProductSimplex(r::Int)
    r >= 2 || throw(ArgumentError("`r` must be ≥ 2 to obtain a valid SHYPS code."))
    H = Matrix{Int}(parity_matrix(Simplex(r)))
    shp = SubsystemHypergraphProduct(H, H)
    return SubsystemHypergraphProductSimplex(r, shp)
end

parity_checks(c::SubsystemHypergraphProductSimplex) = parity_checks(c.shp)

parity_matrix_x(c::SubsystemHypergraphProductSimplex) = parity_matrix_x(c.shp)

parity_matrix_z(c::SubsystemHypergraphProductSimplex) = parity_matrix_z(c.shp)

code_n(c::SubsystemHypergraphProductSimplex) = (2^c.r - 1)^2

code_k(c::SubsystemHypergraphProductSimplex) = c.r^2

distance(c::SubsystemHypergraphProductSimplex) = 2^(c.r - 1)

gauge_generators(c::SubsystemHypergraphProductSimplex) = gauge_generators(c.shp)

function code_g(c::SubsystemHypergraphProductSimplex)
    return code_n(c) - code_k(c) - _stabilizer_rank(c)
end
