"""
SubsystemHypergraphProductSimplex(r::Int)

Constructs the SHYPS(r) code: a Subsystem Hypergraph Product code using the
simplex code `C(r)` parity check matrix as both inputs.

Code parameters `[(2ʳ - 1)², r², 2ʳ⁻¹]`:
- Physical qubits: `n = (2ʳ - 1)²`
- Logical qubits: `k = r²`
- Distance: `d = 2ʳ⁻¹`

Simplex PCM rows are codewords of the dual Hamming code, so they have weight
at most `r+1`. Gauge generators are therefore sparse — gauge weight `O(log n)`,
making this a QLDPC code.

The simplex code's automorphism group `GL_r(2)` admits a universal set of
transversal logical gates on individual logical qubits.

Based on [malcolm2025computing](@cite).

See also: [`SubsystemHypergraphProduct`](@ref), [`Simplex`](@ref)
"""
struct SubsystemHypergraphProductSimplex <: AbstractCSSCode
    r::Int
    shp::SubsystemHypergraphProduct
end

function SubsystemHypergraphProductSimplex(r::Int)
    r >= 2 || throw(ArgumentError("`r` must be ≥ 2 to obtain a valid SHYPS code."))
    H = Matrix{Int}(QECCore.parity_matrix(QECCore.Simplex(r)))
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
