import Base: reinterpret

function reinterpret(::Type{U}, p::PauliOperator) where {U<:Unsigned}
    old = eltype(p.xz)
    total_bytes = length(p.xz) * sizeof(old)
    if total_bytes % sizeof(U) != 0
        throw(
            ArgumentError(
                "cannot reinterpret PauliOperator: backing bytes not divisible by sizeof($(U))",
            ),
        )
    end
    new_xz = Base.reinterpret(U, p.xz)
    nl = length(new_xz)
    if nl % 2 != 0
        throw(
            ArgumentError(
                "cannot reinterpret PauliOperator: resulting backing array length $(nl) is not even (cannot split into X and Z halves).",
            ),
        )
    end
    half_words = nl ÷ 2
    bits_per_half = half_words * sizeof(U) * 8
    if bits_per_half < p.nqubits
        throw(
            ArgumentError(
                lazy"""
cannot reinterpret PauliOperator: resulting backing array of $(nl) $(U) elements would provide $(bits_per_half) bits per half,
which is insufficient to represent $(p.nqubits) qubits. The backing bytes are compatible with $(U), but the resulting layout
would not contain enough bits in each X/Z half to represent the operator.
""",
            ),
        )
    end
    return PauliOperator(p.phase, p.nqubits, new_xz)
end

function reinterpret(::Type{U}, t::Tableau) where {U<:Unsigned}
    rows = size(t.xzs, 1)
    cols = size(t.xzs, 2)
    old = eltype(t.xzs)
    total_bytes = rows * sizeof(old)
    if total_bytes % sizeof(U) != 0
        throw(
            ArgumentError(
                "cannot reinterpret Tableau: backing bytes not divisible by sizeof($(U))",
            ),
        )
    end
    new_rows = (total_bytes) ÷ sizeof(U)
    if new_rows % 2 != 0
        throw(
            ArgumentError(
                "cannot reinterpret Tableau: resulting number of rows $(new_rows) is not even (cannot split into X and Z halves).",
            ),
        )
    end
    half_rows = new_rows ÷ 2
    bits_per_half = half_rows * sizeof(U) * 8
    if bits_per_half < t.nqubits
        throw(
            ArgumentError(
                lazy"""
cannot reinterpret Tableau: resulting backing array of $(new_rows) $(U) rows would provide $(bits_per_half) bits per half,
which is insufficient to represent $(t.nqubits) qubits. The backing bytes are compatible with $(U), but the resulting layout
would not contain enough bits in each X/Z half to represent the tableau.
""",
            ),
        )
    end
    raw = vec(t.xzs)
    new_vec = Base.reinterpret(U, raw)
    new_xzs = reshape(new_vec, new_rows, cols)
    return Tableau(t.phases, t.nqubits, new_xzs)
end

reinterpret(::Type{U}, s::Stabilizer) where {U<:Unsigned} =
    Stabilizer(reinterpret(U, tab(s)))

function reinterpret(::Type{U}, f::PauliFrame) where {U<:Unsigned}
    return PauliFrame(reinterpret(U, f.frame), f.measurements)
end
