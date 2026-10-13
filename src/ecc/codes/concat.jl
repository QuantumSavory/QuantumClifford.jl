"""
`Concat(c₁, c₂)` is a code concatenation of two quantum codes [knill1996concatenated](@cite).

The inner code c₁ and the outer code c₂.
The construction is the following: replace each qubit in code c₂ with logical qubits encoded by code c₁.
The resulting code will have `n = n₁ × n₂` qubits and `k = k₁ × k₂` logical qubits.
"""
struct Concat <: AbstractQECC
    c₁::AbstractQECC
    c₂::AbstractQECC
end

function parity_checks(c::Concat)
    c₁ = c.c₁
    c₂ = c.c₂
    k₁ = code_k(c₁)
    n₁ = code_n(c₁)
    n₂ = code_n(c₂)
    s₁ = code_s(c₁)
    s₂ = code_s(c₂)
    inner_checks = Stabilizer(vcat([embed(n₁ * n₂, 1+(i-1)*n₁:i*n₁, parity_checks(c₁)[j]) for i in 1:n₂ for j in 1:s₁])) # parity checks of c₁ on each qubit of c₂
    h₂ = parity_matrix(c₂)
    phases₂ = phases(parity_checks(c₂))
    h_logx₁ = stab_to_gf2(logx_ops(c₁))
    phases_logx₁ = phases(logx_ops(c₁))
    h_logz₁ = stab_to_gf2(logz_ops(c₁))
    phases_logz₁ = phases(logz_ops(c₁))
    # parity checks of c₂ with qubits replaced with logical qubits of c₁
    outer_check_h = transpose(hcat([vcat(
        kron(h₂[i, 1:end÷2], h_logx₁[j, 1:end÷2]) .⊻ kron(h₂[i, end÷2+1:end], h_logz₁[j, 1:end÷2]), # X part
        kron(h₂[i, 1:end÷2], h_logx₁[j, end÷2+1:end]) .⊻ kron(h₂[i, end÷2+1:end], h_logz₁[j, end÷2+1:end]) # Z part
    ) for i in 1:s₂ for j in 1:k₁]...))
    outer_check_phase = [UInt8(sum(h₂[i, 1:end÷2] * phases_logx₁[j]) + sum(h₂[i, end÷2+1:end] * phases_logz₁[j]) + phases₂[i]) & 0x3 for i in 1:s₂ for j in 1:k₁]
    outer_checks = Stabilizer(outer_check_phase, outer_check_h)
    vcat(inner_checks, outer_checks)
end

function _concat_parity_matrix_xz(c::Concat)
    c₁ = c.c₁
    c₂ = c.c₂
    n₁ = code_n(c₁)
    n₂ = code_n(c₂)
    logicals₁ = MixedDestabilizer(parity_checks(c₁))
    h_logx₁ = stab_to_gf2(logicalxview(logicals₁))[:, 1:n₁]
    h_logz₁ = stab_to_gf2(logicalzview(logicals₁))[:, n₁+1:2n₁]
    hx = vcat(
        kron(I(n₂), parity_matrix_x(c₁)),
        kron(parity_matrix_x(c₂), h_logx₁),
    )
    hz = vcat(
        kron(I(n₂), parity_matrix_z(c₁)),
        kron(parity_matrix_z(c₂), h_logz₁),
    )
    return hx, hz
end

parity_matrix_x(c::Concat) = _concat_parity_matrix_xz(c)[1]

parity_matrix_z(c::Concat) = _concat_parity_matrix_xz(c)[2]

code_n(c::Concat) = code_n(c.c₁) * code_n(c.c₂)

code_k(c::Concat) = code_k(c.c₁) * code_k(c.c₂)

function iscss(c::Concat)
    if iscss(c.c₁)==true && iscss(c.c₂)==true # to distinguish from potentially being nothing
        return true
    end
    return nothing # if c.c₁ or c.c₂ are non-CSS; in this case, `Concat(c₁, c₂)` can still be CSS
end
