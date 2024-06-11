"""
Code concatenation of two quantum codes--the inner code c₁ and the outer code c₂.
"""
struct Concat <: AbstractECC
    c₁::AbstractECC
    c₂::AbstractECC
end

"""The function returns the parity checks of the concatenated code."""
function concat(c₁, c₂)
    k₁ = code_k(c₁)
    n₁ = code_n(c₁)
    n₂ = code_n(c₂)
    s₁ = code_s(c₁)
    s₂ = code_s(c₂)
    inner_checks = Stabilizer(vcat([embed(n₁ * n₂, 1+(i-1)*n₁:i*n₁, parity_checks(c₁)[j]) for i in 1:n₂ for j in 1:s₁]))
    h₂ = parity_matrix(c₂)
    phases₂ = phases(parity_checks(c₂))
    h_logx₁ = stab_to_gf2(logx_ops(c₁))
    phases_logx₁ = phases(logx_ops(c₁))
    h_logz₁ = stab_to_gf2(logz_ops(c₁))
    phases_logz₁ = phases(logz_ops(c₁))
    outer_check_h = transpose(hcat([vcat(
        kron(h₂[i, 1:end÷2], h_logx₁[j, 1:end÷2]) .⊻ kron(h₂[i, end÷2+1:end], h_logz₁[j, 1:end÷2]), # X part
        kron(h₂[i, 1:end÷2], h_logx₁[j, end÷2+1:end]) .⊻ kron(h₂[i, end÷2+1:end], h_logz₁[j, end÷2+1:end]) # Z part
    ) for i in 1:s₂ for j in 1:k₁]...))
    outer_check_phase = [UInt8(sum(h₂[i, 1:end÷2] * phases_logx₁[j]) + sum(h₂[i, end÷2+1:end] * phases_logz₁[j]) + phases₂[i]) & 0x3 for i in 1:s₂ for j in 1:k₁]
    outer_checks = Stabilizer(outer_check_phase, outer_check_h)
    vcat(inner_checks, outer_checks)
end

code_n(c::Concat) = code_n(c.c₁) * code_n(c.c₂)

code_k(c::Concat) = code_k(c.c₁) * code_k(c.c₂)

function parity_checks(c::Concat)
    concat(c.c₁, c.c₂)
end

function iscss(c::Concat)
    iscss(c.c₁) && iscss(c.c₂) # but it is possible that c.c₁ and c.c₂ are non-css but c is css
end
