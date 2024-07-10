import Nemo: matrix_space

struct LPCode <: AbstractECC
    c₁::LiftedCode
    c₂::LiftedCode
    repr::Function

    function LPCode(c₁::LiftedCode, c₂::LiftedCode, repr::Function)
        new(c₁, c₂, repr)
    end
end

function LPCode(c₁::LiftedCode, c₂::LiftedCode)
    LPCode(c₁, c₂, permutation_repr)
end

iscss(::Type{LPCode}) = true

function parity_checks_xz(c::LPCode)
    c.c₁.A[1,1].parent == c.c₂.A[1,1].parent || error("The base rings of the two codes must be the same")
    hx, hz = hgp(c.c₁.A, c.c₂.A)
    hx, hz = lift(c.repr, hx), lift(c.repr, hz)
    return hx, hz
end

parity_checks_x(c::LPCode) = parity_checks_xz(c)[1]

parity_checks_z(c::LPCode) = parity_checks_xz(c)[2]

parity_checks(c::LPCode) = parity_checks(CSS(parity_checks_xz(c)...))

code_n(c::LPCode) = size(c.repr(parent(c.c₁.A[1, 1])(0)), 2) * (size(c.c₁.A, 2) * size(c.c₂.A, 2) + size(c.c₁.A, 1) * size(c.c₂.A, 1))

function code_s(c::LPCode)
    hx, hz = parity_checks_xz(c)
    mod2rank(hx) + mod2rank(hz)
end

code_k(c::LPCode) = code_n(c) - code_s(c)

# HGPCode = LPCode{Bool}
# or from a trivial group ring, `R = PermutationGroupRing(GF(2), 1)`
