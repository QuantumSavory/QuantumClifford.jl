struct LPCode{T} <: AbstractECC
    c₁::LiftedCode{T}
    c₂::LiftedCode{T}
    repr::Union{Function, Nothing} # nothing means using repr from each code respectively

    function LPCode(c₁::LiftedCode{T}, c₂::LiftedCode{T}, repr::Function) where T
        new{T}(c₁, c₂, repr)
    end

    function LPCode(c₁::LiftedCode{T}, c₂::LiftedCode{T}) where T
        new{T}(c₁, c₂, nothing)
    end
end

function LPCode(c₁::LiftedCode{PermGroupRingElem{FqFieldElem}}, c₂::LiftedCode{PermGroupRingElem{FqFieldElem}})
    LPCode(c₁, c₂, permutation_repr)
end

iscss(::Type{LPCode{T}}) where {T <: RingElement} = true # TODO wrong

# In most cases, especially for `PermGroupRingElem`,
# "first product then lift" will be much more efficient,
# which requires the same matrix representation for both codes.
function parity_checks_xz(c::LPCode)
    if !isnothing(c.repr) # "first product then lift"
        c.c₁.A[1,1].parent == c.c₂.A[1,1].parent || error("The base rings of the two codes must be the same")
        hx, hz = hgp(c.c₁.A, c.c₂.A)
        # @show hx, hz
        hx, hz = lift(c.repr, hx), lift(c.repr, hz)
    else # fall back to the slow "first lift then product"
        h₁ = lift(c.c₁.repr, parity_checks(c.c₁))
        h₂ = lift(c.c₂.repr, parity_checks(c.c₂))
        hx, hz =  hgp(h₁, h₂)
    end
    return hx, hz
end

parity_checks_x(c::LPCode) = parity_checks_xz(c)[1]

parity_checks_z(c::LPCode) = parity_checks_xz(c)[2]

parity_checks(c::LPCode) = parity_checks(CSS(parity_checks_xz(c)...))

# HGPCode = LPCode{Bool}
# or from a trivial group ring, `R = PermutationGroupRing(GF(2), 1)`
