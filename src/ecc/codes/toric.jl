"""The Toric code."""
struct Toric <: AbstractECC
    dx::Int
    dz::Int
end

function iscss(::Type{Toric})
    return true
end

code_n(c::Toric) = 2*c.dx*c.dz

function parity_checks_xz(c::Toric)
    hx, hz = hgp(parity_checks(RepCode(c.dz)), parity_checks(RepCode(c.dx)))
    hx[1:end-1,:], hz[1:end-1,:]
end

parity_checks_x(c::Toric) = parity_checks_xz(c)[1]
parity_checks_z(c::Toric) = parity_checks_xz(c)[2]

parity_checks(c::Toric) = parity_checks(CSS(parity_checks_xz(c)...))
