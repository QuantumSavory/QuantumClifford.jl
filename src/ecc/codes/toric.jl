"""The Toric code [kitaev2003fault](@cite).

Illustration of a 2x2 toric code, where qubits are located on the edges:

```
|--1-(Z)-2--|
| (X) 5     6
|--3--|--4--|
|     7     8
|     |     |
```

It is important to note that the toric code has periodic boundary conditions, which means that the top and bottom sides are essentially glued together, as are the left and right sides.

Faces like `(1,3,5,6)` have X checks, and crosses like `(1,2,5,7)` have Z checks.

```jldoctest
julia> parity_checks(Toric(2,2))
+ X_X_XX__
+ _X_XXX__
+ X_X___XX
+ ZZ__Z_Z_
+ ZZ___Z_Z
+ __ZZZ_Z_
```
"""
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
