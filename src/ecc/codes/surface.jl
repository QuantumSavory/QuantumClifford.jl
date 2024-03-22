"""The surface code [fowler2012surface](@cite), which is a tailored toric code with only the nearest neighboring connections in a 2D plane.

Illustration of a 3x2 surface code, where qubits are located on the edges:

|---1--(Z)--2---|---3---|
|  (X)  7       8       o
|---4---|---5---|---6---|
|       o       o       o
|       |       |       |

The surface code has an open boundary condition compared to the toric code. To this end, we remove qubits (denoted by "o") and parity checks on the right and bottom sides.

Faces like (1,4,7) have X checks, and crosses like (1,2,7) have Z checks. Due to the removal of the bottom and right sides, we have some 3-qubit checks on the boundaries.

```julia
julia> parity_checks(Surface(3,2))
+ X__X__X_
+ _X__X_XX
+ __X__X_X
+ ZZ____Z_
+ _ZZ____Z
+ ___ZZ_Z_
+ ____ZZ_Z
```

A 9-qubit surface code, specially tailored to be more compact, with a distance of three: [tomita2014low-distance](@cite).
"""
struct Surface <: AbstractECC
    dx::Int
    dz::Int
end

function iscss(::Type{Surface})
    return true
end

code_n(c::Surface) = 2*c.dx*c.dz - c.dx -c.dz + 1

function parity_checks_xz(c::Surface)
    tc = Toric(c.dx, c.dz)
    hx, hz = parity_checks_xz(tc)
    n = code_n(tc)
    nchecks = c.dx*c.dz - 1
    # remove qubits on the right and bottom sides from the toric code
    removed_qubit_indices = vcat(
        c.dx*(c.dz+1):c.dx:c.dx*(2*c.dz-1), # right side
        c.dx*(2*c.dz-1)+1:2*c.dx*c.dz) # bottom side
    qubit_indices = setdiff(1:n, removed_qubit_indices)
    # also remove the checks on the right and bottom sides
    x_check_indices = setdiff(1:nchecks, nchecks-c.dx+2:nchecks)
    z_check_indices = setdiff(1:nchecks, c.dx:c.dx:nchecks)
    hx[x_check_indices,qubit_indices], hz[z_check_indices,qubit_indices]
end

parity_checks_x(c::Surface) = parity_checks_xz(c)[1]
parity_checks_z(c::Surface) = parity_checks_xz(c)[2]

parity_checks(c::Surface) = parity_checks(CSS(parity_checks_xz(c)...))

struct Surface9 <: AbstractECC end

function iscss(::Type{Surface9})
    return true
end

code_n(c::Surface9) = 9

parity_checks(c::Surface9) = S"
                            Z__Z_____
                            _ZZ_ZZ___
                            ___ZZ_ZZ_
                            _____Z__Z
                            XX_XX____
                            _XX______
                            ____XX_XX
                            ______XX_"

parity_checks_x(c::Surface9) = stab_to_gf2(parity_checks(Surface9()))[end-1:end,1:end÷2]
parity_checks_z(c::Surface9) = stab_to_gf2(parity_checks(Surface9()))[1:end-2,end÷2+1:end]
