"""
    $TYPEDEF

The planar surface code refers to the code [kitaev2003fault](@cite) in a 2D lattice with open boundaries.
Illustration of a 3×2 surface code, where qubits are located on the edges:

```
|---1--(Z)--2---|---3---|
|  (X)  7       8       o
|---4---|---5---|---6---|
|       o       o       o
|       |       |       |
```

The surface code has open boundary conditions, unlike the toric code. To this end, we remove qubits (denoted by "o") and parity checks on the right and bottom sides.

Faces like `(1,4,7)` have `X` checks, and crosses like `(1,2,7)` have `Z` checks. Due to the removal of the bottom and right sides, we have some 3-qubit checks on the boundaries.

The parity checks of `Surface(3,2)` are:
```
X__X__X_
_X__X_XX
__X__X_X
ZZ____Z_
_ZZ____Z
___ZZ_Z_
____ZZ_Z
```

More information can be seen in [fowler2012surface](@cite).

### Fields
    $TYPEDFIELDS

"""
struct Surface <: AbstractCSSCode
    """The number of qubits in the `x` direction."""
    dx::Int
    """The number of qubits in the `z` direction."""
    dz::Int
end

code_n(c::Surface) = 2*c.dx*c.dz - c.dx -c.dz + 1

function parity_matrix_xz(c::Surface)
    tc = Toric(c.dx, c.dz)
    hx, hz = parity_matrix_xz(tc)
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

parity_matrix_x(c::Surface) = parity_matrix_xz(c)[1]
parity_matrix_z(c::Surface) = parity_matrix_xz(c)[2]

distance(c::Surface) = min(c.dx, c.dz)
