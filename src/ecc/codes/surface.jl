"""The surface code."""
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

function iscss(::Type{Surface})
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
