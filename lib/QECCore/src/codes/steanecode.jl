# TODO Steane5 qubit code

"""
    Steane7 <: AbstractCSSCode

Steane code [steane1996error](@cite).
"""
struct Steane7 <: AbstractCSSCode end

_steane_mat() = Bool[0 0 0 1 1 1 1;
                    0 1 1 0 0 1 1;
                    1 0 1 0 1 0 1]

parity_matrix_x(c::Steane7) = _steane_mat()
parity_matrix_z(c::Steane7) = _steane_mat()

distance(c::Steane7) = 3
