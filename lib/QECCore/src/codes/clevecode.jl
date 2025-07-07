"""
    $TYPEDEF

A pedagogical example of a quantum error correcting [8,3] code used in [cleve1997efficient](@cite).
"""
struct Cleve8 <: AbstractQECC end

parity_matrix(c::Cleve8) = Bool[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                                1 0 1 0 0 1 0 1 0 0 0 0 1 1 1 1;
                                1 0 1 0 1 0 1 0 0 0 1 1 0 0 1 1;
                                1 0 0 1 0 1 1 0 0 1 0 1 0 1 0 1]