"""
    $TYPEDEF

Shor9 code [shor1995scheme](@cite)is a nine-qubit CSS code that is the first quantum error-correcting code.
"""
struct Shor9 <: AbstractCSSCode end

parity_matrix_x(c::Shor9) = Bool[1 1 1 1 1 1 0 0 0 ;
                                0 0 0 1 1 1 1 1 1 ]

parity_matrix_z(c::Shor9) = Bool[1 1 0 0 0 0 0 0 0;
                                0 1 1 0 0 0 0 0 0;
                                0 0 0 1 1 0 0 0 0;
                                0 0 0 0 1 1 0 0 0;
                                0 0 0 0 0 0 1 1 0;
                                0 0 0 0 0 0 0 1 1]

distance(c::Shor9) = 3