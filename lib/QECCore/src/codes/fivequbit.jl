"""
    $TYPEDEF

`Perfect5` code [lafiamme1996perfect](@cite) is the smallest qubit stabilizer code to correct a single-qubit error.
"""
struct Perfect5 <: AbstractQECC end

parity_matrix(c::Perfect5) = Bool[1 0 0 1 0 0 1 1 0 0;
                                0 1 0 0 1 0 0 1 1 0;
                                1 0 1 0 0 0 0 0 1 1;
                                0 1 0 1 0 1 0 0 0 1]
distance(c::Perfect5) = 3
