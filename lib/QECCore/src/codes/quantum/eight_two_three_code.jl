"""
    $TYPEDEF

The `[[8, 2, 3]]` color code serves as the seed code for constructing Delfosse-Reichardt
generalized `[[8p, 4p − 2, 3]]`codes, as described in [delfosse2020short](@cite).
"""
struct EightTwoThree <: AbstractQECC end

code_n(c::EightTwoThree) = 8

code_k(c::EightTwoThree) = 2

parity_matrix(c::EightTwoThree) = Bool[0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0;
                                       1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0;
                                       0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1;
                                       0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0;
                                       0  1  1  0  0  1  1  0  0  0  1  1  0  0  1  1;
                                       0  0  1  1  0  0  1  1  0  1  0  1  0  1  0  1]
 
distance(c::EightTwoThree) = 3
