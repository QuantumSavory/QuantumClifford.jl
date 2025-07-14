"""
    $TYPEDEF

The `[[8, 2, 3]]` color code, also known as the smallest interesting color code,
serves as the seed code for constructing Delfosse-Reichardt generalized `[[8p, 4p − 2, 3]]`
codes, as described in [delfosse2020short](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/stab_8_3_2).
"""
struct SmallestColorCode <: AbstractQECC end

code_n(c::SmallestColorCode) = 8

code_k(c::SmallestColorCode) = 2

parity_matrix(c::SmallestColorCode) = Bool[0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0;
                                           1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0;
                                           0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1;
                                           0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0;
                                           0  1  1  0  0  1  1  0  0  0  1  1  0  0  1  1;
                                           0  0  1  1  0  0  1  1  0  1  0  1  0  1  0  1]
 
distance(c::SmallestColorCode) = 3
