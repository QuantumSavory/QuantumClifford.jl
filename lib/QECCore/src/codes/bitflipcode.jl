"""
    $TYPEDEF

`Bitflip3` is a three-qubit bit-flip code that corrects single-qubit bit-flip error and does not detect any phase-flip errors.
"""
struct Bitflip3 <: AbstractQECC end

parity_matrix(c::Bitflip3) = hcat(zeros(Bool,2,3), parity_matrix(RepCode(3))[1:2,:])

"""
    $TYPEDEF

`Phaseflip3` is a three-qubit phase-flip code that corrects single-qubit phase-flip error and does not detect any bit-flip errors.
"""
struct Phaseflip3 <: AbstractQECC end

parity_matrix(c::Phaseflip3) = hcat(parity_matrix(RepCode(3))[1:2,:], zeros(Bool,2,3))
