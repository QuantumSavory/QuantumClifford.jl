"""A pedagogical example of a quantum error correcting [8,3] code used in [cleve1997efficient](@cite)."""
struct Cleve8 <: AbstractECC end

code_n(c::Cleve8) = 8
code_k(c::Cleve8) = 3

parity_checks(c::Cleve8) = S"XXXXXXXX
                             ZZZZZZZZ
                             XIXIZYZY
                             XIYZXIYZ
                             XZIYIYXZ"

function encoding_circuit(c::Cleve8)
    c1 = sCNOT(1,4)
    c2 = sCNOT(2,4)
    c3 = sCNOT(3,4)
    h1 = sHadamard(5)
    z1 = sZ(5)
    h2 = sHadamard(6)
    h3 = sHadamard(7)
    z3 = sZ(7)
    h4 = sHadamard(8)
    z4 = sZ(8)
    first_part = [c1,c2,c3,h1,z1,h2,h3,z3,h4,z4]

    c1 = sZCX(5, 1)
    c2 = sZCX(5, 2)
    c3 = sZCY(5, 3)
    c4 = sZCZ(5, 4)
    c5 = sZCZ(5, 6)
    column1 = [c1,c2,c3,c4,c5] # 1st non null column of Zstar

    c1 = sZCX(6, 1)
    c2 = sZCY(6, 2)
    c3 = sZCY(6, 4)
    c4 = sZCZ(6, 5)
    c5 = sZCZ(6, 7)
    column2 = [c1,c2,c3,c4,c5]

    c1 = sZCX(7, 1)
    c2 = sZCY(7, 3)
    c3 = sZCX(7, 4)
    c4 = sZCZ(7, 5)
    c5 = sZCZ(7, 8)
    column3 = [c1,c2,c3,c4,c5]

    c1 = sZCY(8, 2)
    c2 = sZCX(8, 3)
    c3 = sZCX(8, 4)
    c4 = sZCZ(8, 5)
    c5 = sZCZ(8, 6)
    column4 = [c1,c2,c3,c4,c5]

    return vcat(first_part, column1, column2, column3, column4)
end