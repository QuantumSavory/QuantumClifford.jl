struct Shor9 <: AbstractECC end

code_n(c::Shor9) = 9

parity_checks(c::Shor9) = S"ZZ_______
                            _ZZ______
                            ___ZZ____
                            ____ZZ___
                            ______ZZ_
                            _______ZZ
                            XXXXXX___
                            ___XXXXXX"

function encoding_circuit(c::Shor9)
    c1 = sCNOT(1,4)
    c2 = sCNOT(1,7)

    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    c3 = sCNOT(1,2)
    c4 = sCNOT(4,5)
    c5 = sCNOT(7,8)

    c6 = sCNOT(1,3)
    c7 = sCNOT(4,6)
    c8 = sCNOT(7,9)

    # XXX: The extra sHadamard(1) at the start is due to a popular mismatch in
    # conventions for which logical operator is the X one and which is the Z one
    return [sHadamard(1),c1,c2,h1,h2,h3,c3,c4,c5,c6,c7,c8]
end

distance(c::Shor9) = 3
