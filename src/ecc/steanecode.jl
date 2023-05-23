# TODO Steane5 qubit code

struct Steane7 <: AbstractECC end

parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

code_n(c::Steane7) = 7

function encoding_circuit(c::Steane7)
    sc1 = sCNOT(1,2)
    sc2 = sCNOT(1,3)

    sh1 = sHadamard(5)
    sh2 = sHadamard(6)
    sh3 = sHadamard(7)

    sc3 = sCNOT(7,4)
    sc4 = sCNOT(7,2)
    sc5 = sCNOT(7,1)
    sc6 = sCNOT(6,4)
    sc7 = sCNOT(6,3)
    sc8 = sCNOT(6,1)
    sc9 = sCNOT(5,4)
    sc10 = sCNOT(5,3)
    sc11 = sCNOT(5,2)

    return [sc1,sc2,sh1,sh2,sh3,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10,sc11]
end

distance(c::Steane7) = 3
