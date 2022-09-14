struct Steane5 <: AbstractECC end
struct Steane7 <: AbstractECC end

#Steane 5 qubit code
parity_checks(c::Steane5) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

code_n(c::Steane5) = 5

#Encoding circuit ----------------------------------
#https://www.researchgate.net/figure/Encoding-circuit-for-the-five-qubit-code-a-Circuit-to-encode-the-logical-minus-state_fig1_337273308

function encoding_circuit(c::Steane5) #black, white
    
    z1 = sZ(1)
    h1 = sHadamard(3)
    h2 = sHadamard(4)
    is1 = sInvPhase(1)
    c1 = sCNOT(3,5)
    c2 = sCNOT(4,2)
    h3 = sHadamard(2)
    c3 = sCNOT(4,5)
    c4 = sCNOT(2,1)
    is2 = sInvPhase(3)
    s2 = sPhase(4)
    is3 = sInvPhase(5)
    s2 = sPhase(1)
    s3 = sPhase(2)
    z2 = sZ(3)
    c5 = sCNOT(5,1)
    h4 = sHadamard(5)
    c6 = sCNOT(5,2)

    return  [z1,h1,h2,is1,c1,c2,h3,c3,c4,is2,s2,is3,s2,s3,z2,c5,h4,c6] 
    
end

#----------------------------------------------------------------

distance(c::Steane5) = 3

isdegenerate(c::Steane5) = false 

##----------------------------------------------------------------------------------------------------------------------------

#steane 7 qubit code
parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

code_n(c::Steane7) = 7

#Encoding circuit ----------------------------------
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
#----------------------------------------------------------------

distance(c::Steane7) = 3

isdegenerate(c::Steane7) = false 
