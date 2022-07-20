struct Shor9 <: AbstractECC end

"""The number of physical qubits in a code."""
code_n(c::Shor9) = 9

"""Parity check tableau of a code."""
parity_checks(c::Shor9) = S"ZZ_______
                            _ZZ______
                            ___ZZ____
                            ____ZZ___
                            ______ZZ_
                            _______ZZ
                            XXXXXX___
                            ___XXXXXX"

parity_matrix(c::Shor9) = stab_to_gf2(parity_checks(c::Shor9))

syndrome_circuit(c::Shor9) = #measurement circuit rows / 1 per row: conditional gates
#simple type: straight forward (naive), fault tolerant (3 types) -Neil, STeane, SHor

#Enconding circuit ----------------------------------
c1 = sCNOT(1,4)
c2 = sCNOT(1,7)

h1 = sHadamard(1)
h2 = sHadamard(4)
h3 = sHadamard(7)

c3 = sCNOT(4,5)
c4 = sCNOT(4,6)
c5 = sCNOT(7,8)
c6 = sCNOT(7,9) 

 """
Encoding shor circuit 
                  ┌───┐                           
q_1: ──■────■─────┤ H ├──■─────■─────
       |    |     └───┘┌─┴──┐  |     
q_2: ──────────────────┤CNOT├────────
       |    |          └────┘┌─┴──┐  
q_3: ────────────────────────┤CNOT├──
     ┌─┴──┐ |     ┌───┐      └────┘   
q_4: ┤CNOT├───────┤ H ├──■─────■─────
     └────┘ |     └───┘┌─┴──┐  |      
q_5: ──────────────────┤CNOT├────────
            |          └────┘┌─┴──┐  
q_6: ────────────────────────┤CNOT├──
          ┌─┴──┐ ┌───┐       └────┘  
q_7: ─────┤CNOT├─┤ H ├──■─────■──────
          └────┘ └───┘┌─┴──┐  |      
q_8: ─────────────────┤CNOT├─────────
                      └────┘┌─┴──┐                          
q_9: ───────────────────────┤CNOT├───
                            └────┘   
"""

encoding_circuit(c::Shor9) = [c1,c2,h1,h2,h3,c3,c4,c5,c6]
#----------------------------------------------------------------

code_k(c::Shor9) = 1

code_s(c::Shor9) = 8 

rate(c::Shor9) = 1/8 

distance(c::Shor9) = 3 

logx_ops(c::Shor9) = P"XXXXXXXXX"
                       
logz_ops(c::Shor9) = P"ZZZZZZZZZ"

isdegenerate(c::Shor9) = true 
