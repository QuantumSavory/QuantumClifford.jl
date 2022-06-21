#currently in the works

module Codes
abstract type Code end

struct Shor_code <: Code
end

"""documents"""

function rate end

function logicalqubits end

function physicalqubits end

function codedistance end

#function H(code::Shor_code) end

function QuantumClifford.MixedDestabilizer(code::Shor_code)

    #pauli matrices of the code 
    Stabilizer([P"X",P"Y",P"Z"])

    Stabilizer(Bool[0 1;
                    1 0],
                Bool[0 -i;
                    i 1],
                Bool[1 0;
                    0 -1])
end

...

function rate(code::Shor_code) return 1//9 end
    N= 9 #n qubits 

           #Step 1   #Step 2 #Step 3     #Step 4 #S 5  #S 6  #Step 7  #Step 8 #Step 9  #S 10 #S 11
    """
                      ┌───┐              ┌─────┐              ┌────┐   ┌───┐                ┌────┐                
    q_1: ──■────■─────┤ H ├──■─────■─────┤     ├───■─────■────┤CNOT├───┤ H ├────■───────■───┤CNOT├
           |    |     └───┘┌─┴──┐  |     |     | ┌─┴──┐  |    └─|──┘   └───┘    |       |   └─|──┘
    q_2: ──────────────────┤CNOT├────────┤     ├─┤CNOT├─────────■─────────────────────────────────
           |    |          └────┘┌─┴──┐  |     | └────┘┌─┴──┐   |               |       |     |
    q_3: ────────────────────────┤CNOT├──┤     ├───────┤CNOT├───■─────────────────────────────────
         ┌─┴──┐ |     ┌───┐      └────┘  |     |       └────┘ ┌────┐   ┌───┐  ┌─┴──┐    |     |       
    q_4: ┤CNOT├───────┤ H ├──■─────■─────┤     ├───■─────■────┤CNOT├───┤ H ├──┤CNOT├──────────■───
         └────┘ |     └───┘┌─┴──┐  |     |     | ┌─┴──┐  |    └─|──┘   └───┘  └────┘    |     |  
    q_5: ──────────────────┤CNOT├────────┤  E  ├─┤CNOT├─────────■─────────────────────────────────
                |          └────┘┌─┴──┐  |     | └────┘┌─┴──┐   |                       |     |
    q_6: ────────────────────────┤CNOT├──┤     ├───────┤CNOT├───■─────────────────────────────────
              ┌─┴──┐ ┌───┐       └────┘  |     |       └────┘┌────┐    ┌───┐          ┌─┴──┐  |
    q_7: ─────┤CNOT├─┤ H ├──■─────■──────┤     ├───■─────■───┤CNOT├────┤ H ├──────────┤CNOT├──■───
              └────┘ └───┘┌─┴──┐  |      |     | ┌─┴──┐  |   └─|──┘    └───┘          └────┘
    q_8: ─────────────────┤CNOT├─────────┤     ├─┤CNOT├────────■──────────────────────────────────
                          └────┘┌─┴──┐   |     | └────┘┌─┴──┐  |                           
    q_9: ───────────────────────┤CNOT├───┤     ├───────┤CNOT├──■──────────────────────────────────
                                └────┘   └─────┘       └────┘
    """

    #Step 0
    # initial_state = one(Stabilizer, N)

    #Step 1: 1st set of  CNOT gates
    c1 = sCNOT(1,4)
    c2 = sCNOT(1,7)

    #Step 2: 1st set of  Hadamard gates
    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    #Step 3: 2nd set of  CNOT gates
    c3 = sCNOT(4,5)
    c4 = sCNOT(4,6)
    c5 = sCNOT(7,8)
    c6 = sCNOT(7,9)   
    
    #Step 4: Error
    #X: Bit flip error
    #Z: Phase flip error
    for qubit in range(N) #check this
        X = S"X"
        Z = S"Z"
        X(qubit)
        Z(qubit)

    #Step 5: 4th set of  CNOT GATES
    c7 = sCNOT(1,2)
    c8 = sCNOT(4,5)
    c9 = sCNOT(7,8)

    #Step 6: 5th set of  CNOT gates
    c10 = sCNOT(1,2)
    c11 = sCNOT(4,5)
    c12 = sCNOT(7,8)
    
    #Step 7: 1st set of  Toffoli gates
    #are toffoli gates represented by sCCNOT?
    cc1 = sCCNOT(2,3,1)
    cc2 = sCCNOT(5,6,4)
    cc3 = sCCNOT(9,8,7)

    #Step 8: 2nd set of  Haramard gates
    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    #Step 9: 6th set of  CNOT gates
    c13 = sCNOT(1,4)

    #Step 10: 7th set of CNOT gates
    c14 = sCNOT(1,7)

    #Step 11: 2nd set of Toffoli gates
    #Final gates
    cc4 = sCCNOT(7,4,1)

# This circuit performs a depolarization at rate `epsilon` to all qubits,
circuit = [c2,c2,h1,h2,h3,c3,c4,c5,c6,X,Z,c7,c8,c9,c10,c11,c12,cc1,cc2,cc3,h1,h2,h3,,c13,c14,cc4]

end # module