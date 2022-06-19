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

function H(code::Shor_code) end

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
                      ┌───┐              ┌────┐              ┌────┐   ┌───┐                ┌────┐                
    q_1: ──■────■─────┤ H ├──■─────■─────┤    ├───■─────■────┤CNOT├───┤ H ├────■───────■───┤CNOT├
           |    |     └───┘┌─┴──┐  |     |    | ┌─┴──┐  |    └─|──┘   └───┘    |       |   └─|──┘
    q_2: ──────────────────┤CNOT├────────┤    ├─┤CNOT├─────────■─────────────────────────────────
           |    |          └────┘┌─┴──┐  |    | └────┘┌─┴──┐   |               |       |     |
    q_3: ────────────────────────┤CNOT├──┤    ├───────┤CNOT├───■─────────────────────────────────
         ┌─┴──┐ |     ┌───┐      └────┘  |    |       └────┘ ┌────┐  ┌───┐   ┌─┴──┐    |     |       
    q_4: ┤CNOT├───────┤ H ├──■─────■─────┤    ├───■─────■────┤CNOT├──┤ H ├───┤CNOT├──────────■───
         └────┘ |     └───┘┌─┴──┐  |     |    | ┌─┴──┐  |    └─|──┘  └───┘   └────┘    |     |  
    q_5: ──────────────────┤CNOT├────────┤  E ├─┤CNOT├─────────■─────────────────────────────────
                |          └────┘┌─┴──┐  |    | └────┘┌─┴──┐   |                       |     |
    q_6: ────────────────────────┤CNOT├──┤    ├───────┤CNOT├───■─────────────────────────────────
              ┌─┴──┐ ┌───┐       └────┘  |    |       └────┘┌────┐  ┌───┐            ┌─┴──┐  |
    q_7: ─────┤CNOT├─┤ H ├──■─────■──────┤    ├───■─────■───┤CNOT├──┤ H ├────────────┤CNOT├──■───
              └────┘ └───┘┌─┴──┐  |      |    | ┌─┴──┐  |   └─|──┘  └───┘            └────┘
    q_8: ─────────────────┤CNOT├─────────┤    ├─┤CNOT├────────■──────────────────────────────────
                          └────┘┌─┴──┐   |    | └────┘┌─┴──┐  |                           
    q_9: ───────────────────────┤CNOT├───┤    ├───────┤CNOT├──■──────────────────────────────────
                                └────┘   └────┘       └────┘
    """

    #Step 0
    # initial_state = one(Stabilizer, N)

    #Step 1
    c1 = sCNOT(1,4)
    c2 = sCNOT(1,7)

    #Step 2
    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    #Step 3
    c3 = sCNOT(4,5)
    c4 = sCNOT(4,6)
    c5 = sCNOT(7,8)
    c6 = sCNOT(7,9)   
    
    #Step 4
    #error

    #Step 5
    c7 = sCNOT(1,2)
    c8 = sCNOT(4,5)
    c9 = sCNOT(7,8)

    #Step 6
    c10 = sCNOT(1,2)
    c11 = sCNOT(4,5)
    c12 = sCNOT(7,8)    


end # module

