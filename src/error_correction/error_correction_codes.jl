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

    """
                      ┌───┐              ┌────┐              ┌────┐   ┌───┐                ┌────┐                
    q_0: ──■────■─────┤ H ├──■─────■─────┤    ├───■─────■────┤CNOT├───┤ H ├────■───────■───┤CNOT├
           |    |     └───┘┌─┴──┐  |     |    | ┌─┴──┐  |    └─|──┘   └───┘    |       |   └─|──┘
    q_1: ──────────────────┤CNOT├────────┤    ├─┤CNOT├─────────■─────────────────────────────────
           |    |          └────┘┌─┴──┐  |    | └────┘┌─┴──┐   |               |       |     |
    q_2: ────────────────────────┤CNOT├──┤    ├───────┤CNOT├───■─────────────────────────────────
         ┌─┴──┐ |     ┌───┐      └────┘  |    |       └────┘ ┌────┐  ┌───┐   ┌─┴──┐    |     |       
    q_3: ┤CNOT├───────┤ H ├──■─────■─────┤    ├───■─────■────┤CNOT├──┤ H ├───┤CNOT├──────────■───
         └────┘ |     └───┘┌─┴──┐  |     |    | ┌─┴──┐  |    └─|──┘  └───┘   └────┘    |     |  
    q_4: ──────────────────┤CNOT├────────┤ E  ├─┤CNOT├─────────■─────────────────────────────────
                |          └────┘┌─┴──┐  |    | └────┘┌─┴──┐   |                       |     |
    q_5: ────────────────────────┤CNOT├──┤    ├───────┤CNOT├───■─────────────────────────────────
              ┌─┴──┐ ┌───┐       └────┘  |    |       └────┘┌────┐  ┌───┐            ┌─┴──┐  |
    q_6: ─────┤CNOT├─┤ H ├──■─────■──────┤    ├───■─────■───┤CNOT├──┤ H ├────────────┤CNOT├──■───
              └────┘ └───┘┌─┴──┐  |      |    | ┌─┴──┐  |   └─|──┘  └───┘            └────┘
    q_7: ─────────────────┤CNOT├─────────┤    ├─┤CNOT├────────■──────────────────────────────────
                          └────┘┌─┴──┐   |    | └────┘┌─┴──┐  |                           
    q_8: ───────────────────────┤CNOT├───┤    ├───────┤CNOT├──■──────────────────────────────────
                                └────┘   └────┘       └────┘
    """

    # initial_state = one(Stabilizer, N)
    g1 = sCNOT(1,4)
    g2 = sCNOT(1,7)

    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    g3= 


end # module

