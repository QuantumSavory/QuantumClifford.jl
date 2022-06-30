<<<<<<< HEAD

export Shorcircuit

function Shorcircuit()#Codes::Shorcode)
    println("I work") #test

    #=
    N= 9 #n qubits 
           #Step 1   #Step 2 #Step 3     #Step 4 #S 5  #S 6  #Step 7  #Step 8 #Step 9  #S 10 #S 11
    
=======
#currently in the works
using QuantumClifford import X, Z

module Codes end
abstract type Code end

export Shorcode, Shorcircuit

struct Shorcode <: Code end

"""documents"""

function rate end

function logicalqubits end

function physicalqubits end

function codedistance end

"""end documents"""

function rate(code::Shor_code) return 1//9 end

function Shorcircuit(code::Shor_code)
    N= 9 #n qubits 

           #Step 1   #Step 2 #Step 3     #Step 4 #S 5  #S 6  #Step 7  #Step 8 #Step 9  #S 10 #S 11
    """
>>>>>>> 4921213 (error library made aware of existing codes)
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
<<<<<<< HEAD
    =#

    #Step 0
    #initial_state = one(Stabilizer, N) #CHECK THIS
=======
    """

    #Step 0
<<<<<<< HEAD
    # initial_state = one(Stabilizer, N)
>>>>>>> 4921213 (error library made aware of existing codes)
=======
    initial_state = one(Stabilizer, N) #CHECK THIS
>>>>>>> 322da12 (comment)

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
<<<<<<< HEAD
    c6 = sCNOT(7,9)
    
    #Step 4: Error
    single_x(9,1) #X: Bit flip error
    single_z(9,1) #Z: Phase flip error
    #check x and z have the right def (a,b)
=======
    c6 = sCNOT(7,9)   
    
    #Step 4: Error
<<<<<<< HEAD
    #X: Bit flip error
    #Z: Phase flip error
    for qubit in range(N) #check this
        #X = S"X"
        #Z = S"Z"
        X(qubit)
        Z(qubit)
>>>>>>> 4921213 (error library made aware of existing codes)
=======
    single_x(9,1) #X: Bit flip error
    single_z(9,1) #Z: Phase flip error
    #check x and z have the right def (a,b)
>>>>>>> f69433d (bit and phase flip error def)

    #Step 5: 4th set of  CNOT GATES
    c7 = sCNOT(1,2)
    c8 = sCNOT(4,5)
    c9 = sCNOT(7,8)

    #Step 6: 5th set of  CNOT gates
    c10 = sCNOT(1,2)
    c11 = sCNOT(4,5)
    c12 = sCNOT(7,8)
<<<<<<< HEAD

    #Step 7: 1st set of  Toffoli gates
    #are toffoli gates represented by sCCNOT?
    #MethodError
    cc11 = sCNOT(2,1)
    cc12 = sCNOT(3,1)
    cc21 = sCNOT(5,4)
    cc22 = sCNOT(6,4)
    cc31 = sCNOT(8,7)
    cc32 = sCNOT(9,7)
=======
    
    #Step 7: 1st set of  Toffoli gates
    #are toffoli gates represented by sCCNOT?
    cc1 = sCCNOT(2,3,1)
    cc2 = sCCNOT(5,6,4)
    cc3 = sCCNOT(9,8,7)
>>>>>>> 4921213 (error library made aware of existing codes)

    #Step 8: 2nd set of  Haramard gates
    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    #Step 9: 6th set of  CNOT gates
    c13 = sCNOT(1,4)

    #Step 10: 7th set of CNOT gates
    c14 = sCNOT(1,7)

<<<<<<< HEAD

    #Step 11: 2nd set of Toffoli gates
    #Final gates
    cc4 = sCNOT(4,1)
    cc4 = sCNOT(7,1)

    # This circuit performs a depolarization at rate `epsilon` to all qubits,
    #circuit = [c2,c2,h1,h2,h3,c3,c4,c5,c6,single_x,single_z,c7,c8,c9,c10,c11,c12,cc11,cc12,cc21,cc22,cc31,cc32,h1,h2,h3,c13,c14,cc41,cc42]

end #Shorcircuit
=======
    #Step 11: 2nd set of Toffoli gates
    #Final gates
    cc4 = sCCNOT(7,4,1)

<<<<<<< HEAD
# This circuit performs a depolarization at rate `epsilon` to all qubits,
circuit = [c2,c2,h1,h2,h3,c3,c4,c5,c6,X,Z,c7,c8,c9,c10,c11,c12,cc1,cc2,cc3,h1,h2,h3,c13,c14,cc4]
<<<<<<< HEAD

end # module
>>>>>>> 4921213 (error library made aware of existing codes)
=======
>>>>>>> 0d82741 (testing error line 3 expected "end" got "code")
=======
    # This circuit performs a depolarization at rate `epsilon` to all qubits,
    circuit = [c2,c2,h1,h2,h3,c3,c4,c5,c6,single_x,single_z,c7,c8,c9,c10,c11,c12,cc1,cc2,cc3,h1,h2,h3,c13,c14,cc4]

    end #end Shorcircuit

<<<<<<< HEAD
end #module
>>>>>>> df1b849 (reformating shor mod, struct & function hyerarchy)
=======
end #function
>>>>>>> 322da12 (comment)
