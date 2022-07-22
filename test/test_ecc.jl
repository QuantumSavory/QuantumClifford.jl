#Testing file ecc library




#Syndrome circuit testing-------------------------------------

circuit = naive_syndrome_circuit(s,i) #TO DO: other syndrome circuits
test_state = random_stabilizer(n)
results_direct = project!(copy(test_state), s[i]) # directly 
circuit_state = test_state⊗S"Z" # state on n+1 qubits

for gate in circuit; apply!(circuit_state,gate); end;
    results_circuit = projectZ!(circuit_state,n+1) # measure the ancilla qubit

#--------------------------------------------------------------

#Enconding circuit testing-------------------------------------



#--------------------------------------------------------------


