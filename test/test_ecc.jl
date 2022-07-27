

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

#Stabilizer testing -------------------------------------------

function stab_looks_good(s)
    c = canonicalize!(copy(s))
    nrows, ncols = size(c)
    all((c.phases .== 0x0) .| (c.phases .== 0x2)) || return false
    H = stab_to_gf2(c)
    good_indices = reduce(|,H,dims=(1,))
    good_indices = good_indices[1:end÷2] .| good_indices[end÷2+1:end]
    colsok = ncols>nrows || all(good_indices) # TODO, this can be stricter
    colsok || return false
    good_indices = reduce(|,H,dims=(2,))
    rowsok = nrows>ncols || all(good_indices) # TODO, this can be stricter
    rowsok || return false
    check_allrowscommute(c) || return false
    return true
end

#----------------------------------------------------------------