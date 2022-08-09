using Revise 
using BenchmarkTools
using QuantumClifford
using QuantumClifford.ECC: Shor9, Bitflip3, Steane7, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate
#ENTER+ALT+SHIFT : SQUARE DIV COMMENTS

#Shor9 run ----------------------------
code = Shor9()

parity_checks(code)

#code_s(code) #MethodError: no method matching parity_checks(::Shor9)
encoding_circuit(code)
code_n(code)

typeof(parity_checks(code))
#println(parity_checks(code))

code_s(code) #MethodError: no method matching parity_checks(::Shor9)
#println("Code s:", code_s(code))
code_k(code)
#println(code_s(code))
#println(code_k(code))
rate(code)

distance(code)
logx_ops(code)
logz_ops(code)
isdegenerate(code)

#println("Length parity checks:", length(parity_checks(code)))
println("Naive syndrome Bit flip:\n")
code2 = Bitflip3()
println(naive_syndrome_circuit(code2))
println(length(parity_checks(code2)))


#x = (parity_matrix(code)) #Bool
#naive_syndrome_circuit(code)

#Measuring
measured_states = []
some_code_tableau = parity_checks(code)
circuit = encoding_circuit(code)
println("Parity checks: \n",some_code_tableau)
##---------------------------------------------------------------

println("Rows in parity checks:" ,length(parity_checks(code)), "\n")
code2 = Bitflip3()
println("Bit flip rows in parity checks:" , code_n(code2))

naive_syndrome_circuit(code)

##---------------------------------------------------------------
# I will now measure the syndrome by using the circuit
test_state1 = copy(some_code_tableau)⊗S"X" # add the ancillary qubit
counter = 1
for gate in encoding_circuit(code)
    println("Gate ", counter, ": \n")
    println(apply!(test_state1,gate))
    println("\n")
    counter = counter+1 #should there be 11 gates??
    #append!(measured_states, apply!(test_state1,gate)) 
end

#measurement_result = project!(measured_states)

##---------------------------------------------------------------
# And I will also measure it without constructing a circuit
#no_circuit_measurement = project!(copy(some_code_tableau), some_code_tableau[1])
