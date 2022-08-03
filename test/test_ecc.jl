using Revise
using BenchmarkTools
using QuantumClifford
using QuantumClifford.Main.ECC: Shor9,code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate
#ENTER+ALT+SHIFT : SQUARE DIV COMMENTS

##Shor9 run ----------------------------
code = Shor9()
parity_checks(code)
encoding_circuit(code)
code_n(code)
#=
code_s(code) #MethodError: no method matching parity_checks(::Shor9)
code_k(code)
rate(code)
=#
#=
distance(code)
logx_ops(code)
logz_ops(code)
isdegenerate(code)
=#
#x = (parity_matrix(code)) #Bool
naive_syndrome_circuit(code)

#Measuring
measured_states = []
some_code_tableau = parity_checks(code)
circuit = encoding_circuit(code)
# I will now measure the syndrome by using the circuit
test_state1 = copy(some_code_tableau)⊗S"X" # add the ancillary qubit
for gate in encoding_circuit(code)
    println(apply!(test_state1,gate))
    #append!(measured_states, apply!(test_state1,gate)) 
end

#measurement_result = project!(measured_states)

# And I will also measure it without constructing a circuit
no_circuit_measurement = project!(copy(some_code_tableau), some_code_tableau[1])

