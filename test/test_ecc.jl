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


##---------------------------------------------------------------
##---------------------------------------------------------------
#TESTING HAMMING CODE GENERATORS
##---------------------------------------------------------------
##---------------------------------------------------------------

using LinearAlgebra

#Input desired length - TEMPORARY
block_length = 7
message_length = 4

#Change type input
#=
block_len = convert(Number, block_length)
message_len = convert(Number, message_length)
Bit = convert(Number, Bits)
=#
@assert block_len > message_len

r = block_len - message_len

parity_cols = [ 2^(x-1) for x in range(1,r) ]

data_cols = findall( [ !(x in parity_cols) for x in range(1,block_len) ] )

parity_masks = [ (1 << x) for x in range(0,length(parity_cols)) ]

# First make the whole table, then excise parity columns later
#A_init = zeros( Bit, r, block_len ) #SOME SORT OF ERROR HERE - TEST
A_init = zeros( r, block_len )

for row = 1:size(A_init,1)
    for column = 1:size(A_init,2)
        # Is this column protected by the parity bit in question?
        A_init[row,column] = ( column & parity_masks[row] ) > 0
    end
end

# And now excise the parity columns
A = zeros( r, message_len )
for d = 1:length(data_cols)
    A[:,d] = A_init[:,data_cols[d]]
end

println("A:",A) #empty entry
r = size(A,1)
block_len = 2^r - 1
message_len = 2^r - r - 1
data = A
println("data",data)

parity = I(r)
println("parity",parity)

# Each parity column 'n' will be placed in column 2^(n-1)
parity_cols = [ 2^(x-1) for x in range(1,r) ]

# Data columns are all the non-parity columns
data_cols = findall( [ !(x in parity_cols) for x in range(1,block_len) ] )
H = zeros( r, block_len )

for p = 1:size(parity,2)
    H[:,parity_cols[p]] = parity[:,p]
    println("p:",p)
end

for d = 1:size(data,2)
    H[:,data_cols[d]] = data[:,d]
    println("d:",d)
end
#println(parity)

println(data)

#println("This is H: \n",H)

##------------
using LinearAlgebra

r = size(A,1)
block_len = 2^r - 1
message_len = 2^r - r - 1
parity = A
data = I(message_len)

# Since G is being constructed in its transposed form, we think of parity bits as
# protecting rows, not columns
parity_rows = [ 2^(x-1) for x in range(1,r) ]

# The data rows are all the non-parity rows
data_rows = findall( [ !(x in parity_rows) for x in range(1,block_len) ] )
G = zeros( block_len, message_len )

for p = 1:size(parity,1)
    G[parity_rows[p],:] = parity[p,:]
end

for d = 1:size(data,1)
    G[data_rows[d],:] = data[d,:]
end

println(G)

##-----------------------
using LinearAlgebra

I(4)