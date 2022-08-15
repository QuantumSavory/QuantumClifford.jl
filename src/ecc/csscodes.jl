using Nemo 
using LinearAlgebra
#using Statistics

#structure for CSS codes
struct CSS <: AbstractECC end

#=struct CSS <: AbstractECC

    function css_n end
    function classical_code_G_matrix end
    function classical_code_H_matrix end
    function classic_parity_checks end
    function dual_code_prt1 end
    function dual_code_prt2 end

end=#

#Testing code for CSS code contruction--------------------------------------------------
#[7,4] Hamming code -> Steane's 7
#=
#classical_code_H_matrix
classical_code_H_matrix = [0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]
#classical_code_G_matrix - dual code
classical_code_G_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] #tst
#classical_code_G_matrix(c) = gf2_H_to_G(classical_code_H_matrix) 
=#
##---------------------------------------------------------------------------------

##Hamming code generator-------------------------------------------------------------=#
include("hammingcodegenerator.jl") #hamming code generator
import .Hamming
using .Hamming
using .Hamming: construct_A, construct_H_from_A, construct_G_from_A

#=
#function construct_A( block_len::SizeType, message_len::SizeType )
function construct_A()
    #Input desired length - TEMPORARY
    block_length = 12
    message_length = 10

    #Change type input
    block_len = convert(Int64, block_length)
    message_len = convert(Int64, message_length)

    @assert block_len > message_len
    # r is the number of parity bits to be added
    r = block_len - message_len
    # Each parity bit 'n' will be placed in column 2^(n-1)
    # These columns are exised in A, so we need to keep track of them
    parity_cols = [ 2^(x-1) for x in range(1,r) ]
    # Data columns include all columns in A except parity columns
    data_cols = findall( [ !(x in parity_cols) for x in range(1,block_len) ] )
    # Each parity bit 'n' protects the codeword bits in columns where the bit corresponding
    # to 'n' is set. So parity bit 1 protects all odd columns, etc.
    parity_masks = [ (1 << x) for x in range(0,length(parity_cols)) ]
    # First make the whole table, then excise parity columns later
    A_init = zeros( Bit, r, block_len )
    for row = 1:size(A_init,1)
        for column = 1:size(A_init,2)
            # Is this column protected by the parity bit in question?
            A_init[row,column] = ( column & parity_masks[row] ) > 0
        end
    end
    # And now excise the parity columns
    A = zeros( Bit, r, message_len )
    for d = 1:length(data_cols)
        A[:,d] = A_init[:,data_cols[d]]
    end
    return A
end

function construct_G_from_A( A::BitArray )
    # All parameters used to construct A can be retreived from its structure
    r = size(A,1)
    block_len = 2^r - 1
    message_len = 2^r - r - 1
    parity = A
    data = eye( Bit, message_len )
    # Since G is being constructed in its transposed form, we think of parity bits as
    # protecting rows, not columns
    parity_rows = [ 2^(x-1) for x in range(1,r) ]
    # The data rows are all the non-parity rows
    data_rows = find( [ !(x in parity_rows) for x in range(1,block_len) ] )
    G = zeros( Bit, block_len, message_len )
    for p = 1:size(parity,1)
        G[parity_rows[p],:] = parity[p,:]
    end
    for d = 1:size(data,1)
        G[data_rows[d],:] = data[d,:]
    end
    return G
end

    #=
    construct_H_from_A( A::BitArray )
    Create a parity check matrix for a Hamming code of arbitrary size from a
    previously-constructed
    A matrix. It is easiest to pass the return value of construct_A( block_len,
    message_len ) as
    the parameter to this function.
    =#

function construct_H_from_A( A::BitArray )
    # All parameters used to construct A can be retreived from its structure
    r = size(A,1)
    block_len = 2^r - 1
    message_len = 2^r - r - 1
    data = A
    parity = eye( Bit, r )
    # Each parity column 'n' will be placed in column 2^(n-1)
    parity_cols = [ 2^(x-1) for x in range(1,r) ]
    # Data columns are all the non-parity columns
    data_cols = find( [ !(x in parity_cols) for x in range(1,block_len) ] )
    H = zeros( Bit, r, block_len )
    for p = 1:size(parity,2)
        H[:,parity_cols[p]] = parity[:,p]
    end
    for d = 1:size(data,2)
        H[:,data_cols[d]] = data[:,d]
    end
    return H
end
=#
##---------------------------------------------------------------------------------
#=
A = construct_A(7,4)
classical_code_H_matrix = construct_H_from_A(A)
classical_code_G_matrix = construct_G_from_A(A)
=#
#Input desired length - TEMPORARY
block_length = 7
message_length = 4

#Change type input

block_len = convert(Number, block_length)
message_len = convert(Number, message_length)

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

using LinearAlgebra: I
parity = I(r) #WARNING: both LinearAlgebra and QuantumClifford export "I"; uses of it in module ECC must be qualified
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

#----------------------------------------------------------------------------------------------

#----------------CSS code generation ----------------------------------------------------------
#----------Dual code -------------------
#Size
size_row_H = size(H, 1)
size_column_H = size(H, 2)

size_row_G = size(G, 1)
size_column_G = size(G, 2)

#Dual code build
X_zeros = zeros(Int8, size_row_H, size_column_H)
Z_zeros = zeros(Int8, size_row_G, size_column_G)

#Final X & Z matrix
X_matrix = X_zeros
Z_matrix = G
#TODO: overlap X & Z -> Y !!!
hcat(Z_matrix,Z_zeros)
hcat(X_matrix,H)

#-----------Building CSS code ----------------
#RE-SIZE BEFORE NOT AFTER
parity_checks(c::CSS) = Stabilizer(Z_matrix,X_matrix) #READ MANUAL 

code_n(c::CSS) = size(X_matrix, 1) #variable input dependant

parity_matrix(c::CSS) = stab_to_gf2(parity_checks) 

#Encoding circuit ----------------------------------

#encoding_circuit(c::CSS) = [] #TODO -> START SYNDROME CIRCUIT
#----------------------------------------------------------------


code_s(c::CSS) = length(parity_checks(c))

code_k(c::CSS) = code_n(c) - code_s(c)

rate(c::CSS) = code_k(c)/code_s(c)


logx_ops(c::CSS) = P"XXXXXXXXX"

logz_ops(c::CSS) = P"ZZZZZZZZZ"

logy_ops(c::CSS) = P"YYYYYYYYY" 

#naive_syndrome(c::CSS) #Syndrome circuit
