import .ECC

using LinearAlgebra: I
using .ECC
#=
This code was published as an Appendix in the following project by Ryan Quinn:
Applicability of the Julia Programming Language to Forward Error-Correction Coding in Digital Communications Systems
https://dspace.sunyconnect.suny.edu/bitstream/handle/1951/70206/R%20Quinn%20MS%20Project%20document%2020180507.pdf?sequence=1&isAllowed=y
This project was presented to the Department of Computer and Information Sciences, at State University of New York Polytechnic Institute, Utica, noisy
In partial fulfillment of the requirements of the author's Master of Science Degree, in May 2018
=#

struct Codedata <: AbstractECC 
    block_length
    message_length
end

function Hamming(c::Codedata)

    block_len = c.block_length
    message_len = c.message_length
    @assert block_len > message_len

    r = block_len - message_len
    
    parity_cols = [ 2^(x-1) for x in range(1,r) ]
    
    data_cols = findall( [ !(x in parity_cols) for x in range(1,block_len) ] )
    
    parity_masks = [ (1 << x) for x in range(0,length(parity_cols)) ]
    
    # First make the whole table, then excise parity columns later
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
    
    r = size(A,1)
    block_len = 2^r - 1
    message_len = 2^r - r - 1
    data = A
    
    parity = I(r) #WARNING: both LinearAlgebra and QuantumClifford export "I"; uses of it in module ECC must be qualified
    
    # Each parity column 'n' will be placed in column 2^(n-1)
    parity_cols = [ 2^(x-1) for x in range(1,r) ]
    
    # Data columns are all the non-parity columns
    data_cols = findall( [ !(x in parity_cols) for x in range(1,block_len) ] )
    H = zeros( r, block_len )
    
    for p = 1:size(parity,2)
        H[:,parity_cols[p]] = parity[:,p]
    end
    
    for d = 1:size(data,2)
        H[:,data_cols[d]] = data[:,d]
    end

    #--------------

    r = size(A,1)
    block_len = 2^r - 1
    message_len = 2^r - r - 1
    parity = A
    data = I(message_len)
    
    # Since G is being constructed in its transposed form, we think of parity bits as
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
    
    r = size(A,1)
    block_len = 2^r - 1
    message_len = 2^r - r - 1
    parity = A
    data = I(message_len)
    
    # Since G is being constructed in its transposed form, we think of parity bits as
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

    return H,G

end #function
