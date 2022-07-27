#https://dspace.sunyconnect.suny.edu/bitstream/handle/1951/70206/R%20Quinn%20MS%20Project%20document%2020180507.pdf?sequence=1&isAllowed=y

module Hamming

#----------------Generating Hamming codes ----------------------------------------------------
#= TODO:
    Number the bits starting from 1: bit 1, 2, 3, 4, 5, 6, 7, etc.
    Write the bit numbers in binary: 1, 10, 11, 100, 101, 110, 111, etc.
    All bit positions that are powers of two (have a single 1 bit in the binary form of their position) are parity bits: 1, 2, 4, 8, etc. (1, 10, 100, 1000)
    All other bit positions, with two or more 1 bits in the binary form of their position, are data bits.
    Each data bit is included in a unique set of 2 or more parity bits, as determined by the binary form of its bit position.
        Parity bit 1 covers all bit positions which have the least significant bit set: bit 1 (the parity bit itself), 3, 5, 7, 9, etc.
        Parity bit 2 covers all bit positions which have the second least significant bit set: bits 2-3, 6-7, 10-11, etc.
        Parity bit 4 covers all bit positions which have the third least significant bit set: bits 4–7, 12–15, 20–23, etc.
        Parity bit 8 covers all bit positions which have the fourth least significant bit set: bits 8–15, 24–31, 40–47, etc.
        In general each parity bit covers all bits where the bitwise AND of the parity position and the bit position is non-zero.
=#

# INCLUDES
include("types.jl")

# EXPORTS
export Hamming_Params
export construct_A, construct_G_from_A, construct_H_from_A,
construct_R_from_A, construct_hamming_parameters
export encode, decode, encode_stream, decode_stream

# TYPES
type Hamming_Params
A::BitArray
G::BitArray
H::BitArray
R::BitArray
end

# FUNCTIONS

"""
construct_A( block_len::SizeType, message_len::SizeType )
Create a valid A matrix for a Hamming code of arbitrary size. It is the
responsibility of
the caller to ensure that block_len and message_len are valid parameters for
a Hamming code.
Hamming codes are generally named \"Hamming-(block_len, message_len).\"
Common Hamming codes include:
* Hamming-(7,4)
* Hamming-(15,11)
* Hamming-(31,26)
* etc.
"""
function construct_A( block_len::SizeType, message_len::SizeType )
    @assert block_len > message_len
    # r is the number of parity bits to be added
    r = block_len - message_len
    # Each parity bit 'n' will be placed in column 2^(n-1)
    # These columns are exised in A, so we need to keep track of them
    parity_cols = [ 2^(x-1) for x in range(1,r) ]
    # Data columns include all columns in A except parity columns
    data_cols = find( [ !(x in parity_cols) for x in range(1,block_len) ] )
    # Each parity bit 'n' protects the codeword bits in columns where the bit
    corresponding
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

"""
construct_G_from_A( A::BitArray )
Create a genrator matrix for a Hamming code of arbitrary size from a
previously-construct A
matrix. It is easiest to pass the return value of construct_A( block_len,
message_len ) as
the parameter to this function.
"""
function construct_G_from_A( A::BitArray )
    # All parameters used to construct A can be retreived from its structure
    r = size(A,1)
    block_len = 2^r - 1
    message_len = 2^r - r - 1
    parity = A
    data = eye( Bit, message_len )
    # Since G is being constructed in its transposed form, we think of parity
    bits as
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

"""
construct_H_from_A( A::BitArray )
Create a parity check matrix for a Hamming code of arbitrary size from a
previously-constructed
A matrix. It is easiest to pass the return value of construct_A( block_len,
message_len ) as
the parameter to this function.
"""
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


"""
construct_hamming_parameters( block_len::SizeType, message_len::SizeType
)
Construct the A, G, H, and R matrices necessary for Hamming encoding and
decoding.
"""
function construct_hamming_parameters( block_len::SizeType,
    message_len::SizeType )
    A = construct_A( block_len, message_len )
    G = construct_G_from_A( A )
    H = construct_H_from_A( A )
    return Hamming_Params( A, G, H, R )
end