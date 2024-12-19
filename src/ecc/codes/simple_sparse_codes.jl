"""A bicycle (LDPC) code based on [mackay2004sparse](@cite). This code is produced by taking a circulant matrix based on a difference set, reducing the number of rows, and then using it as the base matrix for a CSS code.

Parameters:
$TYPEDFIELDS

```jldoctest Bicycle
julia> Bicycle(4, 2).n
4

julia> Bicycle(4, 2).m
2

julia> parity_checks(Bicycle(4, 2))
+ XXXX
+ ZZZZ

julia> parity_checks(Bicycle(6, 4))
+ XX_X_X
+ X_X_XX
+ ZZ_Z_Z
+ Z_Z_ZZ
```
"""
struct Bicycle
    "Width of array, should be >= 2. Equal to number to number of physical qubits (N)"
    n::Int
    "Height of array, should be >= 2 and a multiple of 2. Equal to number of parity checks (M)"
    m::Int
end

function parity_checks(b::Bicycle)
    m = b.m
    n = b.n
    if m%2 == 1
        throw(DomainError(m, " M should be a multiple for 2 for bicycle codes."))
    end
    if m < 2
        throw(DomainError(m, " M is too small, make it greater than 1."))
    end
    if n%2 == 1
        throw(DomainError(n, " N should be a multiple for 2 for bicycle codes."))
    end
    if n < 2
        throw(DomainError(n, " N is too small, make it greater than 1."))
    end
    bs = bicycle_set_gen(Int(n/2))
    bsc = circ_to_bicycle_h0(bs, Int(n/2))
    while size(bsc)[1] > m/2
        bsc = reduce_bicycle(bsc)
    end
    return parity_checks(CSS(bsc, bsc))
end

"""A unicycle (LDPC) code based on [mackay2004sparse](@cite). The parity check matrix is produced by taking a perfect difference set, turning it into a circulant matrix, and removing the linearly dependent rows.

Has the drawback of only being possible to make at sizes that correspond to existing perfect difference sets. Reference [mackay2004sparse](@cite) to find good sizes of perfect difference sets to use.

Parameters:
$TYPEDFIELDS

```jldoctest Unicycle
julia> typeof(Unicycle(21, [1, 3, 8, 9, 12]))
Unicycle

julia> parity_checks(Unicycle(21, [1, 3, 8, 9, 12]))
+ _X_________X_X____XX_X
+ __X_________X_X____XXX
+ X__X_________X_X____XX
+ XX__X_________X_X____X
 ⋮
+ ____ZZ__Z_________Z_ZZ
+ Z____ZZ__Z_________Z_Z
+ _Z____ZZ__Z_________ZZ
```
"""
struct Unicycle
    "Size of parity check matrix, also max size of generated set array, should be >= 1. Equal to the number of physical qubits."
    N::Int
    "Perfect difference set modulo N. Array of indices that are 'active' checks in the circulant code."
    set::Vector{Int}
end

function parity_checks(u::Unicycle)
    n = u.N
    set = u.set
    usc = circ_to_unicycle_h0(set, n)
    rusc = reduce_unicycle(usc) # reduced unicycle code
    return parity_checks(CSS(rusc, rusc))
end

"""Takes an untrimmed bicycle matrix and removes the row which keeps the spread of the column weights minimal.

Required before the bicycle code can be used.

```jldoctest reduce_bicycle

julia> reduce_bicycle(circ_to_bicycle_h0([1, 2, 4], 7))
6×14 Matrix{Bool}:
 1  1  0  1  0  0  0  1  0  0  0  1  0  1
 0  1  1  0  1  0  0  1  1  0  0  0  1  0
 0  0  0  1  1  0  1  1  0  1  1  0  0  0
 1  0  0  0  1  1  0  0  1  0  1  1  0  0
 0  1  0  0  0  1  1  0  0  1  0  1  1  0
 1  0  1  0  0  0  1  0  0  0  1  0  1  1

```
"""
function reduce_bicycle(H0::Matrix{Bool})
    m, n = size(H0)
    r_i = 0
    std_min = Inf
    for i in 1:m
        t_H0 = vcat(H0[1:i-1, :], H0[i+1:end, :])
        std_temp = std(convert(Array, sum(t_H0, dims = 1)))
        if std_temp < std_min
            std_min = std_temp
            r_i = i
        end
    end
    return vcat(H0[1:r_i-1, :], H0[r_i+1:end, :])
end

"""Takes a list of indices and creates the base of the bicycle matrix.

For example:

```jldoctest circ_to_bicycle_h0

julia> circ_to_bicycle_h0([1, 2, 4], 7)
7×14 Matrix{Bool}:
 1  1  0  1  0  0  0  1  0  0  0  1  0  1
 0  1  1  0  1  0  0  1  1  0  0  0  1  0
 0  0  1  1  0  1  0  0  1  1  0  0  0  1
 0  0  0  1  1  0  1  1  0  1  1  0  0  0
 1  0  0  0  1  1  0  0  1  0  1  1  0  0
 0  1  0  0  0  1  1  0  0  1  0  1  1  0
 1  0  1  0  0  0  1  0  0  0  1  0  1  1
```

See [mackay2004sparse](@cite) for more details
"""
function circ_to_bicycle_h0(circ_indices, n::Int)
    circ_arr = Vector{Bool}(undef, n)
    circ_matrix = Matrix{Bool}(undef, n, n)
    comp_matrix = Matrix{Bool}(undef, n, 2*n)
    for i = 1:n
        if Int(i) in circ_indices
            circ_arr[i] = true
        else
            circ_arr[i] = false
        end
    end
    for i = 1:n
        circ_matrix[i,1:n] = circ_arr
        li = circ_arr[end]
        circ_arr[2:end] = circ_arr[1:end-1]
        circ_arr[1] = li
    end
    comp_matrix[1:n,1:n] = circ_matrix
    comp_matrix[1:n,n+1:2*n] = transpose(circ_matrix)
    return comp_matrix
end

"""Takes an untrimmed unicycle matrix and removes linearly dependent rows.

Required before the unicycle code can be used.

Typical usage:

```jldoctest reduce_unicycle
julia> reduce_unicycle(circ_to_unicycle_h0([1, 2, 4], 7))
4×8 Matrix{Bool}:
 0  0  0  1  1  0  1  1
 1  0  0  0  1  1  0  1
 0  1  0  0  0  1  1  1
 1  0  1  0  0  0  1  1

```
"""
function reduce_unicycle(m::Matrix{Bool})
    rrzz, = residue_ring(ZZ, 2)
    nm = matrix(rrzz, m)
    r = LinearAlgebra.rank(nm)
    i = 1
    while size(m)[1] > r
        tm = vcat(m[1:i-1,:], m[i+1:end,:])
        tr = LinearAlgebra.rank(matrix(rrzz, tm))
        if(tr == r)
            m = tm
            if(size(m)[1] == r)
                break
            end
        else
            i += 1
        end
    end
    return m
end

"""Takes a list of indices and creates the base of the unicycle matrix.

For example:

```jldoctest circ_to_unicycle_h0
julia> circ_to_unicycle_h0([1, 2, 4], 7)
7×8 Matrix{Bool}:
 1  1  0  1  0  0  0  1
 0  1  1  0  1  0  0  1
 0  0  1  1  0  1  0  1
 0  0  0  1  1  0  1  1
 1  0  0  0  1  1  0  1
 0  1  0  0  0  1  1  1
 1  0  1  0  0  0  1  1

```
See [mackay2004sparse](@cite) for more details
"""
function circ_to_unicycle_h0(circ_indices::Vector{Int}, n::Int)
    circ_arr = fill(false, n)
    one_col = transpose(fill(true, n))
    circ_matrix = Matrix{Bool}(undef, n, n)
    comp_matrix = Matrix{Bool}(undef, n, n+1)
    for i = 1:n
        if i in circ_indices
            circ_arr[i] = true
        else
            circ_arr[i] = false
        end
    end
    for i = 1:n
        circ_matrix[i,1:n] = circ_arr
        li = circ_arr[end]
        circ_arr[2:end] = circ_arr[1:end-1]
        circ_arr[1] = li
    end
    comp_matrix[1:n,1:n] = circ_matrix
    comp_matrix[1:n,n+1] = one_col
    return comp_matrix
end

"""
Generates a list of indices to be used in a bicycle code using a search method given a number of qubits (N). This algirithm finds indicies which have the property that no 2 differences between elements occurs more than once. The differences can also loop over the end of the array. For example, the method will never return the indices 1, 3, and 5, as the difference '2' occurs twice. Occurrences are: 3-1 = 2 and 5-3 = 2.

Note: This algorithm can, but is not guaranteed, to find the optimal sets of indices for certain numbers of qubits. For example N = 13 will give you a perfect difference set, the optimal set of indices. We know this is the optimal set as there are no ways to increase the number of unique differences because the number of unique distances has been maximized per definition of a perfect difference set. Optimal indices can be found via a Monte Carlo search or brute force method. Generally, one does not need to worry about small numbers of missing differences as the function often returns the optimal set, but it is hard to prove so for larger qubit numbers.

This algorithm find the indices by iterating a value 'i' over 1 to N and doing the following for each index: The algorithm adds 'i' to the set of indices, checks if it violates the difference set property, keeps it in the set if it does not violate the afformentioned property, and otherwise removes.

```jldoctest bicycle_set_gen
julia> bicycle_set_gen(13)
4-element Vector{Int64}:
  1
  2
  4
 10

```
"""
function bicycle_set_gen(N::Int)
    circ_arr = Int[0]
    diff_arr = Int[]
    # test new elements
    for add_i = (circ_arr[end] + 1):N - 1
        valid = true
        temp_circ_arr = copy(circ_arr)
        temp_diff_arr = []
        push!(temp_circ_arr, add_i)
        for j = 1:size(temp_circ_arr)[1]
            temp_arr = copy(temp_circ_arr)
            # add lesser elements + N to temp_arr
            for k = 1:size(temp_circ_arr)[1]
                if k < j
                    push!(temp_arr, temp_circ_arr[k] + N)
                else
                    break
                end
            end
            # test if new index is valid
            for k = 1:(size(temp_circ_arr)[1] - 2)
                t_diff = (temp_arr[j + k] - temp_arr[j]) % N
                if ((t_diff) in temp_diff_arr)
                    valid = false
                    break
                else
                    push!(temp_diff_arr, t_diff)
                end
            end
            if !valid
                break
            end
        end
        if valid
            circ_arr = copy(temp_circ_arr)
            diff_arr = copy(temp_diff_arr)
        end
    end
    return circ_arr .+ 1
end
