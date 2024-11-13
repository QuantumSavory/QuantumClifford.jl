"""Takes a height and width of matrix and generates a bicycle code to the specified height and width.
Based on codes from [mackay_sparse-graph_2004](@cite)

Parameters:
$TYPEDFIELDS

```jldoctest Bicycle
julia> O = typeof(Bicycle(4, 2))
Bicycle

julia> O = parity_checks(Bicycle(4, 2))
+ XXXX
+ ZZZZ

julia> O = parity_checks(Bicycle(6, 4))
+ XX_X_X
+ X_X_XX
+ ZZ_Z_Z
+ Z_Z_ZZ
```
"""
struct Bicycle
    "width of array, should be >= 2"
    n::Int
    "height of array, should be >= 2 and a multiple of 2"
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
        throw(DomainError(m, " N should be a multiple for 2 for bicycle codes."))
    end
    if n < 2
        throw(DomainError(m, " N is too small, make it greater than 1."))
    end
    bs = bicycle_set_gen(Int(n/2))
    bsc = circ_to_bicycle_h0(bs, Int(n/2))
    while size(bsc)[1] > m/2
        print(bsc)
        bsc = reduce_bicycle(bsc)
    end
    return parity_checks(CSS(bsc, bsc))
end

"""Takes a height and width of matrix and generates a unicycle code to the specified height and width.
Based on codes from [mackay_sparse-graph_2004](@cite)

Parameters:
$TYPEDFIELDS 

Reference [mackay_sparse-graph_2004](@cite) to find good sizes of perfect difference sets to use.

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
    "size of parity check matrix, also max size of generated set array, should be >= 1."
    N::Int
    "array of indices that are 'active' checks in the circulant code"
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

Typical usage:
```jldoctest reduce_bicycle

julia> reduce_bicycle(circ_to_bicycle_h0([1, 2, 4], 7))
6×14 Matrix{Bool}:
 0  0  1  1  0  1  0  1  0  0  0  1  0  1
 0  0  0  1  1  0  1  1  1  0  0  0  1  0
 1  0  0  0  1  1  0  0  1  1  0  0  0  1
 0  1  0  0  0  1  1  1  0  1  1  0  0  0
 1  0  1  0  0  0  1  0  1  0  1  1  0  0
 1  1  0  1  0  0  0  0  0  1  0  1  1  0

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
 0  1  1  0  1  0  0  0  0  0  1  0  1  1
 0  0  1  1  0  1  0  1  0  0  0  1  0  1
 0  0  0  1  1  0  1  1  1  0  0  0  1  0
 1  0  0  0  1  1  0  0  1  1  0  0  0  1
 0  1  0  0  0  1  1  1  0  1  1  0  0  0
 1  0  1  0  0  0  1  0  1  0  1  1  0  0
 1  1  0  1  0  0  0  0  0  1  0  1  1  0

```

See [mackay_sparse-graph_2004](@cite) for more details"""
function circ_to_bicycle_h0(circ_indices, n::Int)
    circ_arr = Array{Bool}(undef, n)
    circ_matrix = Matrix{Bool}(undef, n, n)
    comp_matrix = Matrix{Bool}(undef, n, 2*n)
    for i = 1:n
        if Int(i-1) in circ_indices
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
See [mackay_sparse-graph_2004](@cite) for more details"""
function circ_to_unicycle_h0(circ_indices::Array{Int}, n::Int)
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

"""Attempts to generate a list of indices to be used in a bicycle code using a search method"""
function bicycle_set_gen(N::Int)
    circ_arr::Array{Int} = [0]
    diff_arr::Array{Int} = []
    circ_arr[1] = 0
    # test new elements
    for add_i = (circ_arr[end] + 1):N - 1
        valid = true
        temp_circ_arr = copy(circ_arr)
        temp_diff_arr::Array{Int} = []
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
    return circ_arr
end

"""Attempts to generate a list of indices to be used in a bicycle code using a randomized check method

Note: This is very slow for large N"""
function bicycle_set_gen_rand(N::Int, d::Int)
    circ_arr::Array{Int} = [0]
    diff_arr::Array{Int} = []
    atmp_add::Array{Int} = [0]
    circ_arr[1] = 0
    # test new elements
    for i = (circ_arr[end] + 1):(N^2)
        valid = true
        temp_circ_arr = copy(circ_arr)
        temp_diff_arr::Array{Int} = []
        add_i = rand(1: N-1)
        atmp_add = push!(atmp_add, add_i)
        if add_i in circ_arr
            continue
        end
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
            if (size(atmp_add)[1] == N) || (size(circ_arr)[1] == d)
                break
            end
        end
    end
    return circ_arr
end
