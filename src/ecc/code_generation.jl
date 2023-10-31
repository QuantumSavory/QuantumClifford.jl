# using QuantumClifford: Stabilizer
# using QuantumClifford.ECC: AbstractECC
# import QuantumClifford.ECC: parity_checks
# using Statistics:std
# using Nemo: residue_ring, matrix
# using LinearAlgebra: rank

"""Struct for arbitrary CSS error correcting codes.

This struct holds:
    - tab: Boolean matrix with the X part taking up the left side and the Z part taking up the right side
    - stab: Stabilizer of the parity check matrix
    - n: Block length
    - d: Code distance"""
struct CSS <: AbstractECC
    tab::Matrix{Bool}
    stab::Stabilizer
    n::Int
end

"""Takes an untrimmed bicycle matrix and removes the row which keeps the spread of the column weights minimal.

Required before the bicycle code can be used.

Typical usage:
ReduceBicycle(Circ2BicycleH0(array_indices, (block length / 2) ) )"""
function ReduceBicycle(H0::Matrix{Bool})
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
Circ2BicycleH0([1, 2, 4], 7)

See https://arxiv.org/abs/quant-ph/0304161 for more details"""
function Circ2BicycleH0(circ_indices::Array{Int}, n::Int)
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
ReduceUnicycle(Circ2UnicycleH0(array_indices, block length) )"""
function ReduceUnicycle(m::Matrix{Bool})
    r = LinearAlgebra.rank(nm7)
    rrzz = Nemo.residue_ring(Nemo.ZZ, 2)
    for i in 1:size(u7)[1]
        tm = vcat(m[1:i-1,:], m[i+1:end,:])
        tr = LinearAlgebra.rank(Nemo.matrix(rrzz, tm))
        if(tr == r)
            m = tm
            i -= 1
            if(size(m)[1] == r)
                break
            end
        end
    end
    return m
end

"""Takes a list of indices and creates the base of the unicycle matrix.

For example:
Circ2UnicycleH0([1, 2, 4], 7)

See https://arxiv.org/abs/quant-ph/0304161 for more details"""
function Circ2UnicycleH0(circ_indices::Array{Int}, n::Int)
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

function AssembleCSS end

"""Creates a CSS code using the two provided matrices where H contains the X checks and G contains the Z checks."""
function AssembleCSS(H::Matrix{Bool}, G::Matrix{Bool})::CSS
    Hy, Hx = size(H)
    Gy, Gx = size(G)
    comp_matrix = fill(false, (Hy + Gy, Hx + Gx))
    # comp_matrix = Matrix{Bool}(undef, Hy + Gy, Hx + Gx)
    comp_matrix[1:Hy, 1:Hx] = H
    comp_matrix[Hy+1:end, Hx+1:end] = G
    pcm_stab = Stabilizer(fill(0x0, Hy+Gy), GetXTableau(comp_matrix), GetZTableau(comp_matrix))
    return CSS(comp_matrix, pcm_stab, Hx)
    # return comp_matrix
end

"""Creates a CSS code using the provided matrix for the X and Z checks."""
function AssembleCSS(H::Matrix{Bool})::CSS
    return AssembleCSS(H, H)
end

"""Attempts to generate a list of indices to be used in a bicycle code using a search method"""
function BicycleSetGen(N::Int)
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
function BicycleSetGenRand(N::Int, d::Int)
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

"""Takes in a boolean Matrix and returns the parity check tableau as a string of characters.

Note: Only works when the block length for the X and Z checks are the same!"""
function GetCodeTableau(ecc::Matrix{Bool})
    eccx = size(ecc)[2]
    eccy = size(ecc)[1]
    ps::String = ""
    for i = 1:size(ecc)[1]
        for j = 1:(Int(size(ecc)[2]/2))
            if (ecc[i, j] == 0) && (ecc[i, j + Int(eccx / 2)] == 0)
                ps = string(ps, "I")
            elseif (ecc[i, j] == 1) && (ecc[i, j + Int(eccx / 2)] == 0)
                ps = string(ps, "X")
            elseif (ecc[i, j] == 1) && (ecc[i, j + Int(eccx / 2)] == 1)
                ps = string(ps, "Y")
            else
                ps = string(ps, "Z")
           end
        end
        ps = string(ps,"\n")
    end
    return ps
end

"""Takes in a matrix and returns just the X checks portion while keeping the full height of the matrix.

Note: Only works when the block length for the X and Z checks are the same!"""
function GetXTableau(ecc::Matrix{Bool})
   return ecc[1:size(ecc)[1], 1:Int(size(ecc)[2]/2)]
end

"""Takes in a matrix and returns just the Z checks portion while keeping the full height of the matrix.

Note: Only works when the block length for the X and Z checks are the same!"""
function GetZTableau(ecc::Matrix{Bool})
   return ecc[1:size(ecc)[1], Int(size(ecc)[2]/2) + 1:end]
end

"""Takes in a matrix and returns just the X checks portion while keeping the full height of the matrix.

Note: Only works when the block length for the X and Z checks are the same!"""
function GetXTableau(ecc::CSS)
   return GetXTableau(ecc.tab)
end

"""Takes in a matrix and returns just the Z checks portion while keeping the full height of the matrix.

Note: Only works when the block length for the X and Z checks are the same!"""
function GetZTableau(ecc::CSS)
   return GetZTableau(ecc.tab)
end

"""Returns the matrix form of the X and Z checks."""
tableau(c::CSS) = c.tab

"""Returns the stabilizer making up the parity check tableau."""
parity_checks(c::CSS) = c.stab

"""Returns the block length of the code."""
code_n(c::CSS) = c.n