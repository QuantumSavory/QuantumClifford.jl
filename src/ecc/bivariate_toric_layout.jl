using LinearAlgebra
using SparseArrays: SparseMatrixCSC, sparse, spzeros

function cyclic_shift_matrix(l::Int)
    arr = spzeros(Int,l,l)
    for i in 1:l
        arr[i, mod(i,l)+1] = 1
    end  
    return arr
end

pow²(mat::SparseMatrixCSC{Int, Int},exp::Int) = (mat^exp).%2

function order(mat::SparseMatrixCSC{Int}, m::Int, l::Int)
    n = size(mat,1)
    I = sparse(LinearAlgebra.I(n))
    p = mat
    for i in 1:m*l
        if p == I return i end
        p *= mat
    end
    return -1
end

function _toric_layout₁(As::Vector{SparseMatrixCSC{Int64, Int64}}, Bs::Vector{SparseMatrixCSC{Int64, Int64}}, m::Int, l::Int)
    A_ord = [order(AA,m,l) for AA in As]
    B_ord = [order(BB,m,l) for BB in Bs]
    valid_orders = []
    for (i, Ao) in enumerate(A_ord)
        for (j, Bo) in enumerate(B_ord)
            Ao * Bo == m*l && push!(valid_orders, (Ao,Bo,i,j))
        end
    end
    return valid_orders
end

function _toric_layout₂(As::Vector{SparseMatrixCSC{Int64, Int64}}, Bs::Vector{SparseMatrixCSC{Int64, Int64}}, m::Int, l::Int, codes::NTuple{4,Int64})
    mₑ, lₑ, A_idx, B_idx = codes
    (A_idx < 1 || A_idx > length(As)) && throw(DimensionMismatch("Invalid index A_idx: $A_idx"))
    (B_idx < 1 || B_idx > length(Bs)) && throw(DimensionMismatch("Invalid index B_idx: $B_idx"))
    visited = Set{Int}()                            
    v = sparse(As[A_idx])                        
    h = sparse(Bs[B_idx])                      
    zero = zeros(Int, m*l)                    
    zero[1] = 1
    for i in 0:mₑ-1
        tmp = (v^i)*zero
        for j in 0:lₑ-1
            visited = union(visited, Set(findall(j == 0 ? tmp .>0 : h^j*tmp .> 0)))
        end
    end
    return length(visited) == l*m
end

"""Determines if a bivariate bicycle code has a toric layout on a `2ℓ × 2m` grid."""
function bivariate_toric_layout(code::Vector{Int})
    l = code[1]
    m = code[2]
    x = kron(cyclic_shift_matrix(l), LinearAlgebra.I(m))
    y = kron(LinearAlgebra.I(l), cyclic_shift_matrix(m))
    A₁ = pow²(x, code[3])
    A₂ = pow²(y, code[4])
    A₃ = pow²(y, code[5])
    A = (A₁+A₂+A₃).%2
    B₁ = pow²(y, code[6])
    B₂ = pow²(x, code[7])
    B₃ = pow²(x, code[8])
    B = (B₁+B₂+B₃).%2
    As = [A₁*A₂', A₂'*A₁, A₂*A₃', A₃'*A₂, A₁*A₃', A₃'*A₁]
    Bs = [B₁*B₂', B₂'*B₁, B₂*B₃', B₃'*B₂, B₁*B₃', B₃'*B₁]
    selected_codes = []
    valid_codes = _toric_layout₁(As, Bs, m, l)
    for valid_code in valid_codes
        _toric_layout₂(As, Bs, m, l, valid_code) && push!(selected_codes, valid_code)
    end
    orig_layout = [(params[1], params[2], params[3]-1, params[4]-1) for params in selected_codes]
    return orig_layout
end
