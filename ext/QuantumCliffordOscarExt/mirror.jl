# documentation in a bit
struct Mirror{T<:Union{Tuple, GroupElem}, S<:Union{Tuple, GroupElem}} <: AbstractQECC
    G::Union{<:Group, <:FinGenAbGroup}
    A::Vector{T}
    B::Vector{S}
    symmetric::Bool

    function Mirror(G::Union{<:Group, <:FinGenAbGroup}, A::Vector{T}, B::Vector{S}, symmetric::Bool=true) where {T<:Union{Tuple, GroupElem}, S<:Union{Tuple, GroupElem}}
        isempty(A) && error("A must not empty")
        isempty(B) && error("B must not empty")
        new{T,S}(G, A, B, symmetric)
    end
end

function parity_matrix(c::Mirror)
    G, A_vec, B_vec, sym = c.G, c.A, c.B, c.symmetric
    check_G = is_abelian(G)
    A = [a isa Tuple ? G([Int(a_i) for a_i in a]) : a for a in A_vec]
    B = [b isa Tuple ? G([Int(b_i) for b_i in b]) : b for b in B_vec]
    elems = collect(G)
    n = length(elems)
    idx = Dict(elems[i] => i for i in 1:n)
    H = zeros(Int, n, 2n)
    # Symmetric mirror code: S(g) = Z(Ag)X(Bg⁻¹) where Ag = {ag | a ∈ A} and Bg⁻¹ = {bg⁻¹ ∈ B}
    # Asymmetric mirror code: S(g) = Z(Ag)X(g⁻¹B) where Ag = {ag | a ∈ A} and g⁻¹B = {g⁻¹b ∈ B}
    for (i, g) in enumerate(elems)
        if sym
            for b in B
                support_x = check_G ? b+(-g) : b*inv(g)
                haskey(idx, support_x) && (H[i, idx[support_x]] = 1)
            end
        else
            for b in B
                support_x = check_G ? (-g)+b : inv(g)*b
                haskey(idx, support_x) && (H[i, idx[support_x]] = 1)
            end
        end
        for a in A
            support_z = check_G ? a+g : a*g
            haskey(idx, support_z) && (H[i, n+idx[support_z]] = 1)
        end
    end
    H
end
