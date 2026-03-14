# documentation in a bit
struct Mirror <: AbstractQECC
    G::Union{<:Group, <:FinGenAbGroup}
    A::Vector{Tuple}
    B::Vector{Tuple}
    A_elems::Vector{Union{<:GroupElem, <:FinGenAbGroupElem}}
    B_elems::Vector{Union{<:GroupElem, <:FinGenAbGroupElem}}
    symmetric::Bool

    function Mirror(G::Union{<:Group, <:FinGenAbGroup}, A::Vector{T}, B::Vector{T}, symmetric::Bool=true) where T <: Tuple
        isempty(A) && error("A must not empty")
        isempty(B) && error("B must not empty")
        A_elems = [G([Int(a_i) for a_i in a]) for a in A]
        B_elems = [G([Int(b_i) for b_i in b]) for b in B]
        new(G, A, B, A_elems, B_elems, symmetric)
    end
end

function parity_matrix(c::Mirror)
    G, A, B, sym = c.G, c.A_elems, c.B_elems, c.symmetric
    elems = collect(G)
    n = length(elems)
    idx = Dict(elems[i] => i for i in 1:n)
    H = zeros(Int, n, 2n)
    # Symmetric mirror code: S(g) = Z(Ag)X(Bg⁻¹) where Ag = {ag | a ∈ A} and Bg⁻¹ = {bg⁻¹ ∈ B}
    # Asymmetric mirror code: S(g) = Z(Ag)X(g⁻¹B) where Ag = {ag | a ∈ A} and g⁻¹B = {g⁻¹b ∈ B}
    for (i, g) in enumerate(elems)
        if sym
            for b in B
                support_x = b+(-g)
                if haskey(idx, support_x)
                    H[i, idx[support_x]] = 1
                end
            end
        else
            for b in B
                support_x = (-g)+b
                if haskey(idx, support_x)
                    H[i, idx[support_x]] = 1
                end
            end
        end
        for a in A
            support_z = a+g
            if haskey(idx, support_z)
                H[i, n+idx[support_z]] = 1
            end
        end
    end
    return H
end
