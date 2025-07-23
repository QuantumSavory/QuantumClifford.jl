"""The family of `[[2ʲ, 2ʲ - 2j - 2, 4]]` Gottesman codes, also known as a 'Class of Distance Four Codes', as described in [Gottesman's 1997 PhD thesis](@cite gottesman1997stabilizer) and in [gottesman1996class](@cite).

The stabilizer generators of the original `[[2ʲ, 2ʲ - j - 2, 3]]` Gottesman codes are incorporated to create a set of generators for this distance-four code. The resulting stabilizer set, denoted by S, incorporates the following elements: The first two generators are the Pauli-`X` and Pauli-`Z` operators acting on all qubits, represented by `Mₓ` and `Mz`, respectively. The next `j` generators correspond to `M₁` through `Mⱼ`, which are directly inherited from the `[[2ʲ, 2ʲ - j - 2, 3]]` Gottesman  code's stabilizers. This inclusion ensures that `S` retains the inherent distance-three property of the original Gottesman code. The final `j` generators are defined as `Nᵢ = RMᵢR`, where `i` ranges from `1` to `j`. Here, `R` signifies a Hadamard Rotation operation applied to all `2ʲ` qubits, and `Mᵢ` refers to one of the existing generators from the second set `(M₁ to Mⱼ)`. By incorporating the stabilizers of a distance-three code, the constructed set `S` inherently guarantees a minimum distance of three for the resulting distance-four Gottesman code.
"""
struct Gottesman4 <: AbstractECC
    j::Int
    function Gottesman4(j)
        (j >= 3 && j < 21) || error("In `Gottesman4(j)`, `j` must be ≥  3 in order to obtain a valid code and < 21 to remain tractable")
        new(j)
    end
end

code_n(c::Gottesman4) = 2^c.j
code_k(c::Gottesman4) = 2^c.j - 2*c.j - 2
distance(c::Gottesman4) = 4

function parity_checks(c::Gottesman4)
    H₁ = parity_checks(Gottesman(c.j))
    Hⱼ = H₁[3:end]
    for qᵢ in 1:nqubits(Hⱼ)
        apply!(Hⱼ, sHadamard(qᵢ))
    end
    H = vcat(H₁, Hⱼ)
    Stabilizer(H)
end
