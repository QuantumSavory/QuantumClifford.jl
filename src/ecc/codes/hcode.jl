"""
The family of `[[k + 4, k, 2]]` CSS quantum codes that encode an even number `k` logical qubits using `(k + 4)` physical qubits, collectively referred to as 'H codes' as first described by Cody Jones in his 2013 paper [jones2013multilevel](@cite). Each H code is denoted as `Hₙ`, where `n = k + 4` corresponds to the number of physical qubits. All H codes have a distance of two, indicating they can detect a single physical Pauli error.

The construction method of Stabilizer works as follows:
- The first stabilizer generator  S₁ is a product of Pauli-X operators on the first four physical qubits: `S₁ = X₁X₂X₃X₄`.
- The second stabilizer generator S₂ is a product of Pauli-Z operators on the first four physical qubits: `S₂ = Z₁Z₂Z₃Z₄`.
- The third stabilizer generator  S₃ includes Pauli-X operators on the first two qubits and from the fifth to the `n`-th qubit: `S₃ = X₁X₂X₅X₆...Xₙ`.
- The fourth stabilizer generator S₄ includes Pauli-Z operators on the first two qubits and from the fifth to the `n`-th qubit: `S₄ = Z₁Z₂Z₅Z₆...Zₙ`.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_h).
"""
struct HCode <: AbstractECC
    k::Int
    function HCode(k)
        (k >= 2 && iseven(k)) || error("In `HCode(k)`, `k` must be ≥ 2 and even in order to obtain a valid code.")
        new(k)
    end
end

code_n(c::HCode) = c.k + 4
code_k(c::HCode) = c.k
distance(c::HCode) = 2
 
function parity_checks(c::HCode)
    Hx = zeros(Int64, 2, c.k + 4) 
    Hz = zeros(Int64, 2, c.k + 4)
    Hx[1, 1:4] .= 1 # S₁ = X₁X₂X₃X₄
    Hz[1, 1:4] .= 1 # S₂ = Z₁Z₂Z₃Z₄
    Hx[2, 1:2] .= 1 # S₃ = X₁X₂ ...
    Hz[2, 1:2] .= 1 # S₄ = Z₁Z₂ ...
    for c in 5:c.k + 3
        Hx[2, c:c+1] .= 1 # S₃ = ...X₅X₆...Xₙ
        Hz[2, c:c+1] .= 1 # S₄ = ...Z₅Z₆...Zₙ
    end
    return Stabilizer(CSS(Hx, Hz))
end
