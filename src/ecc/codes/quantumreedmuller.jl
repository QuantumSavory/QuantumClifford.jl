"""
The family of [[2ᵐ - 1, 1, 3]] CSS Quantum-Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

Quantum codes are constructed from shortened Reed-Muller codes RM(1, m), by removing the first row and column of the generator matrix Gₘ. Similarly, we can define truncated dual codes RM(m - 2, m) using the generator matrix Hₘ. The quantum Reed-Muller codes QRM(m) derived from RM(1, m) are CSS codes. This means their stabilizer generators split into two sets, Axₘ and Azₘ. 

The Axₘ elements are formed by converting the 1s in the rows of Gm to Xs and the 0s to Is. The Azₘ elements are similarly formed but using the generator matrix of the truncated dual code RM(m - 2, m). Because RM(1, m) is a subset of RM(m - 2, m), Azₘ includes the same operators as Axₘ with Xs replaced by Zs, along with additional operators from the dual code.

In CSS codes, Ax detects z-type errors and Az detects x-type errors. Given that the stabilizers of the quantum code are defined through the generator matrix of the classical code, the minimum distance of the quantum code corresponds to the minimum distance of the dual classical code, which is d = 3, thus it can correct any single qubit error. Since one stabilizer from the original and one from the dual code are removed in the truncation process, the code parameters are [[2ᵐ - 1, 1, 3]].

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller).
"""
struct QuantumReedMuller <: AbstractECC
    m::Int

    function QuantumReedMuller(m)
        if  m < 3 || m > 11
            throw(ArgumentError("Invalid parameters: m must be ≤ 3 and m ≤ 11 in order to valid code."))
        end
        new(m)
    end
end

function parity_checks(c::QuantumReedMuller)
    RM₁₋ₘ = generator(RecursiveReedMuller(1, c.m))
    RM₍ₘ₋₂₎₋ₘ₎ = generator(RecursiveReedMuller(c.m - 2, c.m))
    QRM = CSS(RM₁₋ₘ[2:end, 2:end], RM₍ₘ₋₂₎₋ₘ₎[2:end, 2:end])
    Stabilizer(QRM)
end

code_n(c::QuantumReedMuller) = 2 ^ c.m - 1

code_k(c::QuantumReedMuller) = 1

distance(c::QuantumReedMuller) = 3
