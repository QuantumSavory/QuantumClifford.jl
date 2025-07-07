"""
    QuantumReedMuller <: AbstractCSSCode
    QuantumReedMuller(m)

The family of `[[2ᵐ - 1, 1, 3]]` CSS Quantum-Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

Quantum codes are constructed from shortened Reed-Muller codes `RM(1, m)`, by removing the first row and column of the generator matrix `Gₘ`. Similarly, we can define truncated dual codes `RM(m - 2, m)` using the generator matrix `Hₘ` [anderson2014fault](@cite). The quantum Reed-Muller codes `QRM(m)` derived from `RM(1, m)` are CSS codes. 

Given that the stabilizers of the quantum code are defined through the generator matrix of the classical code, the minimum distance of the quantum code corresponds to the minimum distance of the dual classical code, which is `d = 3`, thus it can correct any single qubit error. Since one stabilizer from the original and one from the dual code are removed in the truncation process, the code parameters are `[[2ᵐ - 1, 1, 3]]`.

You might be interested in consulting [anderson2014fault](@cite) and [campbell2012magic](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller).

See also: [`ReedMuller`](@ref) and [`RecursiveReedMuller`](@ref)

### Fields
- `m`: The log-length of the code.
"""
struct QuantumReedMuller <: AbstractCSSCode
    m::Int
    function QuantumReedMuller(m)
        if  m < 3
            throw(DomainError("Invalid parameters: m must be bigger than 2 in order to have a valid code."))
        end
        new(m)
    end
end

code_n(c::QuantumReedMuller) = 2^c.m - 1

code_k(c::QuantumReedMuller) = 1

distance(c::QuantumReedMuller) = 3

parity_matrix_x(c::QuantumReedMuller) = generator(RecursiveReedMuller(1, c.m))[2:end, 2:end]

parity_matrix_z(c::QuantumReedMuller) = generator(RecursiveReedMuller(c.m-2, c.m))[2:end, 2:end]