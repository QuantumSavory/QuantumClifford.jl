using Test
using QuantumClifford
using QuantumClifford.ECC: AbstractECC, Cleve8, Steane7, Shor9, Bitflip3, Perfect5,
    naive_syndrome_circuit, naive_encoding_circuit, code_n, parity_checks, code_s, code_k

##

@testset "encoding circuits - compare to algebraic construction of encoded state" begin
    # This test verifies that logical measurements on an encoded state match the physical pre-encoded state.
    # This test skips verifying the permutations of qubits during canonicalization are properly undone,
    # i.e. we modify the code we are testing so that the canonicalization does not need any permutations.
    for codeexpr in [
        :(Cleve8()),
        :(Steane7()),
        :(Shor9()),
        :(Perfect5()),
        :(Bitflip3()),
        :(S"Y_"),
        :(S"Z_"),
        :(S"X_"),
        fill(:(random_stabilizer(5,7)), 100)...
        ]

        # pre-process the tableau to remove permutations and negative phases
        code = eval(codeexpr)
        stab, r, s, xperm, zperm = canonicalize_gott!(parity_checks(code))
        code = stab # using this tableau guarantees we do not need to worry about permutations of the qubits

        circ = naive_encoding_circuit(code)
        #display(circ)

        # the state to be encoded (k physical qubits)
        pre_encₖ = one(Stabilizer, code_k(code))

        # n-k ancillary qubits in state zero prepended
        pre_encₙ = one(Stabilizer, code_s(code)) ⊗ pre_encₖ

        # running the encoding circuit
        encodedₙ = mctrajectory!(pre_encₙ, circ)[1] |> canonicalize!

        # creating a valid state purely algebraically
        algebraicₙ = MixedDestabilizer(code)
        algebraicₙ.rank = nqubits(algebraicₙ)
        algebraicₙ = stabilizerview(algebraicₙ) |> canonicalize!

        @test (encodedₙ == algebraicₙ)

        @show affectedqubits.(circ)
        @show circ
        #println("$codeexpr, $(encodedₙ == algebraicₙ)")
    end
end
