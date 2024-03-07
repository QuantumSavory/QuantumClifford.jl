using Test
using QuantumClifford
using QuantumClifford.ECC

##

@testset "encoding circuits - compare to algebraic construction of encoded state" begin
    # This test verifies that logical measurements on an encoded state match the physical pre-encoded state.
    # This test skips verifying the permutations of qubits during canonicalization are properly undone,
    # i.e. we modify the code we are testing so that the canonicalization does not need any permutations.
    for undoperm in [true, false],
        codeexpr in [
        :(Cleve8()),
        :(Steane7()),
        :(Shor9()),
        :(Perfect5()),
        :(N8_K3_D3()),
        :(N16_K10_D3()),
        :(Bitflip3()),
        :(S"Y_"),
        :(S"Z_"),
        :(S"X_"),
        :(CSS([0 1 1 0; 1 1 0 0], [1 1 1 1])),
        :(Toric(3,3)),
        :(Toric(3,6)),
        :(Toric(6,4)),
        :(Toric(8,8)),
        fill(:(random_stabilizer(5,7)), 100)...
        ]

        code = eval(codeexpr)
        if undoperm==false
            # Pre-process the tableau to remove permutations and negative phases.
            # Usually that is handled by `naive_encoding_circuit`, but we just want to check both branches for its `undoperm` kwarg.
            stab, r, s, xperm, zperm = canonicalize_gott!(parity_checks(code))
            code = stab # using this tableau guarantees we do not need to worry about permutations of the qubits
        end

        circ = naive_encoding_circuit(code; undoperm=true)
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

        #println("$codeexpr, $(encodedₙ == algebraicₙ)")
    end
end
