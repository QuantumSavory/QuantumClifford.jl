using QuantumClifford
using QuantumClifford: AbstractOperation

@testset "Syndrome Measurements with mctrajectory!" begin # TODO this is a rather old test that is now done in a few other places, e.g. the ECC module -- double check and consider deleting
    codeˢᵗᵉᵃⁿᵉ = S"Z__Z_ZZ
                   _Z_ZZ_Z
                   __Z_ZZZ
                   X__X_XX
                   _X_XX_X
                   __X_XXX";
    function get_check(code,index) # TODO a fairly ugly piece of code that should just be part of the library
        s,n = size(code)
        @assert index<=s
        op = code[index]
        anc = n+1
        circuit = AbstractOperation[]
        measurement = nothing
        resetgate = nothing
        for qubit in 1:n
            if op[qubit] == (false, true) # it is a Z
                if isnothing(measurement) || measurement==sMZ
                    resetgate = Reset(S"Z",[anc])
                    measurement = sMZ
                else
                    error("can not do mixed syndromes")
                end
                push!(circuit, sCNOT(qubit, anc))
            elseif op[qubit] == (true, false) # it is a X
                if isnothing(measurement) || measurement==sMX
                    resetgate = Reset(S"X",[anc])
                    measurement = sMX
                else
                    error("can not do mixed syndromes")
                end
                push!(circuit, sCNOT(anc, qubit))
            elseif op[qubit] == (true, true) # it is a Y
                error("can not do Y syndromes")
            end
        end
        push!(circuit,measurement(anc,index))
        circuit = AbstractOperation[resetgate,circuit...]
    end
    function get_all_checks(code)
        s,n = size(code)
        vcat([get_check(code,i) for i in 1:s]...)
    end
    # test the code is correct
    state = Register(MixedDestabilizer(codeˢᵗᵉᵃⁿᵉ ⊗ S"Z"),zeros(Bool,6))
    full_circuit = get_all_checks(codeˢᵗᵉᵃⁿᵉ)
    for _ in 1:100
        err = random_pauli(7;realphase=true)
        @test mctrajectory!(apply!(copy(state),(err ⊗ P"I")),full_circuit)[1].bits == comm(err*codeˢᵗᵉᵃⁿᵉ, err)
    end
end
