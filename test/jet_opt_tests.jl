@testitem "JET optimization: projective measurements" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    state = one(MixedDestabilizer, 2)

    JET.@test_opt target_modules=(QuantumClifford,) projectY!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectY!(copy(state), 1; phases=false)
    JET.@test_opt target_modules=(QuantumClifford,) projectZ!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectZ!(copy(state), 1; phases=false)
end
