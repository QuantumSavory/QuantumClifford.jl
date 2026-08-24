@testset "Aqua" begin
    using Aqua
    # Persistent-task testing cannot develop packages with local source dependencies.
    Aqua.test_all(QuantumClifford; persistent_tasks=false)
end
