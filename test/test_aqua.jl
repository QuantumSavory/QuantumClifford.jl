@testitem "Aqua" tags=[:aqua] begin
    using Aqua
    Aqua.test_all(QuantumClifford)
end
root_project_path, found = Aqua.root_project_toml(Base.PkgId(Base.UUID("0525e862-1e90-11e9-3e4d-1b39d7109de1"), "QuantumClifford"))
Aqua.precompile_wrapper(root_project_path,10,quote end)