@testitem "Aqua" tags=[:aqua] begin
    using Aqua
    using QECCore
    Aqua.test_all(QECCore)
end
