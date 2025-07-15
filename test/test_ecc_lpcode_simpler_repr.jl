@testitem "ECC LPCode representations of commutative algebras" tags=[:ecc] begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC
    using QECCore

    import Hecke
    const QuantumCliffordHeckeExt = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)

    include("test_ecc_base.jl")
    codes = vcat(LP04, LP118, test_gb_codes, test_bb_codes, test_mbb_codes, test_coprimeBB_codes, test_hcubic_codes, test_honeycomb_color_codes);

    for c in codes
        @test QECCore.parity_matrix_xz(c) == QuantumCliffordHeckeExt._parity_matrix_xz_if_comm_algebra(c)
    end

end
