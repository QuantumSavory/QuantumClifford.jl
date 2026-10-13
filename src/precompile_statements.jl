# Generated from one successful clean trace of the executable documentation
# examples on Julia 1.12.6 (x86_64-linux-gnu). The trace was filtered with
# sortprecompile.py revision 4ea3c813d9a5adc494d84368e08825f5d765a0a6
# using `-m 50.000001`. Timing comments preserve the largest observation for
# each exact statement within that trace.
# Recompilations and signatures tied to documentation machinery, optional
# dependencies, generated names, compiler internals, or modules not bound in
# QuantumClifford are excluded. QECCore and QuantumInterface names are expressed
# through equivalent bindings already available from QuantumClifford.
# These observations select statements; they are not stable benchmarks.

using PrecompileTools: @setup_workload

@setup_workload begin
    #= 1610.4 ms =# precompile(Tuple{typeof(QuantumClifford.ECC.evaluate_decoder), typeof(QuantumClifford.ECC.TableDecoder(QuantumClifford.ECC.QECCore.Shor9())), QuantumClifford.ECC.CommutationCheckECCSetup, Int64})
    #= 1034.5 ms =# precompile(Tuple{Type{QuantumClifford.ECC.TableDecoder}, QuantumClifford.ECC.QECCore.Shor9})
    #=  463.3 ms =# precompile(Tuple{typeof(QuantumClifford.pftrajectories), QuantumClifford.PauliFrame{QuantumClifford.Stabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, LinearAlgebra.Adjoint{UInt64, Array{UInt64, 2}}}}, Array{Bool, 2}}, Array{QuantumClifford.CompactifiedGate, 1}})
    #=  322.1 ms =# precompile(Tuple{typeof(QuantumClifford.canonicalize_clip!), QuantumClifford.Stabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}})
    #=  316.4 ms =# precompile(Tuple{typeof(Base.alignment), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, QuantumClifford.sCNOT})
    #=  224.4 ms =# precompile(Tuple{typeof(QuantumClifford.petrajectories), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, Array{QuantumClifford.AbstractOperation, 1}})
    #=  214.4 ms =# precompile(Tuple{typeof(QuantumClifford.applywstatus!), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.BellMeasurement})
    #=  174.7 ms =# precompile(Tuple{typeof(QuantumClifford.applywstatus!), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.VerifyOp})
    #=  164.9 ms =# precompile(Tuple{typeof(QuantumClifford.tensor), QuantumClifford.Stabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.Stabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}})
    #=  147.4 ms =# precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:max_order,), Tuple{Int64}}, typeof(QuantumClifford.applybranches), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.BellMeasurement})
    #=  123.9 ms =# precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:keep_result, :phases), Tuple{Bool, Base.Val{true}}}, typeof(QuantumClifford.project_cond!), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, Int64, Base.Val{QuantumClifford.isX}, Base.Val{(false, true)}})
    #=  121.0 ms =# precompile(Tuple{typeof(QuantumClifford.ECC.physical_ECC_circuit), QuantumClifford.Stabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.ECC.ShorSyndromeECCSetup})
    #=  115.4 ms =# precompile(Tuple{typeof(Base.alignment), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, QuantumClifford.sSWAP})
    #=  109.4 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, String, QuantumClifford.PauliMeasurement{Array{UInt8, 0}, Array{UInt64, 1}}})
    #=  105.2 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, String, QuantumClifford.sCNOT})
    #=   97.6 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, String, QuantumClifford.sHadamard})
    #=   96.3 ms =# precompile(Tuple{typeof(QuantumClifford.ECC.naive_encoding_circuit), QuantumClifford.Stabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}})
    #=   80.8 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, String, QuantumClifford.ConditionalGate})
    #=   80.6 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, String, QuantumClifford.NoiseOp{QuantumClifford.UnbiasedUncorrelatedNoise{Float64}, 2}})
    #=   72.8 ms =# precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:trajectories,), Tuple{Int64}}, typeof(QuantumClifford.mctrajectories), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, Array{QuantumClifford.AbstractOperation, 1}})
    #=   67.3 ms =# precompile(Tuple{typeof(Base.alignment), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, QuantumClifford.sCPHASE})
    #=   66.4 ms =# precompile(Tuple{typeof(Base.alignment), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, QuantumClifford.sZCX})
    #=   63.9 ms =# precompile(Tuple{typeof(QuantumClifford.ECC.batchdecode), typeof(QuantumClifford.ECC.TableDecoder(QuantumClifford.ECC.QECCore.Shor9())), Array{Bool, 2}})
    #=   63.7 ms =# precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:max_order,), Tuple{Int64}}, typeof(QuantumClifford.applybranches), QuantumClifford.Register{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.NoiseOp{QuantumClifford.UnbiasedUncorrelatedNoise{Float64}, 7}})
    #=   61.8 ms =# precompile(Tuple{typeof(QuantumClifford.canonicalize_rref!), QuantumClifford.Stabilizer{QuantumClifford.Tableau{Base.SubArray{UInt8, 1, Array{UInt8, 1}, Tuple{Base.UnitRange{Int64}}, true}, Base.SubArray{UInt64, 2, Array{UInt64, 2}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.UnitRange{Int64}}, true}}}, Array{Int64, 1}})
    #=   61.2 ms =# precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:max_order,), Tuple{Int64}}, typeof(QuantumClifford.petrajectories), QuantumClifford.Register{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, Array{QuantumClifford.AbstractOperation, 1}})
    #=   61.0 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, QuantumClifford.sXCX})
    #=   60.9 ms =# precompile(Tuple{typeof(QuantumClifford.canonicalize_rref!), QuantumClifford.Stabilizer{QuantumClifford.Tableau{Base.SubArray{UInt8, 1, Array{UInt8, 1}, Tuple{Base.UnitRange{Int64}}, true}, Base.SubArray{UInt64, 2, Array{UInt64, 2}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.UnitRange{Int64}}, true}}}, Base.UnitRange{Int64}})
    #=   57.3 ms =# precompile(Tuple{typeof(QuantumClifford.random_stabilizer), Int64, Int64})
    #=   52.7 ms =# precompile(Tuple{typeof(Base.show), Base.IOContext{Base.GenericIOBuffer{Memory{UInt8}}}, Array{Union{QuantumClifford.sMX, QuantumClifford.sMY, QuantumClifford.sMZ}, 1}})
    #=   51.6 ms =# precompile(Tuple{typeof(QuantumClifford.applywstatus!), QuantumClifford.MixedDestabilizer{QuantumClifford.Tableau{Array{UInt8, 1}, Array{UInt64, 2}}}, QuantumClifford.sHadamard})
end
