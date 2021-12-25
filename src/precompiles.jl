const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(*),SingleQubitOperator,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.18761368
    Base.precompile(Tuple{typeof(*),sX,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.077348635
    Base.precompile(Tuple{typeof(*),sSWAP,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.05619128
    Base.precompile(Tuple{typeof(*),sCNOT,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.05606607
    let fbody = try __lookup_kwbody__(which(mul_left!, (SubArray{UInt64, 1, Matrix{UInt64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{UInt64, 1, Matrix{UInt64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,typeof(mul_left!),SubArray{UInt64, 1, Matrix{UInt64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{UInt64, 1, Matrix{UInt64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},))
        end
    end   # time: 0.0534866
    Base.precompile(Tuple{typeof(*),sHadamard,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.052111782
    Base.precompile(Tuple{typeof(*),sPhase,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.047713652
    Base.precompile(Tuple{typeof(*),sZ,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.039525334
    Base.precompile(Tuple{typeof(*),sY,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.03883976
    let fbody = try __lookup_kwbody__(which(canonicalize_gott!, (Stabilizer{Vector{UInt8}, Matrix{UInt64}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,typeof(canonicalize_gott!),Stabilizer{Vector{UInt8}, Matrix{UInt64}},))
        end
    end   # time: 0.026974404
    Base.precompile(Tuple{Type{Stabilizer},Matrix{Bool}})   # time: 0.026579907
    Base.precompile(Tuple{typeof(xbit),PauliOperator{Array{UInt8, 0}, Vector{UInt64}}})   # time: 0.025276214
    let fbody = try __lookup_kwbody__(which(gott_standard_form_indices, (SubArray{UInt64, 2, Matrix{UInt64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Int64,Int64,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,typeof(gott_standard_form_indices),SubArray{UInt64, 2, Matrix{UInt64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Int64,Int64,))
        end
    end   # time: 0.023716416
    Base.precompile(Tuple{typeof(zbit),PauliOperator{Array{UInt8, 0}, Vector{UInt64}}})   # time: 0.021749053
    Base.precompile(Tuple{typeof(getindex),Stabilizer{Vector{UInt8}, Matrix{UInt64}},Int64})   # time: 0.005580251
    Base.precompile(Tuple{typeof(getindex),PauliOperator{Array{UInt8, 0}, Vector{UInt64}},Vector{Int64}})   # time: 0.00542647
    Base.precompile(Tuple{typeof(*),sId1,Stabilizer{Vector{UInt8}, Matrix{UInt64}}})   # time: 0.004167271
    Base.precompile(Tuple{typeof(unsafe_bitfindnext_),Vector{UInt64},Int64})   # time: 0.002361816
    Base.precompile(Tuple{typeof(colpermute!),Stabilizer{Vector{UInt8}, Matrix{UInt64}},Vector{Int64}})   # time: 0.002167833
    Base.precompile(Tuple{Core.kwftype(typeof(mul_left!)),NamedTuple{(:phases,), Tuple{Bool}},typeof(mul_left!),Stabilizer{Vector{UInt8}, Matrix{UInt64}},Int64,Int64})   # time: 0.001602084
end

_precompile_()