# Here be dragons...

using SumTypes
using InteractiveUtils: subtypes

"""An intermediary when we want to create a new concrete type in a macro."""
struct SymbolicDataType
    name::Symbol
    types#::Core.SimpleVector
    fieldnames
    originaltype
    originaltype_parameterized
end
_header(s) = s
_header(s::SymbolicDataType) = s.name
_symbol(s) = Symbol(s)
_symbol(s::SymbolicDataType) = s.name
_types(s) = s.types
_fieldnames(s) = fieldnames(s)
_fieldnames(s::SymbolicDataType) = s.fieldnames
_originaltype(s) = s
_originaltype(s::SymbolicDataType) = s.originaltype
_originaltype_parameterized(s) = s
_originaltype_parameterized(s::SymbolicDataType) = s.originaltype_parameterized

"""
```
julia> make_variant(sCNOT)
:(sCNOT(::Int64, ::Int64))
```
"""
function make_variant(type::Union{DataType,SymbolicDataType})
    Expr(:call, _symbol(type), [:(::$t) for t in _types(type)]...)
end

"""
```
julia> make_variant_deconstruct(sCNOT, :apply!, (:s,))
:(sCNOT(q1, q2) => apply!(s, sCNOT(q1, q2)))
```
"""
function make_variant_deconstruct(type::Union{DataType,SymbolicDataType}, call, preargs=(), postargs=())
    variant = Expr(:call, _symbol(type), _fieldnames(type)...)
    original = :(($(_originaltype_parameterized(type)))($(_fieldnames(type)...)))
    #:($variant => begin $(Expr(:call, call, preargs..., original, postargs...)); nothing end) # useful when you are searching for type instabilities due to inconsistent output types for a method (usually also pointing to a method not following the conventions of the API)
    :($variant => $(Expr(:call, call, preargs..., original, postargs...)))
end

"""
```
julia> make_sumtype([sCNOT])
quote
    @sum_type CompactifiedGate :hidden begin
        sCNOT(::Int64, ::Int64)
    end
end
```
"""
function make_sumtype(concrete_types)
    return quote
        @sum_type CompactifiedGate :hidden begin
            $([make_variant(t) for t in concrete_types if isa(t, DataType) || isa(t, SymbolicDataType)]...)
        end
    end
end

"""
```
julia> make_sumtype_method([sCNOT], :apply!, (:s,))
quote
    function QuantumClifford.apply!(s, g::CompactifiedGate)
        @cases g begin
            sCNOT(q1, q2) => apply!(s, sCNOT(q1, q2))
        end
    end
end
"""
function make_sumtype_method(concrete_types, call, preargs=(), postargs=())
    return quote
        function QuantumClifford.$call($(preargs...), g::CompactifiedGate, $(postargs...))
            @cases g begin
                $([make_variant_deconstruct(t, call, preargs, postargs) for t in concrete_types if isa(t, DataType) || isa(t, SymbolicDataType)]...)
                #_ => @error "something wrong is happening when working with $(g) -- you are probably getting wrong results, please report this as a bug" # this being present ruins some safety guarantees, but it is useful for debugging
            end
        end
    end
end

"""
```
julia> make_sumtype_variant_constructor(sCNOT)
:(CompactifiedGate(g::sCNOT) = begin
    (CompactifiedGate').sCNOT(g.q1, g.q2)
end)
```
"""
function make_sumtype_variant_constructor(type)
    if isa(type, DataType) || isa(type, SymbolicDataType)
        return :( CompactifiedGate(g::$(_header(type))) = CompactifiedGate'.$(_symbol(type))($([:(g.$n) for n in _fieldnames(type)]...)) )
    else
        return :() # this is taken care of by a default constructor that also warns about the failure to compactify
    end
end

genericsupertypeparams(t) = :body ∈ propertynames(t) ? genericsupertypeparams(t.body) : t

"""Returns a tuple of all concrete subtypes and all UnionAll non-abstract subtypes of a given type."""
function get_all_subtypes(type)
    if !isabstracttype(type)
        if isa(type, DataType)
            isbitstype(type) || @debug "$type will be problematic during compactification"
            return [type], []
        elseif isa(type, UnionAll)
            return [], [type]
        else
            @error "The gate compiler has encountered a type that it can not handle: $type. The QuantumClifford library should continue functioning, but potentially at degraded performance. Please report this as a performance bug."
        end
    else
        return Iterators.flatten.(zip(get_all_subtypes.(subtypes(type))...))
    end
end

module_of_type(t::UnionAll) = genericsupertypeparams(t).name.module
module_of_type(t::DataType) = t.name.module

function make_all_sumtype_infrastructure_expr(t::DataType, callsigs)
    concrete_types, unionall_types = get_all_subtypes(t)
    concrete_types = collect(Any, concrete_types)
    concrete_types = Any[t for t in concrete_types if module_of_type(t)==QuantumClifford]
    unionall_types = Any[t for t in unionall_types if module_of_type(t)==QuantumClifford]
    concretifier_workarounds_types = [] # e.g. var"ClassicalXOR_{2}" generated as a workaround for providing a concrete type for ClassicalXOR{N}
    concretifier_additional_constructors = [] # e.g. CompactifiedGate(g::ClassicalXOR{2}) = CompactifiedGate'.var"ClassicalXOR_{2}"(g.bits, g.store)
    for ut in unionall_types
        names, generated_concretetypes, generated_variant_constructors = concretifier(ut)
        append!(concretifier_workarounds_types, generated_concretetypes)
        append!(concrete_types, names)
        append!(concretifier_additional_constructors, generated_variant_constructors)
        push!(concrete_types, ut) # fallback
    end
    sumtype = make_sumtype(concrete_types)
    @debug "compiling a total of $(length(concrete_types)) concrete types"
    constructors = make_sumtype_variant_constructor.(concrete_types)
    methods = [make_sumtype_method(concrete_types, call, preargs, postargs) for (call, preargs, postargs) in callsigs]
    modulename = gensym(:CompactifiedGate)
    return quote
        #module $(modulename)
        #using QuantumClifford
        #import QuantumClifford: CompactifiedGate, # todo
        $(concretifier_workarounds_types...)
        $(sumtype.args...) # defining the sum type
        $(constructors...) # creating constructors for the sumtype which turn our concrete types into instance of the sum type
        $(concretifier_additional_constructors...) # creating constructors for the newly generated "workaround" concrete types
        :(CompactifiedGate(g::AbstractOperation) = (@warn "The operation is of a type that can not be unified, defaulting to slower runtime dispatch" typeof(g); return g) )
        $(methods...)
        #end
    end
end

function concrete_typeparams(t)
    @debug "The gate compiler is not able to concretify the type $t. Define a `concrete_typeparams` method for this type to improve performance."
    return ()
end

function concretifier(t)
    names = []
    generated_concretetypes = []
    generated_variant_constructors = []
    for typeparams in concrete_typeparams(t)
        name = Symbol(t,"{",typeparams,"}")
        parameterized_type = t{typeparams...}
        ftypes = parameterized_type.types
        fnames = fieldnames(t)
        push!(names, SymbolicDataType(name, ftypes, fnames, t, parameterized_type))
        push!(generated_concretetypes, :(
            struct $(name)
                $([:($n::$t) for (n,t) in zip(fnames,ftypes)]...)
            end
        ))
        push!(generated_variant_constructors, make_concretifier_sumtype_variant_constructor(parameterized_type, name))
    end
    return names, generated_concretetypes, generated_variant_constructors
end

function make_concretifier_sumtype_variant_constructor(parameterized_type, variant_name)
    return :( CompactifiedGate(g::$(parameterized_type)) = CompactifiedGate'.$(variant_name)($([:(g.$n) for n in _fieldnames(parameterized_type)]...)) )
end

function make_all_sumtype_infrastructure()
    make_all_sumtype_infrastructure_expr(AbstractOperation,
        [
            (:apply!, (:(s::Register),), ()),
            (:applywstatus!, (:(s::Register),), ()),
            (:apply!, (:(s::PauliFrame),), ()),
            (:applywstatus!, (:(s::PauliFrame),), ()),
            (:affectedqubits, (), ()),
            (:affectedbits, (), ()),
        ]
    ) |> eval
end

"""
Convert a list of gates to a more optimized "sum type" format which permits faster dispatch.

Generally, this should be called on a circuit before it is used in a simulation.
"""
function compactify_circuit(circuit)
    return CompactifiedGate.(circuit)
end


##
# `concrete_typeparams` annotations for the parameteric types we care about
##

function concrete_typeparams(t::Type{ClassicalXOR})
    return 2:16
end

function concrete_typeparams(t::Type{NoiseOp})
    return [
        [(UnbiasedUncorrelatedNoise{Float64}, i) for i in 1:8];
        [(PauliNoise{Float64}, i) for i in 1:8];
    ]
end


# XXX This has to happen after defining all the `concrete_typeparams` methods

make_all_sumtype_infrastructure()
