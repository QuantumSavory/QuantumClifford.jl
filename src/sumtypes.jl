using SumTypes
using InteractiveUtils: subtypes

"""
```
julia> make_variant(sCNOT)
:(sCNOT(::Int64, ::Int64))
```
"""
function make_variant(type::DataType)
    Expr(:call, Symbol(type), [:(::$t) for t in type.types]...)
end

"""
```
julia> make_variant_deconstruct(sCNOT, :apply!, (:s,))
:(sCNOT(q1, q2) => apply!(s, sCNOT(q1, q2)))
```
"""
function make_variant_deconstruct(type::DataType, call, preargs=(), postargs=())
    variant = Expr(:call, Symbol(type), fieldnames(type)...)
    original = :(($type)($(fieldnames(type)...)))
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
            $([make_variant(t) for t in concrete_types if isa(t, DataType)]...)
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
                $([make_variant_deconstruct(t, call, preargs, postargs) for t in concrete_types if isa(t, DataType)]...)
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
    if isa(type, DataType)
        return :( CompactifiedGate(g::$(type)) = CompactifiedGate'.$(Symbol(type))($([:(g.$n) for n in fieldnames(type)]...)) )
    else
        return :( CompactifiedGate(g::$(type)) = (@warn "The operation is of a type that can not be unified, defaulting to slower runtime dispatch" typeof(g); return g) )
    end
end


function make_all_sumtype_infrastructure_expr(concrete_types, callsigs)
    sumtype = make_sumtype(concrete_types)
    constructors = make_sumtype_variant_constructor.(concrete_types)
    methods = [make_sumtype_method(concrete_types, call, preargs, postargs) for (call, preargs, postargs) in callsigs]
    return quote
        $(sumtype.args...)
        $(constructors...)
        $(methods...)
    end
end

function get_all_concrete_subtypes(type)
    if !isabstracttype(type)
        return [type]
    else
        return vcat(get_all_concrete_subtypes.(subtypes(type))...)
    end
end

module_of_type(t::UnionAll) = module_of_type(t.body)
module_of_type(t::DataType) = t.name.module

function make_all_sumtype_infrastructure_expr(t::DataType, callsigs)
    concrete_types = get_all_concrete_subtypes(t)
    non_experimental_concrete_types = [t for t in concrete_types if module_of_type(t)==QuantumClifford]
    make_all_sumtype_infrastructure_expr(non_experimental_concrete_types, callsigs)
end

function make_all_sumtype_infrastructure()
    make_all_sumtype_infrastructure_expr(AbstractOperation,
        [
            (:apply!, (:(s::Register),), ()),
            (:applywstatus!, (:(s::Register),), ()),
            (:apply!, (:(s::PauliFrame),), ()),
            (:applywstatus!, (:(s::PauliFrame),), ()),
        ]
    ) |> eval
end

make_all_sumtype_infrastructure()

"""
Convert a list of gates to a more optimized "sum type" format which permits faster dispatch.

Generally, this should be called on a circuit before it is used in a simulation.
"""
function compactify_circuit(circuit)
    return CompactifiedGate.(circuit)
end
