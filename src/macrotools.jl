using MacroTools

"""Turns `f(Val(x))` into `x ? f(Val(true)) : f(Val(false))` in order to avoid dynamic dispatch

See [discourse discussion](https://discourse.julialang.org/t/allocations-due-to-boolean-keyword-arguments-how-to-avoid-them/87654)"""
macro valbooldispatch(expr, bools...)
    for bool in bools
        true_branch  = MacroTools.postwalk(x->@capture(x, Val($bool)) ? :(Val(true )) : x, expr)
        false_branch = MacroTools.postwalk(x->@capture(x, Val($bool)) ? :(Val(false)) : x, expr)
        expr = quote
            if $bool
                $true_branch
            else
                $false_branch
            end
        end
    end
    esc(expr)
end
