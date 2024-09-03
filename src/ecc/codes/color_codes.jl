abstract type ColorCode <: AbstractECC end
abstract type TriangularCode <: ColorCode end

struct Triangular4_8_8 <: TriangularCode
    d::Int
    function Triangular4_8_8(d)
        if d%2!=1
            throw(DomainError("only odd distance trianglular color codes are allowed.\nRefer to https://arxiv.org/abs/1108.5738"))
        end
        return new(d)
    end
end

Triangular4_8_8() = Triangular4_8_8(3) # smallest d

function parity_checks(c::Triangular4_8_8)
    # Note that this is the same as the Steane7 code
    if c.d == 3
        return S"___XXXX
                _XX__XX
                X_X_X_X
                ___ZZZZ
                _ZZ__ZZ
                Z_Z_Z_Z"
    else
        throw("Parity checks for a distance $(c.d) code of type $(typeof(c)) are not currently implemented.")
    end
end

function iscss(::ColorCode)
    return true
end

"""From https://arxiv.org/abs/1108.5738 Fig. 2's caption:"""
code_n(c::Triangular4_8_8) = 0.5*c.d^2+c.d-0.5 |> Int