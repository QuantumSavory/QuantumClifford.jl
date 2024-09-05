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
    matrix = nothing
    if c.d == 3 # Note that this is the same as the Steane7 code
        matrix = Matrix{Bool}([
            [1 1 1 1 0 0 0];
            [1 0 1 0 1 1 0]; 
            [0 0 1 1 0 1 1]])
    elseif c.d== 5
        matrix = Matrix{Bool}([
            [1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
            [1 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0];
            [0 0 1 1 0 1 1 0 0 1 1 0 0 1 1 0 0];
            [0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0];
            [0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0];
            [0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 1];
            [0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0];
            [0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0]])
    end

    if isnothing(matrix)
        throw("Parity checks for a distance $(c.d) code of type $(typeof(c)) are not currently implemented.")
    else
        num_checks, qubits = size(matrix)
        return Stabilizer(vcat(matrix,zeros(Bool,num_checks,qubits)), vcat(zeros(Bool,num_checks,qubits), matrix))
    end
end

function prototype_parity_checks(c::Triangular4_8_8)
    n = code_n(c)
    num_checks = (n-1)/2 |> Int

    checks = zeros(Bool, num_checks, n)

    num_layers = (c.d-1)/2 |> Int

    i = 1
    checks_written = 0
    for layer in 1:num_layers
        # TODO convert half 8-gons from previous layer into full 8-gons

        # red 4-gons
        for j in 0:(layer-1)
            checks[checks_written+1,i+j*2] = 1
            checks[checks_written+1,i+1+j*2] = 1
            checks[checks_written+1,i+2*layer+j*2] = 1
            checks[checks_written+1,i+2*layer+1+j*2] = 1
            checks_written += 1
        end

        if layer%2 == 1
            # green half 8-gons on left side
            checks[checks_written+1,i] = 1
            checks[checks_written+1,i+2*layer] = 1
            checks[checks_written+1,i+4*layer] = 1
            checks[checks_written+1,i+4*layer+1] = 1
            checks_written += 1

            
        else
            # blue half 8-gon on right side
            checks[checks_written+1,i+1+(layer-1)*2] = 1
            checks[checks_written+1,i+2*layer+1+(layer-1)*2] = 1

            # when d=5,9,13,... the final row qubits is indexed slightly differently.
            offset = 0
            if layer == num_layers
                offset = 1
            end

            checks[checks_written+1,(i+4*layer)+layer*2-offset] = 1
            checks[checks_written+1,(i+4*layer)+1+layer*2-offset] = 1
            checks_written += 1
        end

        # blue/green half 8-gons on the bottom
        for j in 0:(layer-1)
            checks[checks_written+1,i+2*layer+j*2] = 1
            checks[checks_written+1,i+2*layer+1+j*2] = 1

            offset = 0
            if layer == num_layers && num_layers%2==0
                offset = 1
            end

            checks[checks_written+1,(2*layer+1)+i+2*layer+j*2 - offset] = 1
            checks[checks_written+1,(2*layer+1)+i+2*layer+1+j*2 - offset] = 1

            checks_written += 1
        end

        i += 4*layer  
    end

    return checks
end

function iscss(::ColorCode)
    return true
end

"""From https://arxiv.org/abs/1108.5738 Fig. 2's caption:"""
code_n(c::Triangular4_8_8) = 0.5*c.d^2+c.d-0.5 |> Int
code_k(c::TriangularCode) = 1