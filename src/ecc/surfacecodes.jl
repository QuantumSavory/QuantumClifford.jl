import .ECC

using LinearAlgebra
using .ECC

struct Surface <: AbstractECC 
    lx
    ly
end 

function calculating_n(c::Surface)
    if c.lx == c.ly
        n = 2*c.lx*(c.lx+1)
    elseif c.lx > c.ly
        n = (2*c.lx*(c.lx+1)) - ((c.lx-c.ly)+1)
    else
        n = (2*c.ly*(c.ly+1)) - ((c.ly-c.lx)+1)
    end

    return n
end

code_n(c::Surface) = calculating_n(c)

#Parity checks ----------------------------------

function plaquette_to_qubit_indices(rown,columnn) 
    
    q1 = [rown,columnn]
    q2 = [rown+1,columnn]
    q3 = [rown+1,columnn+1]
    q4 = [rown+2,columnn]
    return q1,q2,q3,q4

end


function grid_index_to_linear_index(c,q)::Int8 #to test with big grids
    a,b = q 
    @show(q)
    if a<= 2
        l = c.lx*(a-1) +b 
    elseif a%2 == 0
        l = c.lx*(a-1) +b +a/2 -1
    else
        l = c.lx*(a-1) +b +(a-1)/2 
    end
    return l
    @show(l)
end

function checks_s(c::Surface) 
    rown = c.lx*2 + 1
    columnn = c.ly +1
    Z = zeros(code_n(c),code_n(c))
    X = zeros(code_n(c),code_n(c))
    stab_index = 1

    for i1 in range(1,rown)
        for i2 in range(1,columnn)
            if (i1+2 <= (rown)) && (i2+1< columnn)
                q1, q2, q3, q4 = plaquette_to_qubit_indices(i1,i2)
                # q1 is a tuple
                Q1 = grid_index_to_linear_index(c,q1) # Q1 is an Int
                Q2 = grid_index_to_linear_index(c,q2)
                @show(Q2)
                Q3 = grid_index_to_linear_index(c,q3)
                Q4 = grid_index_to_linear_index(c,q4)
                
                if ((i1%2 != 0) && (i2%2 != 0)) || ((i1%2 == 0) && (i2%2 == 0))
                    Z[stab_index, Q1] = 1 
                    Z[stab_index, Q2] = 1 
                    Z[stab_index, Q3] = 1 
                    Z[stab_index, Q4] = 1 

                elseif ((i1%2 != 0) && (i2%2 == 0)) || ((i1%2 == 0) && (i2%2 != 0))
                    X[stab_index, Q1] = 1 
                    X[stab_index, Q2] = 1 
                    X[stab_index, Q3] = 1 
                    X[stab_index, Q4] = 1 
                end

                stab_index += 1
            end
        end
    end
    

    #making X & Z into bool
    Z = !=(0).(Z)
    X = !=(0).(X) 

    return Stabilizer(X,Z) #appears twice!

end #function

parity_checks(c::Surface) = checks_s(c)

#-----------------------------------------------------

parity_matrix(c::Surface) = stab_to_gf2(parity_checks(c))

distance(c::Surface) = min(c.lx,c.ly)

code_k(c::Surface) = 1
