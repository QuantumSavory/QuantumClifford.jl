import .ECC

using LinearAlgebra
using .ECC

struct Toric <: AbstractECC 
    l
end 

#Qubit n -------------------------------------

code_n(c::Toric) = 2*c.l*(c.l+1)

#Parity checks ----------------------------------

function plaquette_to_qubit_indices_Z(rown,columnn) 
    
    q1 = [rown,columnn]
    q2 = [rown+1,columnn]
    q3 = [rown+1,columnn+1]
    q4 = [rown+2,columnn]
    return q1,q2,q3,q4

end

function plaquette_to_qubit_indices_X_q4(row::Int64,column::Int64,j::Int64) 
    
    #middle cross(s)
    q1 = [row,column]
    q2 = [row+1,column]
    q3 = [row+1,column-1]
    q4 = [row+2,column]
    return q1,q2,q3,q4

end

function plaquette_to_qubit_indices_X_q2(row::Int64,column::Int64,j::Int64) 
    
    #special cases: corners
    if (row == 1 && column == 1) || (row == 2*j && column == 1) 
        q1 = [row,column]
        q2 = [row+1,column]
        return q1,q2

    elseif (row == 1) && (column == j)
        q1 = [row,column]
        q2 = [row+1,column+1]
        return q1,q2   
    
    elseif (row == j*2 -2) && (column == j)
        q1 = [row,column]
        q2 = [row+1,column]
        return q1,q2 
    return q1,q2 

    end
end

function grid_index_to_linear_index_toric(c,q)::Int8 #to test with big grids
    a,b = q 
    
    if a<= 2
        l = c.l*(a-1) +b 
    elseif a%2 == 0
        l = c.l*(a-1) +b +a/2 -1
    else
        l = c.l*(a-1) +b +(a-1)/2 
    end
    return l
end

function checks_t(c::Toric) 
    Z = zeros(code_n(c),code_n(c))
    X = zeros(code_n(c),code_n(c))
    z_stab_index = 1
    x_stab_index = 1

    #Z checks
    for i1 in range(start=1,step=2,length=c.l)
        for i2 in range(1,c.l)
            if i1+2 <= (code_n(c)) && i1%2  != 0
                q1, q2, q3, q4 = plaquette_to_qubit_indices_Z(i1,i2)
                # q1 is a tuple
                Q1 = grid_index_to_linear_index_toric(c,q1) # Q1 is an Int
                Q2 = grid_index_to_linear_index_toric(c,q2)
                Q3 = grid_index_to_linear_index_toric(c,q3)
                Q4 = grid_index_to_linear_index_toric(c,q4)
                
                Z[z_stab_index, Q1] = 1 
                Z[z_stab_index, Q2] = 1 
                Z[z_stab_index, Q3] = 1 
                Z[z_stab_index, Q4] = 1 
                z_stab_index += 1
            end
        end
    end
    
    #X checks
    for i1 in range(start=1,length=c.l)
        for i2 in range(1,c.l)        

            if ((i1 == 1 && i2 == 1) || (i1 == 2*c.l+1 && i2 == 1)) || ((i1 == 1) && (i2 == c.l)) 
                q1, q2 = plaquette_to_qubit_indices_X_q2(i1,i2,c.l) 
                # q1 is a tuple
                Q1 = grid_index_to_linear_index_toric(c,q1) # Q1 is an Int
                Q2 = grid_index_to_linear_index_toric(c,q2)

                X[x_stab_index, Q1] = 1 
                X[x_stab_index, Q2] = 1 

            elseif (i1 %2 == 0) && (i2 <= c.l)   && (i2 > 1) 

                q1, q2, q3, q4 = plaquette_to_qubit_indices_X_q4(i1,i2,c.l) 
    
                Q1 = grid_index_to_linear_index_toric(c,q1) # Q1 is an Int
                Q2 = grid_index_to_linear_index_toric(c,q2)
                Q3 = grid_index_to_linear_index_toric(c,q3)
                Q4 = grid_index_to_linear_index_toric(c,q4)
                X[x_stab_index, Q1] = 1 
                X[x_stab_index, Q2] = 1 
                X[x_stab_index, Q3] = 1 
                X[x_stab_index, Q4] = 1 

            end
            
            x_stab_index +=1
        end
    end

    @show(size(Z))
    @show(size(X))

    #making X & Z into bool
    Z = !=(0).(Z)
    X = !=(0).(X) 

    return Stabilizer(X,Z) #appears twice!

end #function

parity_checks(c::Toric) = checks_t(c)
#-----------------------------------------------------

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

distance(c::Toric) = sqrt(c.n)

code_k(c::Toric) = 1