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

    return (q1,q2,q3,q4)
end

function plaquette_to_qubit_indices_X(rown,columnn,j) 
    
    #special cases: corners
    if (rown == 1 && column == 1) || (rown == 2*j && column == 1) 
        q1 = [rown,columnn]
        q2 = [rown+1,columnn]
        return q1,q2
    elseif (row == 1) && (column == j)
        q1 = [rown,columnn]
        q2 = [rown+1,columnn+1]
        return q1,q2        
    elseif (rown %2 == 0) && (columnn!= j) && (columnn!= 1) && (rown < j^2)
        q1 = [rown,columnn]
        q2 = [rown+1,columnn]
        q3 = [rown+1,columnn-1]
        q4 = [rown+2,columnn]
        return q1,q2,q3,q4
    end
    
end

function grid_index_to_linear_index(q) 
    a,b = q 
    
    if a%2 != 0
        l = a*4+b
    elseif a%2 == 0 
        l = 3+a*3+b 
    end

    return l
end

function checks(c::Toric) 

    #Z checks
    for i1 in range(1,2,c.l)
        for i2 in range(1,c.l)
            if i2+2 <= (c.l^2 + 1)
                q1, q2, q3, q4 = plaquette_to_qubit_indices_Z(i1,i2)
                # q1 is a tuple
                Q1 = grid_index_to_linear_index(q1) # Q1 is an Int
                #Z[stab_index, Q1] = true
            end
        end
    end

    #X checks
    for i1 in range(1,2,c.l)
        #=
        if (i1 == 1 && i2 == 1) || (i1 == 2*j && i2 == 1) || (i1 == 1) && (i2 == j)
            q1, q2 == plaquette_to_qubit_indices_X(i1,i2,c.l) #may not work for corners
        else
            q1, q2, q3, q4 == plaquette_to_qubit_indices_X(i1,i2,c.l) #may not work for corners
        end
        =#
        q1, q2, q3, q4 == plaquette_to_qubit_indices_X(i1,i2,c.l) #may not work for corners
        # q1 is a tuple
        Q1 = grid_index_to_linear_index(q1) # Q1 is an Int
        #X[stab_index, Q1] = true
    end

    #return Stabilizer(X,Z)

end #function

parity_checks(c::Toric) = checks(c)
#-----------------------------------------------------

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

distance(c::Toric) = sqrt(c.n)

code_k(c::Toric) = 1