#Apparently Julia doesn't require importing things here? Only in the main module file. Okay 

#return the eigenvalues of the the symbolic matrix 

#NOTE: Currently can only handle the following 
"""
  - 2x2 matrices

  - Will write explicit code for 

  3x3 matrices and maybe 4x4 

  - 5x5 and up is not possible in the explicit sense, but borrowing the path of matlab/mathematica 
  of returning the polynomial equation that yields the eigenvalue is possible. 

  - Solution in specialized cases will also be developed, especially if one can know that 
  that solution comes from jacobians or etc. Mainly developing this for conservation laws which are generally
  fairly structured.
"""
function eigenvals(A::AbstractMatrix{<:Union{Symbolic,Num}})
  (m,n) = size(A)
  #We dispatch into specialzied funcitons based on the size of the matrix
  if (m ==! n)
    error("This function requires a square matrix")
  end

  if (m == 2)
    return eigenvals_2x2(A)
  else
    error("This function only works for 2x2 matrices ")
  end
end

function eigenvals_2x2(::AbstractMatrix{<:Union{Symbolic,Num}})
  (m,n) = size(A)
  if (m == 2 && n == 2)
    error("This function is only for 2x2 matrices")
  end
  #Trace, Determinant
  T = A[1,1] + A[2,2]
  D = A[1,1]*A[2,2]+A[1,2]*A[2,1]

  L₁ = T/2 + sqrt(T^2/4-D)
  L₂ = T/2 - sqrt(T^2/4-D)
  return [L_1,L_2]
end

function eigenvecs(A::AbstractMatrix{<:Union{Symbolic,Num}})
  (m,n) = size(A)
  #We dispatch into specialzied funcitons based on the size of the matrix
  if (m ==! n)
    error("This function requires a square matrix")
  end

  if (m == 2)
    return eigenvecs_2x2(A)
  else
    error("This function only works for 2x2 matrices ")
  end
end

function eigenvecs_2x2(A::AbstractMatrix{<:Union{Symbolic,Num}})
  (m,n) = size(A)
  if (m == 2 && n == 2)
    error("This function is only for 2x2 matrices")
  end
  L = eigenvals_2x2(A)
  c = simplify(A[2,1]) #need to check if this is zero 

  #these could be switched based on the more stable choice later 
  if (c ==! 0)
    return [L[1]-A[2,2] c; L[2]-A[2,2] c]
  else
    b = simplify(A[1,2])
    if b ==! 0
      return [b L[1]-A[1,1];b  L[2]-A[1,1]]
    else 
      #this should have the same type as the input.... 
      return I{Num}
  end
end