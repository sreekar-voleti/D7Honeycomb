import numpy as np
from params import Vcf

N_neibs = 1

#The following changes are made to ALL matrices (both rows and columns)

def onebodymat(A):

    #The basis they gave us the matrix in is different so we need to interchange some rows and columns
    ## Their basis: {xy , yz , 3z2-1 , xz , x2-y2}
    ## Our basis:   {yz , xz , xy , x2-y2 , 3z2-1}

    # Interchanging rows
    A[[0,1],:] = A[[1,0],:] # After this : {yz , xy , 3z2-1 , xz , x2-y2}
    A[[1,3],:] = A[[3,1],:] # After this : {yz , xz , 3z2-1 , xy , x2-y2}
    A[[2,3],:] = A[[3,2],:] # After this : {yz , xz , xy , 3z2-1 , x2-y2}
    A[[3,4],:] = A[[4,3],:] # After this : {yz , xz , xy , x2-y2 , 3z2-1} (target)

    # Interchamging columns (same sequence as above, but with columns)
    A[:,[0,1]] = A[:,[1,0]] # After this : {yz , xy , 3z2-1 , xz , x2-y2}
    A[:,[1,3]] = A[:,[3,1]] # After this : {yz , xz , 3z2-1 , xy , x2-y2}
    A[:,[2,3]] = A[:,[3,2]] # After this : {yz , xz , xy , 3z2-1 , x2-y2}
    A[:,[3,4]] = A[:,[4,3]] # After this : {yz , xz , xy , x2-y2 , 3z2-1} (target)

#Extra set of sign changes on the ROWS only of the hopping matrices before shuffling rows and columns

def fix_basis_and_signs(A):

    #Flip the signs of the ROWS on the hopping matrices, taking the relative sign into account
    ## Flip before shuffling the rows and columns

    A[0,:] *= -1
    A[3,:] *= -1

    onebodymat(A)

# HOPPING MATRICES

# 1. Connecting Vec ( 0.000 0.000 -1.381  )

t_12 = np.array([
[  0.0017 ,  0.0084 ,  0.0000 , -0.0193 ,  0.0000 ],
[  0.0084 ,  0.0404 ,  0.0000 ,  0.0000 , -0.0193 ],
[  0.0000 ,  0.0000 , -0.0562 ,  0.0000 ,  0.0000 ],
[ -0.0193 ,  0.0000 ,  0.0000 ,  0.0404 , -0.0084 ],
[  0.0000 , -0.0193 ,  0.0000 , -0.0084 ,  0.0017 ]
])

onebodymat(t_12)

# Testing Matrices

# print(np.linalg.eigh(t_12)[0])
