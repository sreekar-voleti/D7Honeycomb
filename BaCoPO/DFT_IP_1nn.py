import numpy as np
from params import Vcf

N_neibs = 3

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

# 1. Connecting Vec ( -1.000 , 0.000 , 0.056 )

t1_12 = np.array([
[  0.0987 ,  0.2255 , -0.0194 , -0.0309 , -0.0137 ],
[ -0.2255 , -0.2091 ,  0.0435 ,  0.0040 , -0.0496 ],
[  0.0194 ,  0.0435 , -0.0103 , -0.0031 , -0.0874 ],
[ -0.0309 , -0.0040 ,  0.0031 ,  0.0612 , -0.0033 ],
[  0.0137 , -0.0496 , -0.0874 ,  0.0033 ,  0.0115 ]
])

fix_basis_and_signs(t1_12)

# 2. Connecting Vec ( 0.500 , -0.866 , 0.056 )

t2_12 = np.array([
[  0.1217  , -0.2215 , -0.0096 ,  0.0265 ,  0.0213 ],
[  0.2215  , -0.1827 , -0.0315 , -0.0309 ,  0.0908 ],
[  0.0096  , -0.0315 ,  0.0161 ,  0.0059 , -0.0761 ],
[  0.0265  ,  0.0309 , -0.0059 ,  0.0711 , -0.0043 ],
[ -0.0213  ,  0.0908 , -0.0761 ,  0.0043 , -0.0084 ]
])

fix_basis_and_signs(t2_12)

# 3. Connecting Vec ( 0.500 , 0.866 , 0.056 )

t3_12 = np.array([
[ -0.0363 , -0.0024 , -0.0519 , -0.0106 , -0.0057 ],
[  0.0024 ,  0.1162 ,  0.0031 , -0.0458 , -0.0152 ],
[  0.0519 ,  0.0031 , -0.2813 , -0.0217 ,  0.1553 ],
[ -0.0106 ,  0.0458 ,  0.0216 ,  0.0748 , -0.0078 ],
[  0.0057 , -0.0152 ,  0.1553 ,  0.0078 , -0.1642 ]
])

fix_basis_and_signs(t3_12)

# Testing Matrices

# print(np.linalg.eigh(t1_12)[0])
# print(np.linalg.eigh(t2_12)[0])
# print(np.linalg.eigh(t3_12)[0])
