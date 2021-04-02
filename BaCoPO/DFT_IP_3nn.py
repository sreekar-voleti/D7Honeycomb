import numpy as np
from params import Vcf

# ONSITE MATRIX

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

# Connecting Vec -1.000 -1.732 0.056

tA_12 = np.array([
[  0.0197 ,  0.0014 ,  0.0042 , -0.0189 , -0.0236 ],
[ -0.0014 ,  0.0243 , -0.0044 , -0.0360 , -0.0084 ],
[ -0.0042 , -0.0044 , -0.0275 ,  0.0141 ,  0.0087 ],
[ -0.0189 ,  0.0360 , -0.0141 , -0.1278 ,  0.0055 ],
[  0.0236 , -0.0084 ,  0.0087 , -0.0055 , -0.0463 ]
])

fix_basis_and_signs(tA_12)

# Connecting Vec -1.000 1.732 0.056

tB_12 = np.array([
[  0.0342 , -0.0222 , 0.0116 ,  0.0292 , -0.0151 ],
[  0.0222 , -0.0062 , 0.0014 , -0.0117 ,  0.0200 ],
[ -0.0116 ,  0.0014 , 0.0121 , -0.0358 ,  0.0298 ],
[  0.0292 ,  0.0117 , 0.0358 , -0.0059 ,  0.0659 ],
[  0.0151 ,  0.0200 , 0.0298 , -0.0659 ,  0.0809 ]
])

fix_basis_and_signs(tB_12)

# Connecting Vec 2.000 0.000 0.056

tC_12 = np.array([
[  0.0276 ,  0.0208 ,  0.0233 , -0.0056 ,  0.0066 ],
[ -0.0208 , -0.0118 , -0.0054 , -0.0024 , -0.0211 ],
[ -0.0233 , -0.0054 , -0.0072 ,  0.0348 ,  0.0172 ],
[ -0.0056 ,  0.0024 , -0.0348 , -0.0038 , -0.0639 ],
[ -0.0066 , -0.0211 ,  0.0172 ,  0.0639 ,  0.1013 ]
])

fix_basis_and_signs(tC_12)

# Testing Matrices

# =============================================================================
# print(np.linalg.eigh(tA_12)[0])
# print(np.linalg.eigh(tB_12)[0])
# print(np.linalg.eigh(tC_12)[0])
# =============================================================================
