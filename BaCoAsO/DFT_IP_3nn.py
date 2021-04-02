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
[  0.0283 , -0.0255 ,  0.0005 , -0.0298 ,  0.0060 ],
[  0.0255 , -0.0177 , -0.0026 , -0.0075 , -0.0102 ],
[ -0.0005 , -0.0026 ,  0.0251 , -0.0377 ,  0.0394 ],
[ -0.0298 ,  0.0075 ,  0.0377 , -0.0013 ,  0.0658 ],
[ -0.0060 , -0.0102 ,  0.0394 , -0.0658 ,  0.0653 ]
])

fix_basis_and_signs(tA_12)

# Connecting Vec -1.000 1.732 0.056

tB_12 = np.array([
[  0.0242 ,  0.0243 , -0.0151 ,  0.0084 , -0.0129 ],
[ -0.0243 , -0.0226 ,  0.0094 , -0.0107 ,  0.0141 ],
[  0.0151 ,  0.0094 ,  0.0051 ,  0.0363 ,  0.0260 ],
[  0.0084 ,  0.0107 , -0.0363 , -0.0019 , -0.0661 ],
[  0.0129 ,  0.0141 ,  0.0260 ,  0.0661 ,  0.0856 ]
])

fix_basis_and_signs(tB_12)

# Connecting Vec 2.000 0.000 0.056

tC_12 = np.array([
[  0.0027 ,  0.0007 , -0.0015 ,  0.0134 ,  0.0086 ],
[ -0.0007 ,  0.0061 ,  0.0000 ,  0.0137 ,  0.0077 ],
[  0.0015 ,  0.0000 , -0.0258 ,  0.0078 ,  0.0149 ],
[  0.0134 , -0.0137 , -0.0078 , -0.1327 ,  0.0042 ],
[ -0.0086 ,  0.0077 ,  0.0149 , -0.0042 , -0.0646 ]
])

fix_basis_and_signs(tC_12)

# Testing Matrices

# =============================================================================
# print(np.linalg.eigh(tA_12)[0])
# print(np.linalg.eigh(tB_12)[0])
# print(np.linalg.eigh(tC_12)[0])
# =============================================================================
