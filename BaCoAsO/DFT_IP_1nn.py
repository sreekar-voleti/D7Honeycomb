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
[ -0.0442 , -0.0009 ,  0.0476 ,  0.0076 , -0.0052 ],
[  0.0009 ,  0.0961 , -0.0033 ,  0.0558 ,  0.0123 ],
[ -0.0476 , -0.0033 , -0.1699 , -0.0221 ,  0.1650 ],
[  0.0076 , -0.0558 ,  0.0221 ,  0.0240 , -0.0129 ],
[  0.0052 ,  0.0123 ,  0.1650 ,  0.0129 , -0.1668 ]
])

fix_basis_and_signs(t1_12)

# 2. Connecting Vec ( 0.500 , -0.866 , 0.056 )

t2_12 = np.array([
[  0.0849 ,  0.1913 ,  0.0034 ,  0.0711 , -0.0005 ],
[ -0.1913 , -0.1581 , -0.0147 , -0.0391 ,  0.0531 ],
[ -0.0034 , -0.0147 ,  0.0038 , -0.0034 , -0.0649 ],
[  0.0711 ,  0.0391 ,  0.0034 ,  0.0226 , -0.0043 ],
[  0.0005 ,  0.0531 , -0.0649 ,  0.0043 ,  0.0416 ]
])

fix_basis_and_signs(t2_12)

# 3. Connecting Vec ( 0.500 , 0.866 , 0.056 )

t3_12 = np.array([
[  0.0883 , -0.1892 ,  0.0012 , -0.0676 , -0.0240 ],
[  0.1892 , -0.1539 ,  0.0095 , -0.0133 , -0.0808 ],
[ -0.0012 ,  0.0095 ,  0.0176 ,  0.0115 , -0.0560 ],
[ -0.0676 ,  0.0133 , -0.0115 ,  0.0250 ,  0.0087 ],
[  0.0240 , -0.0808 , -0.0560 , -0.0087 ,  0.0292 ]
])

fix_basis_and_signs(t3_12)

# Testing Matrices

# print(np.linalg.eigh(t1_12)[0])
# print(np.linalg.eigh(t2_12)[0])
# print(np.linalg.eigh(t3_12)[0])
