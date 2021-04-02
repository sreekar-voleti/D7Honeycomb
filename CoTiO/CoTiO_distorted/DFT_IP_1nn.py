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
    # A[:,0] *= -1
    # A[:,3] *= -1

    onebodymat(A)

# HOPPING MATRICES

# 1. Connecting Vec ( 1.000 , 0.000 , 0.056 ) (Y bond)

t1_12 = np.array([
[ -0.0793 ,  0.0026 ,  0.0051 , -0.0194 ,  0.0515 ],
[ -0.0026 ,  0.0012 , -0.0134 , -0.0100 , -0.0243 ],
[ -0.0051 , -0.0134 , -0.1146 , -0.0026 ,  0.0811 ],
[ -0.0194 ,  0.0100 ,  0.0026 ,  0.0427 , -0.0268 ],
[ -0.0515 , -0.0243 ,  0.0811 ,  0.0268 , -0.1734 ]
])

fix_basis_and_signs(t1_12)

# 2. Connecting Vec ( -0.500 , 0.866 , 0.056 ) (Z bond)

t2_12 = np.array([
[  0.0744 ,  0.1259 ,  0.0154 ,  0.0691 , -0.0405 ],
[ -0.1259 , -0.0779 ,  0.0287 , -0.0214 , -0.0006 ],
[ -0.0154 ,  0.0287 ,  0.0376 ,  0.0123 , -0.0271 ],
[  0.0691 ,  0.0214 , -0.0123 ,  0.0991 , -0.0119 ],
[  0.0405 , -0.0006 , -0.0271 ,  0.0119 , -0.0364 ]
])

fix_basis_and_signs(t2_12)

# 3. Connecting Vec ( -0.500 , -0.866 , 0.056 ) (X bond)

t3_12 = np.array([
[  0.0684 , -0.1233 ,  0.0073 , -0.0236 , -0.0575 ],
[  0.1233 , -0.0803 , -0.0415 , -0.0411 , -0.0163 ],
[ -0.0073 , -0.0415 , -0.0093 , -0.0433 , -0.0641 ],
[ -0.0236 ,  0.0411 ,  0.0433 ,  0.0588 ,  0.0143 ],
[  0.0575 , -0.0163 , -0.0641 , -0.0143 , -0.0333 ]
])

fix_basis_and_signs(t3_12)

# Testing Matrices

# print(np.linalg.eigh(t1_12)[0])
# print(np.linalg.eigh(t2_12)[0])
# print(np.linalg.eigh(t3_12)[0])
