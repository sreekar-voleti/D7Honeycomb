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

# 1. Connecting Vec ( 1.000 , 0.000 , 0.056 )

t1_12 = np.array([
[ -0.0065 ,  0.0356 , -0.0483 ,  0.0277 ,  0.0266 ],
[ -0.0356 , -0.1549 ,  0.0816 ,  0.0730 , -0.0244 ],
[  0.0483 ,  0.0816 , -0.1181 , -0.0505 ,  0.0344 ],
[  0.0277 , -0.0730 ,  0.0505 , -0.0751 , -0.0288 ],
[ -0.0266 , -0.0244 ,  0.0344 ,  0.0288 , -0.1139 ]
])

fix_basis_and_signs(t1_12)

# 2. Connecting Vec ( -0.500 , 0.866 , 0.056 )

t2_12 = np.array([
[ -0.0165 ,  0.0620 , -0.0316 , 0.0385 ,  0.0376 ],
[ -0.0620 ,  0.0325 , -0.0702 , 0.0408 , -0.0340 ],
[  0.0316 , -0.0702 , -0.0936 , 0.0636 , -0.0520 ],
[  0.0385 , -0.0408 , -0.0636 , 0.1900 , -0.0517 ],
[ -0.0376 , -0.0340 , -0.0520 , 0.0517 , -0.0706 ]
])

fix_basis_and_signs(t2_12)

# 3. Connecting Vec ( -0.500 , -0.866 , 0.056 )

t3_12 = np.array([
[  0.1741 , -0.0695 , -0.0397 , -0.0671 , -0.0265 ],
[  0.0695 , -0.0249 , -0.1176 ,  0.0101 , -0.0283 ],
[  0.0397 , -0.1176 , -0.0595 , -0.0729 ,  0.0275 ],
[ -0.0671 , -0.0101 ,  0.0729 ,  0.0077 , -0.0069 ],
[  0.0265 , -0.0283 ,  0.0275 ,  0.0069 , -0.0389 ]
])

fix_basis_and_signs(t3_12)

# Testing Matrices

print( np.array_equal(t1_12 , t1_12.T) and np.array_equal(t2_12 , t2_12.T) and np.array_equal(t3_12 , t3_12.T) )

print(np.linalg.eigh(t1_12)[0])
print(np.linalg.eigh(t2_12)[0])
print(np.linalg.eigh(t3_12)[0])
