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

# 1. Connecting Vec ( -1.000 0.000 -1.587   )

t1_11 = np.array([
[ -0.0109 ,  0.0201 ,  0.0184 , -0.0169 , -0.0022 ],
[  0.0016 , -0.0137 , -0.0085 ,  0.0071 ,  0.0056 ],
[  0.0017 ,  0.0062 ,  0.0013 , -0.0099 , -0.0019 ],
[  0.0104 , -0.0295 , -0.0036 ,  0.0109 ,  0.0275 ],
[ -0.0117 ,  0.0188 , -0.0051 , -0.0542 ,  0.0071 ]
])

fix_basis_and_signs(t1_11)

# 2. Connecting Vec ( 1.000 0.000 1.587   )

t4_11 = t1_11.T

# 3. Connecting Vec ( 0.500 -0.866 -1.587 )

t2_11 = np.array([
[ -0.0027 , -0.0013 ,  0.0179 ,  0.0102 ,  0.0031 ],
[  0.0063 ,  0.0052 ,  0.0090 ,  0.0049 ,  0.0085 ],
[  0.0024 ,  0.0079 , -0.0153 , -0.0224 , -0.0022 ],
[  0.0093 , -0.0146 ,  0.0534 ,  0.0174 ,  0.0232 ],
[ -0.0087 , -0.0278 ,  0.0160 , -0.0180 , -0.0098 ]
])

fix_basis_and_signs(t2_11)

# 4. Connecting Vec ( -0.500 0.866 1.587 )

t5_11 = t2_11.T

# 5. Connecting Vec ( 0.500 0.866 -1.587 )

t3_11 = np.array([
[ -0.0216 , -0.0065 , -0.0140 ,  0.0042 , -0.0129 ],
[ -0.0167 ,  0.0070 , -0.0049 , -0.0046 ,  0.0106 ],
[  0.0264 ,  0.0079 ,  0.0144 , -0.0155 ,  0.0063 ],
[  0.0005 ,  0.0131 ,  0.0208 , -0.0089 ,  0.0358 ],
[  0.0276 ,  0.0124 ,  0.0234 , -0.0308 ,  0.0038 ]
])

fix_basis_and_signs(t3_11)

# 4. Connecting Vec ( -0.500 -0.866 1.587 )

t6_11 = t3_11.T

# def constructhopper(T):
#     return np.block([[np.zeros((5,5) , dtype=complex) , T ],[T.T , np.zeros((5,5) , dtype = complex)]])
#
# print("\n" , np.linalg.eigh(constructhopper(t1_11))[0])
# print("\n" , np.linalg.eigh(constructhopper(t2_11))[0])
# print("\n" , np.linalg.eigh(constructhopper(t3_11))[0])
