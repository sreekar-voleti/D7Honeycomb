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
[ -0.0011 , -0.0084 ,  0.0615 ,  0.0132 ,  0.0109 ],
[ -0.0156 , -0.0059 ,  0.0406 , -0.0241 , -0.0067 ],
[  0.0192 ,  0.0284 , -0.0280 ,  0.0066 ,  0.0109 ],
[  0.0085 ,  0.0007 ,  0.0316 ,  0.0049 ,  0.0113 ],
[  0.0056 , -0.0099 ,  0.0243 ,  0.0144 ,  0.0049 ]
])

onebodymat(t1_11)

# 2. Connecting Vec ( 1.000 0.000 1.587   )

t4_11 = t1_11.T

# 3. Connecting Vec ( 0.500 -0.866 -1.587 )

t2_11 = np.array([
[  0.0106 , -0.0229 , -0.0518 , -0.0015 ,  0.0011 ],
[ -0.0190 ,  0.0124 ,  0.0071 , -0.0112 ,  0.0066 ],
[ -0.0190 , -0.0085 , -0.0280 , -0.0279 ,  0.0112 ],
[ -0.0048 ,  0.0136 , -0.0510 , -0.0133 ,  0.0079 ],
[ -0.0042 ,  0.0048 ,  0.0411 , -0.0001 , -0.0067 ]
])

onebodymat(t2_11)

# 4. Connecting Vec ( -0.500 0.866 1.587 )

t5_11 = t2_11.T

# 5. Connecting Vec ( 0.500 0.866 -1.587 )

t3_11 = np.array([
[ -0.0038 , -0.0029 , -0.0097 , -0.0068 , -0.0041 ],
[ -0.0057 , -0.0079 , -0.0477 , -0.0019 ,  0.0028 ],
[ -0.0002 , -0.0199 , -0.0280 ,  0.0213 , -0.0221 ],
[ -0.0010 ,  0.0229 ,  0.0194 ,  0.0070 ,  0.0211 ],
[ -0.0093 ,  0.0100 , -0.0654 ,  0.0199 ,  0.0076 ]
])

onebodymat(t3_11)

# 4. Connecting Vec ( -0.500 -0.866 1.587 )

t6_11 = t3_11.T

# def constructhopper(T):
#     return np.block([[np.zeros((5,5) , dtype=complex) , T ],[T.T , np.zeros((5,5) , dtype = complex)]])
#
# print("\n" , np.linalg.eigh(constructhopper(t1_11))[0])
# print("\n" , np.linalg.eigh(constructhopper(t2_11))[0])
# print("\n" , np.linalg.eigh(constructhopper(t3_11))[0])
