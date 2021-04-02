import numpy as np
from params import Vcf

# ONSITE MATRIX

#The following changes are made to ALL matrices (both rows and columns)

def fixbasis(A):

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

A_Co1 = np.array([
[ -3.4155 ,  0.0008 ,  0.0157 ,  0.0017 , -0.0900 ],
[  0.0008 , -3.5049 ,  0.0008 , -0.0717 ,  0.0034 ],
[  0.0157 ,  0.0008 , -3.1860 ,  0.0249 ,  0.5051 ],
[  0.0017 , -0.0717 ,  0.0249 , -2.2912 , -0.0136 ],
[ -0.0900 ,  0.0034 ,  0.5051 , -0.0136 , -2.5757 ]
])

A_Co2 = A_Co1.copy()

A_Co2[0,:] *= -1 # xy row
A_Co2[:,0] *= -1 # xy column
A_Co2[3,:] *= -1 # yz row
A_Co2[:,3] *= -1 # yz column

fixbasis(A_Co1)
fixbasis(A_Co2)

print(A_Co1)

CFMat = np.zeros((5,5) , dtype = complex)

CFMat[3][3] += Vcf
CFMat[4][4] += Vcf
