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

#They only gave us one A matrix
A = np.array([
[ -1.6965 ,  0.0834 ,  0.0348 , -0.0423 ,  0.2397 ],
[  0.0834 , -1.4492 ,  0.4427 ,  0.0502 ,  0.2318 ],
[  0.0348 ,  0.4427 , -0.8755 ,  0.1956 , -0.0963 ],
[ -0.0423 ,  0.0502 ,  0.1956 , -1.6728 , -0.2138 ],
[  0.2397 ,  0.2318 , -0.0963 , -0.2138 , -0.7623 ]
])

fixbasis(A)

print(A)
