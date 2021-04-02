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

A_Co2 = np.array([
[ -1.0968 ,  0.0008 , -0.0103 , -0.0040 , -0.0926 ],
[  0.0008 , -1.1642 , -0.0039 , -0.0312 , -0.0063 ],
[ -0.0103 , -0.0039 , -0.8765 ,  0.0044 ,  0.4380 ],
[ -0.0040 , -0.0312 ,  0.0044 , -0.1429 , -0.0033 ],
[ -0.0926 , -0.0063 ,  0.4380 , -0.0033 , -0.4133 ]
])

A_Co1 = A_Co2.copy()

A_Co1[0,:] *= -1 # xy row
A_Co1[:,0] *= -1 # xy column
A_Co1[3,:] *= -1 # yz row
A_Co1[:,3] *= -1 # yz column

fixbasis(A_Co1)
fixbasis(A_Co2)

print(A_Co1)

evals , evecs = np.linalg.eigh(A_Co1)

def newA(d):
    g1 = (evals[0]+evals[1])/2.
    g2 = evals[2] + d[0]
    g3 = (evals[3]+evals[4])/2. + d[1]
    return g1*(np.outer(evecs[:,0],evecs[:,0])+np.outer(evecs[:,1],evecs[:,1])) + g2*(np.outer(evecs[:,2],evecs[:,2])) + g3*(np.outer(evecs[:,3],evecs[:,3])+np.outer(evecs[:,4],evecs[:,4]))

A_Co1_new = newA([-0.03 , -0.1035]) #Played with the A matrices so that the spectrum matches somewhat to the Neutron data

# print(A_Co1)
CFMat = np.zeros((5,5) , dtype = complex)

CFMat[3][3] += Vcf
CFMat[4][4] += Vcf
