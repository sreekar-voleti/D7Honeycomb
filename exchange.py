import numpy as np
import sys
import matplotlib.pyplot as plt
from expvals import alpha

def extract3x3(H):
    a = (H[1][2]+H[2][1]+H[0][3]+H[3][0])
    b = (H[1][2]-H[0][3]+H[2][1]-H[3][0])
    c = 2 * (H[0][0]-H[1][1])
    d = -1.0j *(H[1][2]-H[0][3]-H[2][1]+H[3][0])
    e = 2 * (H[2][0]+H[0][2])
    f = -2.0j * (H[2][0]-H[0][2])
    g = -1.0j * (H[2][1]+H[3][0]-H[1][2]-H[0][3])
    h = 2.0 * (H[0][1]+H[1][0])
    p = -2.0j * (H[1][0]-H[0][1])
    return np.array([[a,d,e],[g,b,f],[h,p,c]])

def localcouplingsZ(smat):
    J = (smat[0][0] + smat[1][1])/2.
    A = (smat[0][0] - smat[1][1])/2.
    K = smat[2][2] - J
    Gamma = (smat[0][1] + smat[1][0])/2.
    Gammap1 = (smat[0][2] + smat[2][0])/2.
    Gammap2 = (smat[1][2] + smat[2][1])/2.
    return np.array([J , A , K , Gamma , Gammap1 , Gammap2])

def localcouplingsY(smat):
    J = (smat[0][0] + smat[2][2])/2.
    A = (- smat[0][0] + smat[2][2])/2.
    K = smat[1][1] - J
    Gamma = (smat[2][0] + smat[0][2])/2.
    Gammap1 = (smat[1][2] + smat[2][1])/2.
    Gammap2 = (smat[0][1] + smat[1][0])/2.
    return np.array([J , A , K , Gamma , Gammap1 , Gammap2])

def localcouplingsX(smat):
    J = (smat[1][1] + smat[2][2])/2.
    A = (smat[1][1] - smat[2][2])/2.
    K = smat[0][0] - J
    Gamma = (smat[2][1] + smat[1][2])/2.
    Gammap1 = (smat[1][0] + smat[0][1])/2.
    Gammap2 = (smat[0][2] + smat[2][0])/2.
    return np.array([J , A , K , Gamma , Gammap1 , Gammap2])

def globalcouplingsZ(smat):
    Jxy = (smat[0][0] + smat[1][1])/2.
    Jz = smat[2][2]
    A = (smat[0][0] - smat[1][1])/2.
    B = (smat[0][1] + smat[1][0])/2.
    C = (smat[0][2] + smat[2][0])/2.
    F = (smat[1][2] + smat[2][1])/2.
    return np.array([JXY , JZ, A, B, C, F])

def globalcouplingsY(smat):
    Jxy = (smat[0][0] + smat[2][2])/2.
    Jz = smat[1][1]
    A = (smat[0][0] - smat[1][1])/2.
    B = (smat[0][1] + smat[1][0])/2.
    C = (smat[0][2] + smat[2][0])/2.
    F = (smat[1][2] + smat[2][1])/2.
    return np.array([JXY , JZ, A, B, C, F])

phiv = 0.59

def rotZ(alpha):
    return np.array([[np.cos(alpha) , -np.sin(alpha) , 0],[np.sin(alpha) , np.cos(alpha) , 0],[0 , 0 , 1]])

g_to_l = np.array([
[0.707107, 0., -0.707107],
[-0.408248, 0.816497, -0.408248],
[0.57735, 0.57735, 0.57735]
])

#The passive transformations are done assuming that the Sz is aligned with the vector that the theta and phi indicate

def passive_transformation(smat , trans):
    return trans.dot(smat).dot(trans.T)

Material = sys.argv[1]
Bond_Type = sys.argv[2]
Neighbour = sys.argv[3]
U_index = sys.argv[4]

U = 2.0 + int(U_index) * 0.3

if Bond_Type == 'IP':

    Neibs = [3 , 6 , 3 , 6]

    Matset = np.zeros( (Neibs[int(Neighbour)-1] , 4 , 4) , dtype = complex )

    for k in range(Neibs[int(Neighbour)-1]):
        Matset[k] += np.load(Material+'/'+Material+'_NP/Hameffs/IP_'+Neighbour+'nn_t'+str(k+1)+'_U='+str(U)+'eV.npy')

elif Bond_Type == 'OP':

    Neibs = [1 , 6]

    Matset = np.zeros( (Neibs[int(Neighbour)-1] , 4 , 4) , dtype = complex )

    for k in range(Neibs[int(Neighbour)-1]):
        Matset[k] += np.load(Material+'/'+Material+'_NP/Hameffs/OP_'+Neighbour+'nn_t'+str(k+1)+'_U='+str(U)+'eV.npy')

WhatDo = sys.argv[5]

if WhatDo == 'evs':
    for k in range(Neibs[int(Neighbour)-1]):
        print("\n" , np.linalg.eigh(Matset[k])[0] )


elif WhatDo == 'spinmats':
    print("Global Correction angle: " , np.round(alpha , 3))
    print("Spin Matrices in Khalliulin's GLOBAL basis: ")
    for k in range(Neibs[int(Neighbour)-1]):
        smat = extract3x3(Matset[k])*1e3 #Spin matrix in meV
        smat_g = passive_transformation( smat , rotZ(alpha) ) #Passively Rotating these by alpha to match Khalliulin's basis
        print("\n" , np.round( smat_g , 3))

elif WhatDo == 'globalcouplings':
    print("Using Z bond in plane:")
    x = globalcouplingsZ( rotmat( rotmat( extract3x3(Matset[2])*1e3 , phiv) , 2*np.pi/3))
    y = globalcouplingsZ( rotmat( rotmat( extract3x3(Matset[0])*1e3 , phiv) , -2*np.pi/3))
    z = globalcouplingsZ( rotmat( extract3x3(Matset[1])*1e3 , phiv) )
    avgcouplings = (x+y+z)/3.
    print("Jxy = " , np.round(avgcouplings[0],3))
    print("Jz = " , np.round(avgcouplings[1],3))
    print("A = " , np.round(avgcouplings[2],3))
    print("B = " , np.round(avgcouplings[3],3))
    print("C = " , np.round(avgcouplings[4],3))
    print("F = " , np.round(avgcouplings[5],3))

elif WhatDo == 'rot':
    print( np.round( rotmat(extract3x3(Matset[1]) * 1e3 , 2*np.pi/3) , 3 ) )
    for i in range(len(Matset)):
        print("\n" , np.round(extract3x3(Matset[i]) , 6) * 1e3 )

elif WhatDo == 'localmatrices':
    print("Global Correction angle: " , np.round(alpha , 3))
    print("Spin Matrices in Khalliulin's LOCAL basis: ")
    for i in range(len(Matset)):
        smat = extract3x3(Matset[k])*1e3 #Spin matrix in meV
        smat_g = passive_transformation( smat , rotZ(alpha) ) #Passively Rotating these by alpha to match Khalliulin's basis
        smat_l = passive_transformation( smat_g , g_to_l )
        print("\n" , np.round( smat_l , 3))

elif WhatDo == 'localcouplings':
    y_couplings = localcouplingsY(transform_to_local( extract3x3(Matset[0]) * 1e3 ))
    z_couplings = localcouplingsZ(transform_to_local( extract3x3(Matset[1]) * 1e3 ))
    x_couplings = localcouplingsX(transform_to_local( extract3x3(Matset[2]) * 1e3 ))
    avgcouplings = (x_couplings + y_couplings + z_couplings)/3.
    print("J = " , np.round(avgcouplings[0],3))
    print("A = " , np.round(avgcouplings[1],3))
    print("K = " , np.round(avgcouplings[2],3))
    print("Gamma = " , np.round(avgcouplings[3],3))
    print("Gammap1 = " , np.round(avgcouplings[4],3))
    print("Gammap2 = " , np.round(avgcouplings[5],3))
