import numpy as np
import sys
from params import dist

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

def printcouplings(smat):
    Jxy = np.round( (smat[0][0] + smat[1][1])/2. , 5)
    Jz = np.round( smat[2][2] , 5)
    A = np.round( (smat[0][0] - smat[1][1])/2. , 5)
    B = np.round( (smat[0][1] + smat[1][0])/2. , 5)
    C = np.round( (smat[0][2] + smat[2][0])/2. , 5)
    F = np.round( (smat[1][2] + smat[2][1])/2. , 5)
    Dx = np.round( (smat[0][1] - smat[1][0])/2. , 5)
    Dy = np.round( (smat[2][0] - smat[0][2])/2. , 5)
    Dz = np.round( (smat[1][2] - smat[2][1])/2. , 5)

    print("JXY = " , Jxy)
    print("JZ = " , Jz)
    print("A = " , A)
    print("B = " , B)
    print("C = " , C)
    print("F = " , F)
    print("Dx = " , Dx)
    print("Dy = " , Dy)
    print("Dz = " , Dz)

def storecouplings(smat):
    Jxy = (smat[0][0] + smat[1][1])/2.
    Jz = smat[2][2]
    A = (smat[0][0] - smat[1][1])/2.
    B = (smat[0][1] + smat[1][0])/2.
    C = (smat[0][2] + smat[2][0])/2.
    F = (smat[1][2] + smat[2][1])/2.
    Dx = (smat[0][1] - smat[1][0])/2.
    Dy = (smat[2][0] - smat[0][2])/2.
    Dz = (smat[1][2] - smat[2][1])/2.

    return np.array([ Jxy , Jz , A , B , C , F , Dx , Dy , Dz])

def rotmat(spinmat , alpha):
    return np.transpose(np.array([[np.cos(alpha) , np.sin(-alpha),0],[np.sin(alpha),np.cos(alpha),0],[0,0,1]])).dot(spinmat).dot(np.array([[np.cos(alpha) , np.sin(-alpha),0],[np.sin(alpha),np.cos(alpha),0],[0,0,1]]))

typ = sys.argv[1]
neib = sys.argv[2]
Uinp = sys.argv[3]

U = 2.0 + int(Uinp) * 0.3

if typ == 'IP':

    if neib == '1':

        N_neib = 3

        Mat1 = np.load('Hameffs/IP_1nn_'+str(dist)+'_t1_U='+str(U)+'eV.npy')
        Mat2 = np.load('Hameffs/IP_1nn_'+str(dist)+'_t2_U='+str(U)+'eV.npy')
        Mat3 = np.load('Hameffs/IP_1nn_'+str(dist)+'_t3_U='+str(U)+'eV.npy')

        Matset = [Mat1 , Mat2 , Mat3]

    elif neib == '2':

        N_neib = 6

        Mat1 = np.load('Hameffs/IP_2nn_'+str(dist)+'_t1_U='+str(U)+'eV.npy')
        Mat2 = np.load('Hameffs/IP_2nn_'+str(dist)+'_t2_U='+str(U)+'eV.npy')
        Mat3 = np.load('Hameffs/IP_2nn_'+str(dist)+'_t3_U='+str(U)+'eV.npy')
        Mat4 = np.load('Hameffs/IP_2nn_'+str(dist)+'_t4_U='+str(U)+'eV.npy')
        Mat5 = np.load('Hameffs/IP_2nn_'+str(dist)+'_t5_U='+str(U)+'eV.npy')
        Mat6 = np.load('Hameffs/IP_2nn_'+str(dist)+'_t6_U='+str(U)+'eV.npy')

        Matset = [Mat1 , Mat2 , Mat3 , Mat4 , Mat5 , Mat6]

    elif neib == '3':

        N_neib = 3

        Mat1 = np.load('Hameffs/IP_3nn_'+str(dist)+'_tA_U='+str(U)+'eV.npy')
        Mat2 = np.load('Hameffs/IP_3nn_'+str(dist)+'_tB_U='+str(U)+'eV.npy')
        Mat3 = np.load('Hameffs/IP_3nn_'+str(dist)+'_tC_U='+str(U)+'eV.npy')

        Matset = [Mat1 , Mat2 , Mat3]

    elif neib == '4':

        N_neib = 6

        Mat1 = np.load('Hameffs/IP_4nn_'+str(dist)+'_ta1_U='+str(U)+'eV.npy')
        Mat2 = np.load('Hameffs/IP_4nn_'+str(dist)+'_ta2_U='+str(U)+'eV.npy')
        Mat3 = np.load('Hameffs/IP_4nn_'+str(dist)+'_ta3_U='+str(U)+'eV.npy')
        Mat4 = np.load('Hameffs/IP_4nn_'+str(dist)+'_tb1_U='+str(U)+'eV.npy')
        Mat5 = np.load('Hameffs/IP_4nn_'+str(dist)+'_tb2_U='+str(U)+'eV.npy')
        Mat6 = np.load('Hameffs/IP_4nn_'+str(dist)+'_tb3_U='+str(U)+'eV.npy')

        Matset = [Mat1 , Mat2 , Mat3 , Mat4 , Mat5 , Mat6]

whatdo = sys.argv[4]

if whatdo == 'evs':

    print("\n Eigenvalues of 4x4 matrices: ")

    # for i in range(len(Matset)):
    #     print("\n" , np.linalg.eigh(Matset[i])[0])

    def constructhopper(T):
        return np.block([[np.zeros((4,4) , dtype=complex) , T ],[T , np.zeros((4,4) , dtype = complex)]])
    #
    print("\n" , np.linalg.eigh( constructhopper(Matset[0]) )[0] )
    print("\n" , np.linalg.eigh( constructhopper(Matset[1]) )[0] )
    print("\n" , np.linalg.eigh( constructhopper(Matset[2]) )[0] )
    print("\n" , np.linalg.eigh( constructhopper(Matset[3]) )[0] )
    print("\n" , np.linalg.eigh( constructhopper(Matset[4]) )[0] )
    print("\n" , np.linalg.eigh( constructhopper(Matset[5]) )[0] )

elif whatdo == 'spinmats':

    print("\n Spin Matrices: ")

    for i in range(len(Matset)):
        print("\n" , np.round(extract3x3(Matset[i]) , 6) * 1e3 )

elif whatdo == 'couplings':

    bond = input("Bond Number: ")

    print("\n Couplings (in meV): \n")

    printcouplings(1e3 * np.real(extract3x3(Matset[int(bond)-1])))

elif whatdo == 'rot':

    print( np.round( rotmat(extract3x3(Matset[0]) * 1e3 , 2*np.pi/3) , 3 ) )

    for i in range(len(Matset)):
        print("\n" , np.round(extract3x3(Matset[i]) , 6) * 1e3 )
