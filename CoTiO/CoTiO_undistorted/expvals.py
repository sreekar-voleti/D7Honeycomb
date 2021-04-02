import binfuncts as bf
import numpy as np
from Amatrices import A
from params import N , L

def mat(A , s , a , b , L , basis):

    proxy = bf.binstate(s, L)
    count = 0

    if a != b:
        if ((proxy[b] == '1') and (proxy[a] == '0')) :

            count += bf.fillcounter(s , b , L)
            s = bf.flip(s , b , L)
            count += bf.fillcounter(s , a , L)
            s = bf.flip(s , a , L)

            return np.searchsorted(basis , s) , (-1)**count * A[a][b]

        else:
            return 0 , 0
    elif a == b:
        if proxy[b] == '1':
            count += bf.fillcounter(s , b , L)
            s = bf.flip(s , b , L)
            count += bf.fillcounter(s , a , L)
            s = bf.flip(s , a , L)

            return np.searchsorted(basis , s) , (-1)**count * A[a][b]

        else:
            return 0 , 0

def mat2(A , s , a , b , L , basis):

    proxy = bf.binstate(s, L)
    count = 0

    if (( proxy[b] == '1' ) and ( proxy[a] == '0' )) :

        count += bf.fillcounter(s , b , L)
        s = bf.flip(s , b , L)
        count += bf.fillcounter(s , a , L)
        s = bf.flip(s , a , L)

        return np.searchsorted(basis , s) , (-1)**count * A[b][a]
    else:
        return s , 0

BasisN = bf.makebasis(L , N)

def LPS(g , i):
    if i == 0:
        return np.kron(bf.Lx , np.identity(2)) + g*np.kron(np.identity(5) , bf.Sx)
    elif i == 1:
        return np.kron(bf.Ly , np.identity(2)) + g*np.kron(np.identity(5) , bf.Sy)
    elif i == 2:
        return np.kron(bf.Lz , np.identity(2)) + g*np.kron(np.identity(5) , bf.Sz)


g = 1 # This is actually g/2 since S is sigma/2

Ax = LPS(g , 0)
Ay = LPS(g , 1)
Az = LPS(g , 2)

BL = len(BasisN)

L1 = np.zeros((BL,BL) , dtype = complex)
L2 = np.zeros((BL,BL) , dtype = complex)
L3 = np.zeros((BL,BL) , dtype = complex)

for i in range(BL):
    for a in range(L):
        for b in range(L):
            boopx = mat(Ax , BasisN[i] , a , b , L , BasisN)
            L1[boopx[0]][i] += boopx[1]
            boopy = mat(Ay , BasisN[i] , a , b , L , BasisN)
            L2[boopy[0]][i] += boopy[1]
            boopz = mat(Az , BasisN[i] , a , b , L , BasisN)
            L3[boopz[0]][i] += boopz[1]

check = L1 - np.conjugate(L1.T)

from params import lam , Jh , Jp , dist , B

U = 3.2
Up = U - 2*Jh

print("For Co1:")

Co = 1

SpecN = bf.SolveSpectrumfull(L , N , BasisN , lam , U , Jh , Jp , Up , -B)

psi_x = 1/np.sqrt(2.) * ( SpecN[1][:,0] + SpecN[1][:,1] )

lx = np.conjugate(psi_x).dot(L1).dot(psi_x)
ly = np.conjugate(psi_x).dot(L2).dot(psi_x)
lz = np.conjugate(psi_x).dot(L3).dot(psi_x)

mu_x = np.real( np.array([lx,ly,lz])/np.linalg.norm(np.array([lx,ly,lz])) )

print("mu_x: " , mu_x)

psi = 1/np.sqrt(2.) * ( SpecN[1][:,0] + 1j * SpecN[1][:,1] )

lx = np.conjugate(psi).dot(L1).dot(psi)
ly = np.conjugate(psi).dot(L2).dot(psi)
lz = np.conjugate(psi).dot(L3).dot(psi)

mu_y = np.real( np.array([lx,ly,lz])/np.linalg.norm(np.array([lx,ly,lz])) )

print("mu_y: " , mu_y)

psi = SpecN[1][:,0]

lx = np.conjugate(psi).dot(L1).dot(psi)
ly = np.conjugate(psi).dot(L2).dot(psi)
lz = np.conjugate(psi).dot(L3).dot(psi)

mu_z = np.real( np.array([lx,ly,lz])/np.linalg.norm(np.array([lx,ly,lz])) )

print("mu_z: " , mu_z)

print(np.dot(mu_x,mu_y))
print(np.dot(mu_x,mu_z))
print(np.dot(mu_z,mu_y))
