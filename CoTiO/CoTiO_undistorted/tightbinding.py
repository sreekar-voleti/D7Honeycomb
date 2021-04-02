import numpy as np
import matplotlib.pyplot as plt
from Amatrices import A
import DFT_IP_1nn as dft
from scipy.spatial.transform import Rotation as R
import binfuncts as bf

lam = 0.1
ons = 1.

cv1 = np.array([1. , 0. , 0.205])
cv2 = np.array([-0.5  , 0.866 , 0.205])
cv3 = np.array([-0.5  , -0.866 , 0.205])

lv1 = cv1 - cv3
lv2 = cv1 - cv2

d1 = cv2 - cv1
d2 = cv3 - cv1

nnd = 1.
a = np.linalg.norm(d1)

t1_12 = dft.t1_12
t2_12 = dft.t2_12
t3_12 = dft.t3_12

A1 = lv1
A2 = lv2

# A1 = np.array([a,0,0])
# A2 = np.array([a/2., np.sqrt(3.)*a/2.,0])

def rotate(vec , theta):

    rotation_axis = np.array([0, 0, 1])
    rotation_vector = np.radians(theta) * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    rotated_vec = rotation.apply(vec)

    return rotated_vec

def vhat(v):
    return v / np.linalg.norm(v)

cv1 =vhat(A1 + A2)*nnd
cv2 = vhat(A2 - 2*A1)*nnd
cv3 = vhat(A1 - 2*A2)*nnd

theta = 240
q = 5
N = 10

t = 1.

if N == 10:
    def FillKHam(Ham , k , kvec , N , ons , lam):
        for i in range(N):
            for j in range(N):
                Ham[k][i][j] += (ons * np.kron(A , np.identity(2)) + lam* ( bf.LS(0) + bf.LS(1) + bf.LS(2) ) )[i][j]
                Ham[k][N + i][N + j] += (ons * np.kron(A , np.identity(2)) + lam * ( bf.LS(0) + bf.LS(1) + bf.LS(2) ) )[i][j]
                Ham[k][i][N + j] += np.kron(t1_12 , np.identity(2))[i][j]
                Ham[k][N + i][j] += np.kron(t1_12.T , np.identity(2))[i][j]
                Ham[k][i][N + j] += np.exp(1j * kvec.dot(-A1)) * np.kron(t2_12 , np.identity(2))[i][j]
                Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(-A1)) * np.kron(t2_12.T , np.identity(2))[i][j]
                Ham[k][i][N + j] += np.exp(1j * kvec.dot(-A2)) * np.kron(t3_12 , np.identity(2))[i][j]
                Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(-A2)) * np.kron(t3_12.T , np.identity(2))[i][j]

if N == 2:
    def FillKHam(Ham , k , kvec , N):
        for i in range(N):
            for j in range(N):
# =============================================================================
#                 Ham[k][i][N + j] += np.exp(1j * kvec.dot(cv1)) * t
#                 Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(cv1)) * t
#                 Ham[k][i][N + j] += np.exp(1j * kvec.dot(cv2)) * t
#                 Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(cv2)) * t
#                 Ham[k][i][N + j] += np.exp(1j * kvec.dot(cv3)) * t
#                 Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(cv3)) * t
# =============================================================================
                Ham[k][i][N + j] += t
                Ham[k][N + i][j] += t
                Ham[k][i][N + j] += np.exp(1j * kvec.dot(-A1)) * t
                Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(-A1)) * t
                Ham[k][i][N + j] += np.exp(1j * kvec.dot(-A2)) * t
                Ham[k][N + i][j] += np.exp(- 1j * kvec.dot(-A2)) * t

def b(a1, a2): #RL vector finder
    #Tings we need
    omega = np.linalg.norm(np.cross(a1,a2))
    zhat = np.array([0,0,1])
    ##Obtaining RL vectors
    b1 = 2*np.pi*np.cross(a2,zhat)/omega
    b2 = 2*np.pi*np.cross(zhat,a1)/omega
    return b1,b2

A1p = A1
A2p = A2

B = b(A1p,A2p)

def symmpoints(b1, b2):

    ##Finding the special (high symmetry) points
    Gamma = np.array([0. , 0. , 0.])
    K = (2./3.)*b1 + (1./3.)*b2
    M = (1/2.)*b1

    return (Gamma, K, M)

def plotbands(theta , q , N , ons , lam):

    Ham = np.zeros( (5*q , 2*N , 2*N) , dtype = complex )

    Gamma = rotate( symmpoints(B[0],B[1])[0] , theta )
    K = rotate( symmpoints(B[0],B[1])[1] , theta )
    M = rotate( symmpoints(B[0],B[1])[2] , theta )

    # Gamma to K
    for r in range(2*q):
        FillKHam(Ham , r , r/(2*q-1)*(K-Gamma) , N , ons , lam)

    #K to M
    for r in range(q):
        FillKHam(Ham , 2*q + r , K + r/(q-1)*(M-K) , N , ons , lam)

    #M to Gamma
    for r in range(2*q):
        FillKHam(Ham , 3*q + r , M + r/(2*q-1)*(Gamma-M) , N , ons , lam)

    return np.linalg.eigh(Ham)[0]

H0 = plotbands(0 , q , N , ons , lam)

H120 = plotbands(60 , q , N , ons , lam)

H240 = plotbands(120 , q , N , ons , lam)

a = H0 - H120

b = H120 - H240

c = H240 - H0

cond1 = np.allclose( a , b , atol = 1e-3 )
cond2 = np.allclose( b , c , atol = 1e-3 )
cond3 = np.allclose( c , a , atol = 1e-3 )

print(cond1 & cond2 & cond3)

plt.plot(H0 , color = 'red')
plt.plot(H120 , color ='blue')
plt.plot(H240 , color = 'green')

plt.show()
