import numpy as np
from scipy.special import binom
import binfuncts as bf
import time
from params import L , N , Nsize , Np1size , Nm1size
import itertools as it

BasisN = bf.makebasis(L , N)
BasisNp1 = bf.makebasis(L , N+1)
BasisNm1 = bf.makebasis(L , N-1)


def binvec(arr , m): #This returns an ARRAY
    to_str_func = np.vectorize(lambda x: np.binary_repr(x).zfill(m))
    strs = to_str_func(arr)
    ret = np.zeros(list(arr.shape) + [m], dtype=np.int8)
    for bit_ix in range(0, m):
        fetch_bit_func = np.vectorize(lambda x: x[bit_ix] == '1')
        ret[...,bit_ix] = fetch_bit_func(strs).astype("int8")
    return ret

def fillcounter(s , n , L):
    return np.sum(binvec(s & sum(2**(L-i-1) for i in range(n)) , L)) #Sums each element

def cdagmat(alpha):
    cmat = np.zeros((Np1size,Nsize) , dtype=int)
    cond1 = BasisN & 2**(L-alpha-1) == 0
    indices = np.nonzero(cond1)[0]
    bc1 = BasisN[cond1]
    s1 = bf.flip(bc1 , alpha , L )
    i1 = np.searchsorted(BasisNp1 , s1)
    for k in range(len(indices)):
        cmat[ i1[k] ][ indices[k] ] += (-1)**bf.fillcounter(BasisN[indices[k]] , alpha , L)
    return cmat

def cdagmat2(alpha):
    cmat = np.zeros((Nsize,Nm1size) , dtype=int)
    cond1 = BasisNm1 & 2**(L-alpha-1) == 0
    indices = np.nonzero(cond1)[0]
    bc1 = BasisNm1[cond1]
    s1 = bf.flip(bc1 , alpha , L )
    i1 = np.searchsorted(BasisN , s1)
    for k in range(len(indices)):
        cmat[ i1[k] ][ indices[k] ] += (-1)**bf.fillcounter(BasisNm1[indices[k]] , alpha , L)
    return cmat

def cmat(beta):
    cmat = np.zeros((Nm1size,Nsize) , dtype=int)
    cond1 = BasisN & 2**(L-beta-1) != 0
    indices = np.nonzero(cond1)[0]
    bc1 = BasisN[cond1]
    s1 = bf.flip(bc1 , beta , L )
    i1 = np.searchsorted(BasisNm1 , s1)
    for k in range(len(indices)):
        cmat[ i1[k] ][ indices[k] ] = (-1)**bf.fillcounter(BasisN[indices[k]] , beta , L)
    return cmat

def cmat2(beta):
    cmat = np.zeros((Nsize,Np1size) , dtype=int)
    cond1 = BasisNp1 & 2**(L-beta-1) != 0
    indices = np.nonzero(cond1)[0]
    bc1 = BasisNp1[cond1]
    s1 = bf.flip(bc1 , beta , L )
    i1 = np.searchsorted(BasisN , s1)
    for k in range(len(indices)):
        cmat[ i1[k] ][ indices[k] ] = (-1)**bf.fillcounter(BasisNp1[indices[k]] , beta , L)
    return cmat

# Now we want to store the matrices for all of the alphas and betas

cdags = []

for i in range(L):
    cdags.append(cdagmat(i))

cdags2 = []

for i in range(L):
    cdags2.append(cdagmat2(i))

cs = []

for j in range(L):
    cs.append(cmat(j))

cs2 = []

for j in range(L):
    cs2.append(cmat2(j))

def mateld8d6(t , b1 , b2 , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , cdags , cs):

    ans = 0

    for alpha , beta in it.product(range(0,L) , repeat = 2):
        ans += (-1)**N * t[alpha][beta] * np.einsum('i,ij,j->' , np.conjugate(SpecNp1[1][:,e1]) , cdags[alpha] , SpecN[1][:,b1] ) * np.einsum('i,ij,j->' , np.conjugate(SpecNm1[1][:,e2]) , cs[beta] , SpecN[1][:,b2] )

    return ans

def fullmateld8d6(t , b1 , b2 , b1p , b2p , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , E0 , cdags , cs):

    Ee = SpecNp1[0][e1] + SpecNm1[0][e2]

    return ( np.conjugate(mateld8d6(t , b1p , b2p , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , cdags , cs)) * mateld8d6(t , b1 , b2 , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , cdags , cs) )/(E0-Ee)

def mateld6d8(t , b1 , b2 , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , cdags , cs):

    ans = 0

    for alpha , beta in it.product(range(0,L) , repeat = 2):
        ans += t[alpha][beta] * np.einsum( 'i,ij,j->' , np.conjugate(SpecNm1[1][:,e1]) , cs[alpha] , SpecN[1][:,b1] ) * np.einsum( 'i,ij,j->' , np.conjugate(SpecNp1[1][:,e2]) , cdags[beta] , SpecN[1][:,b2] )

    return ans

def fullmateld6d8(t , b1 , b2 , b1p , b2p , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , E0 , cdags , cs):

    Ee = SpecNm1[0][e1] + SpecNp1[0][e2]

    return ( np.conjugate(mateld6d8(t , b1p , b2p , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , cdags , cs)) * mateld6d8(t , b1 , b2 , e1 , e2 , SpecN , SpecNp1 , SpecNm1 , cdags , cs) )/(E0-Ee)
