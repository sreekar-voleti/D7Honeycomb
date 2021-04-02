import binfuncts as bf
import numpy as np
import multiprocessing as mp #For parallel programming
import DFT_OP_2nn as nnmat #TSD's matrices
import itertools as it #To generate the tuples needed for things
import time #Self explanatory
import sys #To take system argyments at the command line
import transition as tr
from transition import cdags , cs
from params import Vcf , lam , Jh , Jp , B , L , N , Np1 , Nm1

def get_result(result):
    global results
    results.append(result)

def testfn():
    return 1.

#We are varying U in this calculation between 2 and 5 eV.
## The argument will be an array from 0 - 10 ( so we multiply each number by 0.3 )

#First argument in the command line gives the neighbour (see dftmatrices.py to see which of the matrices correspond to which neighbour)

if sys.argv[1] == '1':
    t_12 = np.kron(nnmat.t1_11 , np.identity(2))
elif sys.argv[1] == '4':
    t_12 = np.kron(nnmat.t4_11 , np.identity(2))
elif sys.argv[1] == '2':
    t_12 = np.kron(nnmat.t2_11 , np.identity(2))
elif sys.argv[1] == '5':
    t_12 = np.kron(nnmat.t5_11 , np.identity(2))
elif sys.argv[1] == '3':
    t_12 = np.kron(nnmat.t3_11 , np.identity(2))
elif sys.argv[1] == '6':
    t_12 = np.kron(nnmat.t6_11 , np.identity(2))

U = 2.0 + int(sys.argv[2]) * 0.3
Up = U - 2*Jh

start = time.time()

# Bases for the 3 spaces we are looking at.
BasisN = bf.makebasis(L , N)
BasisNm1 = bf.makebasis(L , Nm1)
BasisNp1 = bf.makebasis(L , Np1)

#Eigensystems for the 3 systems (Distortion of 0.8 times the A matrix gives good agreement with data)
SpecN = bf.SolveSpectrumfull(L , N , BasisN , lam , U , Jh , Jp , Up , B)
SpecNm1 = bf.SolveSpectrumfull(L , Nm1 , BasisNm1 , lam , U , Jh , Jp , Up , B)
SpecNp1 = bf.SolveSpectrumfull(L , Np1 , BasisNp1 , lam , U , Jh , Jp , Up , B)

E0 = 2 * SpecN[0][0] # 2 Sites each at the ground state of the N electron problem

print( np.round( (SpecN[0][0:12] - SpecN[0][0])*1e3 ) )

Hameff = np.zeros((4,4) , dtype = complex)

counter1 = 0

for ma in it.product(range(0, 2), repeat=4): # 16 terms here

    st = time.time()
    pool = mp.Pool(mp.cpu_count())

    results = []

    for el in it.product(range(0,tr.Np1size) , range(0,tr.Nm1size)): # 9450 terms
        pool.apply_async( tr.fullmateld8d6 , args=(t_12 , ma[0] , ma[1] , ma[2] , ma[3] , el[0] , el[1] , SpecN , SpecNp1 , SpecNm1 , E0 , cdags , cs) , callback=get_result )

    for el in it.product(range(0,tr.Nm1size) , range(0,tr.Np1size)): # 9450 terms
        pool.apply_async( tr.fullmateld6d8 , args=(t_12 , ma[0] , ma[1] , ma[2] , ma[3] , el[0] , el[1] , SpecN , SpecNp1 , SpecNm1 , E0 , cdags , cs) , callback=get_result )

    pool.close()
    pool.join()

    ind1 = int(str(ma[0])+str(ma[1]) , 2)
    ind2 = int(str(ma[2])+str(ma[3]) , 2)

    Hameff[ind1][ind2] += sum(results)

    counter1 += 1
    print(counter1 , time.time()-st)

# print(np.round(Hameff,6) )
# print(np.linalg.eigh(Hameff)[0])

#Just add a thing where you create a file holding the array

np.save('Hameffs/OP_2nn_t'+sys.argv[1]+'_U='+str(U)+'eV.npy' , Hameff )
