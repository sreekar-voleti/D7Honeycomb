import sys
import numpy as np
from scipy.special import binom
from params import dist , theta , phi
from Amatrices import A_Co1_new

#change
#change2

#The goal in here is to create a function that creates a basis for ANY filling of d orbitals

# Function to convert a number to a set of 10 bits

def binstate(s , L): #This returns a STRING
    return np.binary_repr(s , width = L)

# To convert a number from binary to int, just go for int(bin , 2)

# Function to determine if a binary string can correspond to a valid d3 states

def isvalid(s , L , Ne): # s is the number, Ne is the number of electrons
    proxy = binstate(s , L)
    val = 0
    for i in range(L):
        val += int(proxy[i])
    if val == Ne:
        return 1
    else:
        return 0

# Function that generates a list of all the numbers which can be used for constructing a basis

def makebasis(L , Ne): # ***MAKE SURE*** Ne <= L otherwise you'll end up in an infinite loop I think

    counter = int(binom(L , Ne)) - 1 # Tells you how many basis elements to expect

    basis = np.array([] , dtype = int) # Where I store the basis number

    i = 0

    while(counter > -1): # I know this is sloppy
        if isvalid(i , L , Ne) == 1 :
            basis = np.append(basis , i)
            counter -= 1 # To ensure the loop breaks once we hit the desired number of basis states.
        i += 1

    return basis

# So now we have a function that will create a basis for any electron filling.

# Function that determines the orbital occupations of a state given the integer associated with the basis element
# This is useful for the intraorbital Hubbard term and Crystal Field

def orbitaloccupation(s , L):
    proxy = binstate(s , L)
    c = []
    for i in range(int(L/2)):
        c.append( int(proxy[2*i]) + int(proxy[2*i+1]) )
    return c

def crystalfield(s , L):
    c = orbitaloccupation(s , L)
    return c[3]+c[4] #Just counts the number of eg electrons

def tinymagfield(s , L):
    proxy = binstate(s , L)
    Nup = 0
    Ndown = 0
    for i in range(int(L/2)):
        Nup += int(proxy[2*i]) # Even sites are associated with spin up
        Ndown += int(proxy[2*i+1]) # Odd sites with spin down
    return (Nup - Ndown)

# Function for intraorbital Hubbard

def intrahubbard(s , L):
    bruh = 0
    for i in range(int(L/2)):
        if orbitaloccupation(s , L)[i] == 2:
            bruh += 1
    return bruh

# Function for interorbital Hubbard

def interhubbard(s , L): #This is a diagonal term in the Hamiltonian, so we are just looking for the product of occupation numbers.
    proxy = orbitaloccupation(s , L)
    mandem = 0
    for i in range(int(L/2)):
        for j in range(int(L/2)):
            if i < j:
                mandem += proxy[i]*proxy[j]
    return mandem

# Function to determine how many PAIRS of parallel spins there are in a given state

def parallelspins(s , L):
    proxy = binstate(s , L)
    Nup = 0
    Ndown = 0
    for i in range(int(L/2)):
        Nup += int(proxy[2*i]) # Even sites are associated with spin up
        Ndown += int(proxy[2*i+1]) # Odd sites with spin down
    return int(binom(Nup,2) + binom(Ndown,2))

# The idea is something like the following:
## There are 10 possible states __ __  __ __  __ __  __ __  __ __
##                              u  d   u  d   u  d   u  d   u  d
##                                1      2      3      4      5
## A spin flip operator flips a spin, so it takes a state in an odd(even) indexed location to an even(odd) one
## The flip operator has 2 orbitals associated with it (alpha and beta). Beta acts first (see Hamiltonian)
## Things to consider : How many filled states (1s in the string) does that spin up annihilation operator need to go through to reach the one it is associated with
### This number is the power to which we must raise (-1) in order to evaluate this matrix element

# Recall the fermion anticommutation relations ({cd_m , c_n} = delta_mn , all others anticommute)

# A function which, given a particular "location" in the binary string, determines how many filled states come "before" it

def fillcounter(s , n , L): # s is the state, n is the index
    proxy = binstate(s , L)
    counter = 0
    for i in range(n): # Loops from 0 to n-1
        if proxy[i] == '1':
            counter += 1
    return counter

# Function for spin flip interaction we need to be able to change bits.

def flip(s , n , L): #Flips the occupation at location n in the bit string
    return s ^ (2 ** (L - n - 1))

# Function that finds the index of the new basis elements. We use the fact that the basis is ordered.

# A function that tells you if the spin flip term acts on it, and if it does, what sign it would come with

def hundflipspin(s , a , b , L , basis , Ntot):
    sav = s # Saving the initial state
    p = binstate(s , L)
    count = 0 # Counting the exponent of (-1)
    if ( int(p[2*a]) + int(p[2*b+1]) == 2 ) and ( int(p[2*a+1]) + int(p[2*b]) == 0 ) :
        #We check the fillcounter BEFORE acting on the state with the creation / annihilation operator
        count += fillcounter(s , 2*b , L)
        s = flip(s , 2*b , L)
        count += fillcounter(s , 2*b+1 , L)
        s = flip(s , 2*b+1 , L)
        count += fillcounter(s , 2*a+1 , L)
        s = flip(s , 2*a+1 , L)
        count += fillcounter(s , 2*a , L)
        s = flip(s , 2*a , L)
    if s != sav:
        return np.searchsorted(basis , s) , (-1)**count
    else:
        return 0 , 0 # This ensures that if the spin flip term doesn't act on the state, then it won't be added to the Hamiltonian (second term is zero)

# A function to find the pair hopping

def pairhop(s , a , b , L , basis , Ntot):
    sav = s
    c = orbitaloccupation(s , L)
    count = 0
    if c[b] == 2 and c[a] == 0:
        count += fillcounter(s , 2*b , L)
        s = flip(s , 2*b , L)
        count += fillcounter(s , 2*b+1 , L)
        s = flip(s , 2*b+1 , L)
        count += fillcounter(s , 2*a+1 , L)
        s = flip(s , 2*a+1 , L)
        count += fillcounter(s , 2*a , L)
        s = flip(s , 2*a , L)
    if s != sav:
        return np.searchsorted(basis , s) , (-1)**count
    else:
        return 0 , 0 # This ensures that if the spin flip term doesn't act on the state, then it won't be added to the Hamiltonian (second term is zero)

## OK COOL, now on to SOC

# This part is only for d orbitals. It seems like a useful move would be to convert the LS matrices into 10x10 matrices

##Spin Matrices

Sx = np.array([[0,1],[1,0]])
Sy = np.array([[0,complex(0,-1)],[complex(0,1),0]])
Sz = np.array([[1,0],[0,-1]])

#L matrices (full)

Lx = np.zeros((5,5), dtype = complex)

Lx[1][2] += 1.j
Lx[2][1] += -1.j
Lx[0][3] += -1.j
Lx[3][0] += 1.j
Lx[0][4] += -complex(0 , np.sqrt(3.))
Lx[4][0] += complex(0,np.sqrt(3.))

Ly = np.zeros((5,5), dtype = complex)

Ly[0][2] += -1.j
Ly[2][0] += 1.j
Ly[1][3] += -1.j
Ly[3][1] += 1.j
Ly[1][4] += complex(0,np.sqrt(3.))
Ly[4][1] += -complex(0,np.sqrt(3.))

Lz = np.zeros((5,5), dtype = complex)

Lz[0][1] += 1.j
Lz[1][0] += -1.j
Lz[3][2] += -2.j
Lz[2][3] += 2.j

#L dot n matrix

Ldotn = np.sin(theta)*np.cos(phi)*Lx + np.sin(theta)*np.sin(phi)*Ly + np.cos(theta)*Lz

Lmat = np.kron(Ldotn , np.identity(2))

# We want to match the orbital + spin structure of the basis.

def LS(i):
    if i == 0:
        return np.kron(Lx , Sx)
    elif i == 1:
        return np.kron(Ly , Sy)
    elif i == 2:
        return np.kron(Lz , Sz)

def soc(s , a , b , basis , L):

    proxy = binstate(s , L)
    count = 0
    matel = 0

    if (( proxy[b] == '1' ) and ( proxy[a] == '0' )) :

        count += fillcounter(s , b , L)
        s = flip(s , b , L)
        count += fillcounter(s , a , L)
        s = flip(s , a , L)

        for i in range(3):
            matel += LS(i)[a][b]

    return np.searchsorted(basis , s) , (-1)**count * 0.5 * matel

def LS2(i , Unitary):
    if i == 0:
        return np.kron(-np.conjugate(Unitary.T).dot(Lx).dot(Unitary) , Sx)
    elif i == 1:
        return np.kron(np.conjugate(Unitary.T).dot(Ly).dot(Unitary) , Sy)
    elif i == 2:
        return np.kron(-np.conjugate(Unitary.T).dot(Lz).dot(Unitary) , Sz)

def soc2(s , a , b , basis , L , Unitary):

    proxy = binstate(s , L)
    count = 0
    matel = 0.

    if (( proxy[b] == '1' ) and ( proxy[a] == '0' )) :

        count += fillcounter(s , b , L)
        s = flip(s , b , L)
        count += fillcounter(s , a , L)
        s = flip(s , a , L)

        for i in range(3):
            matel += LS2(i , Unitary)[a][b]

    return np.searchsorted(basis , s) , (-1)**count * 0.5 * matel


#Asp = np.kron(dist*A + CFMat , np.identity(2))

def ons(s , a , b , basis , L , Asp):

    proxy = binstate(s, L)
    count = 0

    if a != b:
        if ((proxy[b] == '1') and (proxy[a] == '0')) :

            count += fillcounter(s , b , L)
            s = flip(s , b , L)
            count += fillcounter(s , a , L)
            s = flip(s , a , L)

            return np.searchsorted(basis , s) , (-1)**count * Asp[a][b]

        else:
            return 0 , 0
    elif a == b:
        if proxy[b] == '1':
            count += fillcounter(s , b , L)
            s = flip(s , b , L)
            count += fillcounter(s , a , L)
            s = flip(s , a , L)

            return np.searchsorted(basis , s) , (-1)**count * Asp[a][b]

        else:
            return 0 , 0

# Function to solve for a spectrum of filling Ne

def SolveSpectrumt2g(L , Ne , basis , Vcf , lam , U , Jh , Jp , Up , B):

    Nbasis = int(binom(L , Ne))

    # Defining the Hamiltonian

    Ham = np.zeros( (Nbasis , Nbasis) , dtype = complex)

    for i in range(Nbasis):
        #Hubbard U
        Ham[i][i] += U * intrahubbard(basis[i] , L)
        #Hubbard Up
        Ham[i][i] += Up * interhubbard(basis[i] , L)
        #Parallel spins Hund
        Ham[i][i] -= Jh * parallelspins(basis[i] , L)
        #Tiny magnetic field
        Ham[i][i] -= B * tinymagfield(basis[i] , L)
        for a in range(L):
            for b in range(L):
                #SOC
                if a != b: # Since no off diagonals for SOX
                    so = soc(basis[i] , a , b , basis , L, Unitary)
                    Ham[so[0]][i] += lam * so[1]

                #Off-diagonal Interactions
                if a < int(L/2) and b < int(L/2):
                    if a != b:
                        sf = hundflipspin(basis[i] , a , b , L , basis , Nbasis)
                        Ham[sf[0]][i] -= Jh * sf[1]
                        ph = pairhop(basis[i] , a , b , L , basis , Nbasis)
                        Ham[ph[0]][i] += Jp * ph[1]

    return np.linalg.eigh(Ham) #Returns both the eigenenvalues and the eigenvectors

#We take the Co1 atom to be the central one, and flip signs on Co2

def SolveSpectrumfull(L , Ne , basis , lam , U , Jh , Jp , Up , B):

    Asp = np.kron( A_Co1_new + B*Ldotn , np.identity(2) )

    Nbasis = int(binom(L , Ne))

    # Defining the Hamiltonian

    Ham = np.zeros( (Nbasis , Nbasis) , dtype = complex)

    for i in range(Nbasis):
        #Hubbard U
        Ham[i][i] += U * intrahubbard(basis[i] , L)
        #Hubbard Up
        Ham[i][i] += Up * interhubbard(basis[i] , L)
        #Parallel spins Hund
        Ham[i][i] -= Jh * parallelspins(basis[i] , L)
        #Tiny magnetic field
        #Ham[i][i] -= B * tinymagfield(basis[i] , L)
        for a in range(L):
            for b in range(L):
                boop = ons(basis[i] , a , b , basis , L , Asp)
                Ham[boop[0]][i] += boop[1]
                #SOC
                if a != b: # Since no off diagonals for SOX
                    so = soc(basis[i] , a , b , basis , L )
                    Ham[so[0]][i] += lam * so[1]

                #Off-diagonal Interactions
                if a < int(L/2) and b < int(L/2):
                    if a != b:
                        sf = hundflipspin(basis[i] , a , b , L , basis , Nbasis)
                        Ham[sf[0]][i] -= Jh * sf[1] #Hund spin flip
                        ph = pairhop(basis[i] , a , b , L , basis , Nbasis)
                        Ham[ph[0]][i] += Jp * ph[1] #Pair hopping

    return np.linalg.eigh(Ham) #Returns both the eigenenvalues and the eigenvectors
