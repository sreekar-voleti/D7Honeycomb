from scipy.special import binom

#Parameters (All in eV)
Vcf = 0. # 2 eV
lam = 0.08 #0.03 +/- 10%
Jh = 0.7 # 0.7 +/- 10%
Jp = Jh
B = 0.00001 # Tiny magnetic field

dist = 0.65

L = 10
N = 7

Nm1 = N-1
Np1 = N+1

Nsize = int(binom(L , N))
Np1size = int(binom(L , N+1))
Nm1size = int(binom(L , N-1))

# Angles that determine the C3 axis
theta = 1.5707406875049383
phi = 0.6616644973460342
