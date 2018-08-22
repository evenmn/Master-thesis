import numpy as np
import numpy.linalg as LA
from time import clock
from scipy.special import binom


def H(x, n):
    '''Hermite polynomial of n'th degree'''
    
    if n==0:
        return 1
    elif n==1:
        return 2*x
    else:
        return 2*x*H(x,n-1)-2*(n-1)*H(x,n-2)


def NQS_WF(Xa, v, sigma):
    '''Neural Quantum State Wavefunction (NQS-WF)'''
    
    prod = 1
    for i in range(len(v)):
        prod *= (1 + np.exp(-v[i]))

    return np.exp(-np.dot(Xa, Xa)/(2 * sigma * sigma)) * prod

    
def magic_numbers(n, D, S=2):
    '''Magic numbers of a D-dimensional quantum dot.
    Returns total number of particles, N, in n fully 
    occupied orbitals. Starts at 0, following the 
    arithmetic series 
    
        N(n) = S * binomial(n+D, D)
        
    with S as the number of spin configurations.'''
    
    if n < 0:
        return 0
    else:
        return S * binom(n+D, D)
    
    
def magic_numbers_inverse(sum, D):
    '''Given a magic number, 'sum', this function
    returns number of orbitals. To be used when
    checking if orbitals are fully occupied.'''
    
    if D == 2:
        return (-3 + np.sqrt(1 + 4*sum))/2
    elif D == 3:
        roots = np.roots([1, 6, 11, 6-3*sum])
        return np.real(magic_numbers_inverse(8)[-1])


def Slater(Xa, v, sigma, D=2):
    '''Setting up Slater determinant'''
    
    N = len(Xa)/D                            # Needs to be an even number
    
    n_orbitals = magic_numbers_inverse(N, D) # Number of orbitals given N
    
    # Check if the orbitals are fully occupied, otherwise break
    if abs(n_orbitals - int(n_orbitals)) > 0.01:
        raise ValueError("Number of particles needs to be a magic number")
    
    levels = np.linspace(0, n_orbitals, n_orbitals + 1, dtype=int)
    
    D_up = np.array([[H(1,0), H(Xa[0], 1), H(Xa[1], 1)], \
                     [H(1,0), H(Xa[4], 1), H(Xa[5], 1)], \
                     [H(1,0), H(Xa[8], 1), H(Xa[9], 1)]])
    
    print(D_up)
    
    D_dn = np.array([[H(1,0), H(Xa[2], 1), H(Xa[3], 1)], \
                     [H(1,0), H(Xa[6], 1), H(Xa[7], 1)], \
                     [H(1,0), H(Xa[10], 1), H(Xa[11], 1)]])
                     
    print(D_dn)
    
    return LA.det(D_up)*LA.det(D_dn)*NQS_WF(Xa, v, sigma)


start = clock()
print(Slater(np.array([.75,.64,.83,.52,.31,.22,.33,.44,.15,.76,.27,.58]), np.array([1,1]), 1))
end = clock()

print(end-start)
