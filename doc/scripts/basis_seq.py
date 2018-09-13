import numpy as np
from scipy.linalg import lu

def list(N, D):
    '''Returns the index list used in Slater'''
    
    length = int(magic_numbers(N-1, D)/2)
    i_list = np.empty([length, D], dtype=int)
    counter = 0

    # Two dimensions
    if D == 2:
        for i in range(N):
            for s in range(i, N):
                j = s - i
                i_list[counter] = [i, j]
                counter += 1

        
    # Three dimensions
    elif D == 3:
        for i in range(N):
            for j in range(N):
                for s in range(i+j, N):
                    k = s - i - j
                    i_list[counter] = [i, j, k]
                    counter += 1
        
    sum = np.sum(i_list, axis=1)
    return i_list[np.argsort(sum)]
    
    
def factorial(n):
    '''Factorial'''
    if n == 0:
        return 1
    else:
        return factorial(n-1)*n
        
        
def binomial(n, p):
    '''Binomial coefficients'''
    return factorial(n+p)/(factorial(n)*factorial(p))
        
        
def magic_numbers(n, D, S=2):
    '''Magic numbers'''
    return int(S*binomial(n, D))


def H(x, n):
    '''Hermite polynomial of n'th degree'''

    if n == 0:
        return 1
    
    elif n == 1:
        return 2*x
    
    else:
        return 2*x*H(x,n-1)-2*(n-1)*H(x,n-2)


def dH(x, n):
    '''Derivative of Hermite polynomial of n'th degree'''
    if n == 0:
        return 0
    else:
        return 2*n*H(x,n-1)


def matrix(Xa, N, D):
    '''N: Number of fully occupied shells'''
    i_list = list(N, D)
    length = len(i_list)

    A = np.ones([length, length])
    
    for i in range(length):
        for j in range(length):
            for k in range(D):
                A[i,j] *= H(Xa[D*i+k], i_list[j,k])          
              
    return A
    
    
def derivative(Xa, N, D, k):
    '''Derivative of A matrix'''
    i_list = list(N, D)
    length = len(i_list)

    dA = np.zeros([length, length])
    
    # Find relevant row
    row = int(k/D)
    
    # Find indices of relevant row
    a = np.zeros(D, dtype=int)
    l = k%D
    for i in range(D):
        a[i] = k-l+i
    
    # Find matrix
    for i in range(length):
        dA[row, i] = dH(Xa[k], i_list[i, l])
        for j in range(D):
            if a[j] != k:
                dA[row, i] *= H(Xa[a[j]], i_list[i, j])
                
    return dA
    
    
def energy(Xa, N, D, k):
    A = matrix(Xa, N, D)
    dA = derivative(Xa, N, D, k)
    
    return np.trace(np.linalg.inv(A).dot(dA))
    
    
Xa = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
#print(matrix(Xa, 2, 3))
                
#print(derivative(Xa, 2, 3, 5))

print(energy(Xa, 2, 2, 2))   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
