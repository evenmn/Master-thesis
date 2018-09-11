import numpy as np

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


def matrix(Xa, N, D):
    '''N: Number of fully occupied shells'''
    i_list = list(N, D)
    length = len(i_list)
    print(length)
    
    A = np.ones([length, length])
    
    count = 0
    for i in range(length):
        for j in range(length):
            for k in range(len(i_list[0])):
                A[i,j] *= H(Xa[D*i+k], i_list[j,k])
                count += 1
                
    return A
    
Xa = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
print(matrix(Xa, 2, 3))
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
