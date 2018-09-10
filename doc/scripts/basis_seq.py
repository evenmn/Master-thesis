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
    '''Binomial'''
    return factorial(n+p)/(factorial(n)*factorial(p))
        
        
def magic_numbers(n, D, S=2):
    '''Magic numbers'''
    return int(S*binomial(n, D))


print(list(4, 3))
