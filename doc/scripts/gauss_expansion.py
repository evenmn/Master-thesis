import numpy as np
import matplotlib.pyplot as plt

def STO(x, a=1):
    return np.exp(-a*np.abs(x))
    
def GTO(x, a=1):
    return np.exp(-a*x*x)
    
x = np.linspace(-3,3,1000)

# Expansion
def GTO_expansion(x, a, c):
    '''a and c must have the same lengths'''
    
    result = 0
    for i in range(len(a)):
        result += c[i]*GTO(x, a[i])
    
    return result
    
    
# Optimization
iterations = 10
basis_func = 3

a = np.random.random((basis_func,))
c = np.random.random((basis_func,))

#for iter in range(iterations):
    
    

plt.plot(x, STO(x))
plt.plot(x, GTO(x))
plt.show()
