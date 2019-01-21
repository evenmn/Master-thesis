import numpy as np
#from MonteCarlo import VMC

class LocalEnergy:
    def __init__(self, a, b, N, D, w):
        '''Constructor'''
        self.a = a
        self.b = b
        self.N = N
        self.D = D
        self.w = w
    
    
class Kinetic(LocalEnergy):
    def __init__(self, a, b, N, D, w):
        '''Constructor'''
        LocalEnergy.__init__(self, a, b, N, D, w)
        
    def __call__(self, r):
        system = 0
        if system == 0:
            return self.Laplacian(r)
            
    def Laplacian(self, r):
        '''Laplace operator'''
        result = 0
        result += 4 * self.a * self.a * np.sum(np.square(r))
        result -= 2 * self.D * self.a * self.N
        
        return -0.5 * result
        
        
class Potential(LocalEnergy):
    def __init__(self, a, b, N, D, w):
        LocalEnergy.__init__(self, a, b, N, D, w)
        
    def __call__(self, r):
        system = 0
        if system == 0:
            return self.HarmonicOscillator(r)  
        elif system == 1:
            return self.AtomicNucleus(r)

    def HarmonicOscillator(self, r):
        '''Harmonic oscillator potential'''
        return 0.5 * self.w * self.w * np.sum(np.square(r))
        
    def AtomicNucleus(self, r):
        '''Atomic potential'''
        return 0.5 * self.w * self.w * np.reciprocal(R).sum()
    
    
class Interaction(LocalEnergy):
    def __init__(self, a, b, N, D, w):
        LocalEnergy.__init__(self, a, b, N, D, w)
        
    def __call__(self, R):
        inter = 0
        if inter == 0:
            return 0
        if inter == 1:
            return self.Coulomb(R)
            
    def Coulomb(self, R):
        '''Interaction energy'''
        summer = 0
        for i in range(self.N):
            for j in range(i):
                summer += 1/R[i,j]
        return summer
