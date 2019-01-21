import numpy as np

class LocalEnergy:
    def __init__(self, a, b, N, D, w):
        self.a = a
        self.b = b
        self.N = N
        self.D = D
        self.w = w
    
    def Laplacian(self, r):
        '''Laplace operator'''
        result = 0
        result += 4 * self.a * self.a * np.sum(np.square(r))
        result -= 2 * self.D * self.a * self.N
        
        return -0.5 * result
        
    def Interaction(self, R):
        '''Interaction energy'''
        return np.reciprocal(R).sum()
        
        
class Potential(LocalEnergy):
    def __init__(self, a, b, N, D, w):
        LocalEnergy.__init__(self, a, b, N, D, w)

    def HarmonicOscillator(self, r):
        '''Harmonic oscillator potential'''
        return 0.5 * self.w * self.w * np.sum(np.square(r))
        
    def AtomicNucleus(self, r):
        '''Atomic potential'''
        return 0.5 * self.w * self.w * np.reciprocal(R).sum()
