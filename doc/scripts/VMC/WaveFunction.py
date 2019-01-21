import numpy as np

class WaveFunction:
    def __init__(self, a, b, N, D, MC, iterations, w, dx, eta):
        self.a = a
        self.b = b
        self.N = N
        self.D = D
        self.MC = MC
        self.iterations = iterations
        self.w = w
        self.dx = dx
        self.eta = eta

    def Gauss(self, r):
        '''Gaussian function'''
        return np.exp(-2*self.a*np.sum(np.square(r)))
        
    def PadeJastrow(self, R):
        '''Pade-Jastrow factor'''
        for i in range(self.N):
            for j in range(i):
                return R(i,j)/(1+self.b*R(i,j))
