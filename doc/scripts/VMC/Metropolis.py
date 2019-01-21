import numpy as np
from WaveFunction import *

class Metropolis:
    def __init__(self, a, b, N, D, w, x, r, dx):
        self.a = a
        self.b = b
        self.N = N
        self.D = D
        self.w = w
        self.x = x
        self.r = r
        self.dx = dx 
        self.Psi = WaveFunction(a, b, N, D, w)
        
    def __call__(self, nRand, dRand):
        sampling = 0
        if sampling == 0:
            return self.BruteForce(nRand, dRand)
        elif sampling == 1:
            return 

    def BruteForce(self, nRand, dRand):
        '''Brute Force Metropolis algorithm'''
        xNew = self.x.copy()
        xNew[nRand, dRand] = self.x[nRand, dRand] + (np.random.rand(1)-0.5) * self.dx
        rNew = self.Dist(xNew)
        RNew = self.Diff(xNew)
    
        PsiRatio = self.Psi.Gauss(rNew)/self.Psi.Gauss(self.r)
        return xNew, rNew, RNew, PsiRatio
        
    def ImportanceSampling(self, nRand, dRand):
        '''Importance Sampling algorithm'''
        xNew = self.x.copy()
        rNew = self.Dist(xNew)
        RNew = self.Diff(xNew)
        
        PsiRatio = self.Psi.Gauss(rNew)/self.Psi.Gauss(self.r)
        return xNew, rNew, RNew, PsiRatio
        
    def Diff(self, x):
        '''Calculate distance matrix between particles'''
        R = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(i):
                count = 0
                for d in range(self.D):
                    count += (x[i,d] - x[j,d])**2
                R[i,j] = np.sqrt(count)
        return R
        
    def Dist(self, x):
        '''Calculate distance from origin'''
        r = np.zeros(self.N)
        for i in range(self.N):
            DistPerParticle = 0
            for d in range(self.D):
                DistPerParticle += x[i,d]**2
            r[i] = DistPerParticle
        return r
