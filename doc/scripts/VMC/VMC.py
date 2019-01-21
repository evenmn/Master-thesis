import numpy as np
from WaveFunction import *
from LocalEnergy import *
from Optimization import *

class VMC:
    def __init__(self, N, D, MC, iterations, w, dx, eta):
        self.N          = N
        self.D          = D
        self.MC         = MC
        self.iterations = iterations
        self.w          = w
        self.dx         = dx
        self.eta        = eta
        
        self.x = np.random.rand(N, D) - 0.5      # Initialize positions
        self.R = self.Diff(self.x)
        self.r = self.Dist(self.x)
        
        self.a = 0.5                               # Initialize gaussian parameter
        self.b = 1                               # Initialize Pade-Jastrow parameter

        
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
        
    def Iterator(self):
        ''' Parameter update '''
        for iter in range(self.iterations):
            # Declare objects
            EL  = LocalEnergy(self.a, self.b, self.N, self.D, self.w)
            Pot = Potential(self.a, self.b, self.N, self.D, self.w)
            Psi = WaveFunction(self.a, self.b, self.N, self.D, self.MC, self.iterations, self.w, self.dx, self.eta)
            GD  = Optimization(self.N, self.D, self.MC, self.iterations, self.w, self.dx, self.eta)
            
            Energy = 0
            for i in range(self.MC):
                xNew = self.x.copy()
                nRand = np.random.randint(self.N)
                dRand = np.random.randint(self.D)
                
                xNew[nRand, dRand] = self.x[nRand, dRand] + (np.random.rand(1)-0.5) * self.dx
                rNew = self.Dist(xNew)
                RNew = self.Diff(xNew)
                
                PsiRatio = Psi.Gauss(rNew)/Psi.Gauss(self.r)
                
                if(PsiRatio >= np.random.rand(1)):
                    self.x = xNew
                    self.r = rNew
                    self.R = RNew
                    
                Energy += EL.Laplacian(self.r) + Pot.HarmonicOscillator(self.r) #+ EL.Interaction(self.R)
            
            self.a = GD.GradientDescent(self.a, self.r)
            
            print("\n--- Iteration {} ---".format(iter+1))
            print("<E>: ", Energy/self.MC)
            print("a:   ", self.a)
                
        
        
