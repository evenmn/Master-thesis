import numpy as np
from LocalEnergy import *
from Optimization import *
from Metropolis import *

#import psyco; psyco.full()

class VMC():
    def __init__(self, N, D, MC, iterations, w, dx, eta):
        self.N          = N
        self.D          = D
        self.MC         = MC
        self.iterations = iterations
        self.w          = w
        self.dx         = dx
        self.eta        = eta
        
        self.x = np.random.rand(N, D) - 0.5      # Initialize positions
        self.R = np.zeros((N, N))                # Initialize rij
        self.r = np.zeros(N)                     # Initialize ri
        
        self.a = 0.5                             # Initialize gaussian parameter
        self.b = 1                               # Initialize Pade-Jastrow parameter
        
        
    def Iterator(self):
        ''' Parameter update '''
        for iter in range(self.iterations):
            # Declare objects
            EL  = LocalEnergy(self.a, self.b, self.N, self.D, self.w)
            Pot = Potential(self.a, self.b, self.N, self.D, self.w)
            Int = Interaction(self.a, self.b, self.N, self.D, self.w)
            Kin = Kinetic(self.a, self.b, self.N, self.D, self.w)
            Met = Metropolis(self.a, self.b, self.N, self.D, self.w, self.x, self.r, self.dx)
            GD  = Optimization(self.N, self.D, self.w, self.eta)
            
            Energy = 0
            for i in range(self.MC):
                nRand = np.random.randint(self.N)                   # Next particle to move
                dRand = np.random.randint(self.D)                   # Direction to move in
                
                xNew, rNew, RNew, PsiRatio = Met(nRand, dRand)      # Metropolis
                
                if(PsiRatio >= np.random.rand(1)):
                    self.x = xNew
                    self.r = rNew
                    self.R = RNew
                    
                Energy += Kin(self.r) + Pot(self.r) + Int(self.R)
            
            #self.a = GD(self.a, self.r)
            
            print("\n--- Iteration {} ---".format(iter+1))
            print("<E>: ", Energy/self.MC)
            print("a:   ", self.a)
                
        
        
