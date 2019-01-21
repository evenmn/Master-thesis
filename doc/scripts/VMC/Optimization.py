import numpy as np

class Optimization:
    def __init__(self, N, D, w, eta):
        self.N = N
        self.D = D
        self.w = w
        self.eta = eta
        
    def __call__(self, a, r):
        return self.GradientDescent(a, r)

    def GradientDescent(self, a, r):
        return a - self.eta*(3*self.N - 4*a*a*np.sum(np.square(r)))
