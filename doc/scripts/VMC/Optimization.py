import numpy as np

class Optimization:

    def __init__(self, N, D, MC, iterations, w, dx, eta):
        self.N = N
        self.D = D
        self.MC = MC
        self.iterations = iterations
        self.w = w
        self.dx = dx
        self.eta = eta

    def GradientDescent(self, a, r):
        return a - self.eta*(3*self.N - 4*a*a*np.sum(np.square(r)))
