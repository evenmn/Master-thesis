from MonteCarlo import *

FermionsInHO = VMC(N=2, D=2, MC=100000, iterations=100, w=1, dx=0.1, eta=0.001)

FermionsInHO.Iterator()
