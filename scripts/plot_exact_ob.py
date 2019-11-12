import numpy as np
import matplotlib.pyplot as plt

def f(r):
    return r*r*np.exp(-2*r)
    
r = np.linspace(0, 5, 1000)

plt.plot(r, f(r))
plt.show()
