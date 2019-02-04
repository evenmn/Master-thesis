import numpy as np
import matplotlib.pyplot as plt

beta12 = 1
gamma = 1

r12 = np.linspace(0.01, 1, 100)

def J(r12):
    return np.exp((beta12 * r12)/(1 + gamma * r12))
    
plt.plot(r12, J(r12))
plt.show()
