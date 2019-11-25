import numpy as np
import matplotlib.pyplot as plt

plt.style.use("bmh")
plt.rcParams["font.family"] = "Serif"
ax = plt.gca()
ax.set_facecolor('white')

def bias2(x):
    return 0.5*np.exp(x)
    
def variance(x):
    return -np.log(x+0.5)+2
    
def error(x):
    return bias2(x) + variance(x)

x = np.linspace(0, 2, 1000)

plt.text(1.7, 1.4, 'bias$^2$')
plt.text(0.2, 0.9, 'variance')
plt.text(1.2, 3.4, 'Eout')

plt.plot(x, bias2(x))
plt.plot(x, variance(x))
plt.plot(x, error(x))
plt.axvline(x[np.argmin(error(x))])

plt.xlabel("Model complexity")
plt.ylabel("Error")

plt.xticks([])
plt.yticks([])

#plt.show()

import tikzplotlib

tikzplotlib.save("tikz/biasvariance.tex")
