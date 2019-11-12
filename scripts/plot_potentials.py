import numpy as np
import matplotlib.pyplot as plt

def doubleWell(x, b, n):
    return (abs(x)**n - (b/2)**n)**2
    
x = np.linspace(-3,3,1000)

plt.plot(x, 2*doubleWell(x, 2, 1), label="$n=1$")
plt.plot(x, doubleWell(x, 2, 2), label="$n=2$")
plt.plot(x, doubleWell(x, 2, 3), label="$n=3$")
plt.axis([-3.2,3.2,0,8])
plt.xlabel("Radial distance")
plt.ylabel("Energy")
plt.legend(loc='upper center')
plt.grid()
#plt.show()

import matplotlib2tikz
matplotlib2tikz.save("/home/evenmn/Master-thesis/doc/text/pgf/doublewell.tex")
