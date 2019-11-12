import numpy as np
import matplotlib.pyplot as plt

def Hermite(x, n, omegaSqrd = 1):
    if n == 0:
        return 1
    elif n == 1:
        return 2 * omegaSqrd * x
    else:
        return 2 * (omegaSqrd * x * Hermite(x, n-1) - (n-1) * Hermite(x,n-2));

def HermiteFunction(x, n, omega = 1):
    return Hermite(x, n, omega * omega) * np.exp(-0.5 * omega * x * x)

def Potential(x):
    return x*x
    
    
N = 1000
label_size = {'size':'14'}
x = np.linspace(-4,4,N)
y = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
plt.yticks(np.arange(y.min(), y.max(), 1))


plt.plot(x, 0.5*np.ones(N), '--k')
plt.plot(x, HermiteFunction(x, 0)*0.4 + 0.5, '-r', label="n=1")
plt.plot(x, 1.5*np.ones(N), '--k')
plt.plot(x, HermiteFunction(x, 1)*0.4 + 1.5, '-g', label="n=2")
plt.plot(x, 2.5*np.ones(N), '--k')
plt.plot(x, HermiteFunction(x, 2)*0.2 + 2.5, '-y', label="n=3")
plt.plot(x, 3.5*np.ones(N), '--k')
plt.plot(x, HermiteFunction(x, 3)*0.08 + 3.5, '-c', label="n=4")
plt.plot(x, Potential(x), '-b', linewidth=2)
plt.text(-4, 0.6, "n=1", **label_size)
plt.text(-4, 1.6, "n=2", **label_size)
plt.text(-4, 2.6, "n=3", **label_size)
plt.text(-4, 3.6, "n=4", **label_size)
plt.axis([-4.5,4.5,-0.2,5])
plt.grid()
plt.xlabel("Relative distance from center, [r$_0$]",**label_size)
plt.ylabel("Energy in natural units",**label_size)
plt.show()
