import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def exact(r1):
    return np.exp(-r1**2)

# Numerical values
data1 = np.loadtxt("../data/ob_density.dat")
#data2 = np.loadtxt("../data/OB_w_interaction.dat")

label_size = {"size":"14"}


r1 = np.linspace(0, 5, len(data1))
#r2 = np.linspace(0, 3, len(data2))

data2 = exact(r1)
data1 = data1/np.sum(data1)
data2 = data2/np.sum(data2)

sns.set()
data1/=(r1)
plt.plot(r1, data1, '.', markersize=4, label="Interaction")
#plt.plot(r1, data2)

#plt.plot(r2, data2, '-', linewidth=1.0, label="W")
plt.title("6 particles with interaction")
plt.xlabel("r", **label_size)
plt.ylabel(r"$\rho$", **label_size)
#plt.legend(loc="best")
plt.savefig("../plots/ob.png")
plt.show()
