import numpy as np
import matplotlib.pyplot as plt

plt.style.use("bmh")
plt.rcParams["font.family"] = "Serif"
plt.rc('figure', max_open_warning = 0)
ax = plt.gca()
ax.set_facecolor('white')

N = 1000

walks = np.zeros((N, 2))

for i in range(N-1):
    direction = np.random.randint(4)
    if direction == 0:
        move = [+1, 0]
    elif direction == 1:
        move = [-1, 0]
    elif direction == 2:
        move = [0, +1]
    elif direction == 3:
        move = [0, -1]
    walks[i+1] = walks[i] + move
    
plt.plot(walks[:,0], walks[:,1])
plt.xticks([])
plt.yticks([])
plt.axis('equal')
plt.grid()
plt.show()
