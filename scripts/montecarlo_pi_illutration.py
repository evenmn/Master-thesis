import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.style.use("bmh")
plt.rcParams["font.family"] = "Serif"
ax = plt.gca()
ax.set_facecolor('white')

N = 400

rect = patches.Rectangle((-1,-1), 2, 2, linewidth=2, edgecolor='k', facecolor='none')
circ = patches.Circle((0,0), 1, linewidth=2, edgecolor='k', facecolor='none')
ax.add_patch(rect)
ax.add_patch(circ)

for i in range(N):
    x = 2 * np.random.random() - 1
    y = 2 * np.random.random() - 1
    
    if np.linalg.norm([x,y]) <= 1:
        color = (0.203921568627451,0.541176470588235,0.741176470588235)
    else:
        color = (0.650980392156863,0.0235294117647059,0.156862745098039)
    
    plt.plot(x, y, 'o', markersize=3, color=color)

y = 1.07
space = 0.2

s = space/2
plt.arrow(-s, y, -1+s, 0, length_includes_head=True, width=0.01, head_width=0.04, head_length=0.04, color='k')
plt.arrow(+s, y, +1-s, 0, length_includes_head=True, width=0.01, head_width=0.04, head_length=0.04, color='k')
plt.text(-0.065, y-0.025, "$2r$", fontsize=16)

plt.axis('equal')
plt.ylim((-1.1, 1.1))
plt.axis('off')
plt.grid()
plt.show()
    
