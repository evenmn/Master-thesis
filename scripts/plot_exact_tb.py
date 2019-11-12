import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.style.use("bmh")

def exact(r1, r2, w):
    return 2 * np.sqrt(w/np.pi) * np.exp(- w * (r1 * r1 + r2 * r2))
    
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
if __name__ == "__main__":
    N = 1000
    radius = 3
    r = np.linspace(-radius, radius, N)
    
    data = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            data[i, j] = exact(r[i], r[j], 1)
    
    #data /= np.sum(data)
    
    size = 28
    size_ticks = 20
    label_size = {"size":str(size)}
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams.update({'figure.autolayout': True})

    fig, ax = plt.subplots(figsize=(8,6))
    
    img = ax.imshow(data, cmap=plt.cm.jet, extent=[-radius,radius,-radius,radius])
    cbar = fig.colorbar(img, fraction=0.046, pad=0.04) #, format=ticker.FuncFormatter(fmt))
    cbar.set_label(r'$\rho(r_i,r_j)$', rotation=90, labelpad=10, y=0.5, **label_size)
    cbar.ax.tick_params(labelsize=size_ticks)
    
    plt.tight_layout()
    
    ax.set_xlabel("$r_j$", **label_size)
    ax.set_ylabel("$r_i$", **label_size)
    ax.tick_params(labelsize=size_ticks)
    
    tick = [-3, -2, -1, 0, 1, 2, 3]
    ax.set_xticks(tick)
    ax.set_yticks(tick)
    
    plt.grid()
    plt.show()
