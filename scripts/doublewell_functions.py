import numpy as np
import scipy.special
from scipy.integrate import trapz, simps
from scipy.special import factorial
import matplotlib.pyplot as plt

def hermite(n,x,w):
    return scipy.special.hermite(n)(np.sqrt(w)*x)

def HO_function(n,x,w):
    constant = 1.0/np.sqrt(2**n * factorial(n))
    constant *= (w/np.pi)**0.25
    return constant*np.exp(-0.5*w*x**2)*hermite(n,x,w)

def HO_energy(n,w,d):
    return w*(n+d/2)
    
def DW_matrix(n,w,d,alpha):
    H_dw = np.zeros((n,n))
    np.fill_diagonal(H_dw,HO_energy(np.arange(n),w,d) + alpha**2/8.0*w*w)
    for l in range(n):
        psi_l_ho = HO_function(l,x,w)
        for m in range(n):
            psi_m_ho = HO_function(m,x,w)
            psi_l_psi_m = simps(psi_l_ho*np.abs(x)*psi_m_ho,x)
            H_dw[l,m] -= 0.5*w*w*alpha*psi_l_psi_m
    return H_dw
    
def find_C(n,w,d,alpha):
    H_dw = DW_matrix(n,w,d,alpha)
    eps, C = np.linalg.eigh(H_dw)
    print("Ground state energy: ", eps[0])
    return C
    
def print_C_to_file(n,w,d,alpha,path='../data/int1/doublewell/'):
    C = find_C(n,w,d,alpha)
    
    np.savetxt(path+str(w)+'00000w/coeffs.dat',C)
    
def DW_function(k,x,n,w,d,alpha):
    C = find_C(n,w,d,alpha) 
    psi_k = np.zeros(len(x))
    for i in range(n):
        psi_k += C[i,k]*HO_function(i,x,w)
    return psi_k
    
def plot_DW_function(k,x,n,w,d,alpha):
    plt.figure(k)
    plt.plot(x,np.power(DW_function(k,x,n_ho,w,d,alpha),2))
    plt.title("%d"%k)
    plt.xlabel("x")
    plt.ylabel("$\psi(x)$")
    plt.grid()
    plt.show()

if __name__ == '__main__':
    Ngrid = 1000
    Lx = 10
    x = np.linspace(-Lx,Lx,Ngrid+1)
    alpha = 2
    w = 1.0
    d = 1
    n_ho = 12 # number of harmonic oscillator basis functions

    print_C_to_file(n_ho,w,d,alpha)
    plot_DW_function(0,x,n_ho,w,d,alpha)
    #for i in range(n_ho):
    #    plot_DW_function(i,x,n_ho,w,d,alpha)
    
    #print(DW_function(2,np.array([0.35]),n_ho,w,d,alpha))
