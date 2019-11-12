#!/usr/bin/python

import numpy as np

def none(x, der=False):
    '''Pure linear activation'''
    if der:
        return 1
    else:
        return x

def logistic(x, der=False):
    '''Maps the argument x in the interval [0, 1]'''
    if der:
        return  x*(1 - x)
    else:
        return (1 + np.exp(-x))**(-1)

def tanh(x, der=False):
    '''Maps the argument x in the interval [-1, 1]'''
    if der:
        return (np.cosh(x))**(-2) 
    else:
        return np.tanh(x) 
    
def ReLU(x, der=False):
    '''Standard ReLU'''
    if der:
        return np.where(x>0, 1, 0)
    else:
        return np.where(x>0, x, 0)

def ELU(x, der=False, a=1):
    '''ELU'''
    if der:
        return np.where(x<0, a*np.exp(x), x)
    else:
        return np.where(x<0, a*(np.exp(x)-1), x)
        
def Leaky_ReLU(x, der=False, a=0.1):
    '''Leaky ReLU'''
    if der:
        return np.where(x>0, 1, a)
    else:
        return np.where(x<0, a*x, x)
        
def Even_ReLU(x, der=False, a=100):
    '''Even ReLU'''
    if der:
        return np.where(x<0, 1+x/a, 1)
    else:
        return np.where(x<0, x+x*x/(2*a), x)
    


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    plt.style.use("bmh")
    plt.rcParams["font.family"] = "Serif"

    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.gcf().subplots_adjust(left=0.20)

    x = np.linspace(-8, 8, 1000)  
    label_size = {'size':'18'}
    
    ax = plt.gca()
    ax.set_facecolor('white')
    plt.axhline( 0, linestyle = '--', color = 'r', linewidth=0.5)
    plt.axhline( 1, linestyle = '--', color = 'r', linewidth=0.5)
    plt.plot(x, logistic(x), linewidth = 2)

    plt.title(r"Sigmoid",**label_size)
    plt.xlabel("$x$",**label_size)
    plt.ylabel("$f(x)$",**label_size)
    plt.text(-8, 0.9, r'$f(x)=(1 + \exp(-x))^{-1}$',**label_size)
    plt.savefig("../plots/sigmoid.png")
    plt.show()

    ax = plt.gca()
    ax.set_facecolor('white')
    plt.axhline(-1, linestyle = '--', color = 'r', linewidth=0.5)
    plt.axhline( 1, linestyle = '--', color = 'r', linewidth=0.5)
    plt.plot(x, tanh(x), linewidth = 2)

    plt.title(r"Tangens hyperbolicus",**label_size)
    plt.xlabel("$x$",**label_size)
    plt.ylabel("$f(x)$",**label_size)
    plt.text(-8, 0.8, r'$f(x)=tanh(x)$',**label_size)
    plt.savefig("../plots/tanh.png")
    plt.show()
    
    ax = plt.gca()
    ax.set_facecolor('white')
    ReLU_array = np.empty(len(x))
    for i in range(len(x)):
        ReLU_array[i] = ReLU(x[i])
        
    plt.plot(x, ReLU_array, linewidth = 2)
    plt.title(r"ReLU",**label_size)
    plt.xlabel("$x$",**label_size)
    plt.ylabel("$f(x)$",**label_size)
    plt.text(-8, 7, r'$f(x)=max(0,x)$',**label_size)
    plt.savefig("../plots/ReLU.png")
    plt.show()
    
    ax = plt.gca()
    ax.set_facecolor('white')
    LeakyReLU_array = np.empty(len(x))
    for i in range(len(x)):
        LeakyReLU_array[i] = Leaky_ReLU(x[i])
    
    plt.plot(x, LeakyReLU_array, linewidth = 2)
    plt.title(r"Leaky ReLU",**label_size)
    plt.xlabel("$x$",**label_size)
    plt.ylabel("$f(x)$",**label_size)
    plt.text(-8, 7, r'$f(x)=0.1x$ if $x<0$, $x$ else',**label_size)
    plt.savefig("../plots/LeakyReLU.png")
    plt.show()
    
    ax = plt.gca()
    ax.set_facecolor('white')
    ELU_array = np.empty(len(x))
    for i in range(len(x)):
        ELU_array[i] = ELU(x[i])
    
    plt.plot(x, ELU_array, linewidth = 2)
    plt.title(r"ELU",**label_size)
    plt.xlabel("$x$",**label_size)
    plt.ylabel("$f(x)$",**label_size)
    plt.text(-8, 7, r'$f(x)=exp(x)-1$ if $x<0$, $x$ else',**label_size)
    plt.savefig("../plots/ELU.png")
    plt.show()

