import numpy as np
import matplotlib.pyplot as plt

plt.style.use("bmh")
ax = plt.gca()
ax.set_facecolor('white')

x_ = np.linspace(0,16,1000)
y_ = [0, 15, 28, 39, 48, 55, 60, 63, 64, 63, 60, 55, 48, 39, 28, 15, 0]

x = np.array([1, 9, 10, 2, 4, 6, 7, 11, 13, 16])
y = np.array([15, 63, 60, 30, 50, 60, 65, 55, 40, 0])

x_val = np.array([1, 9, 10]) 
y_val = np.array([15, 63, 60]) 

size_label = 20
size_legend = 16
label_size = {'size':str(size_label)}
plt.rcParams["font.family"] = "Serif"

p=6
X = np.zeros((len(x),p+1))
for j in range(p+1):
    for i in range(len(x)):
        X[i,j] = x[i]**j
        
X_new = X.T.dot(X)        
INV = np.linalg.inv(X_new)
beta = INV.dot(X.T.dot(y))
print(beta)

def f(x):
    summer = 0
    for i in range(p+1):
        summer += beta[i]*x**i
    return summer
print(np.sum((f(x)-y)**2)/len(x))
print(np.sum((f(x_val)-y_val)**2)/len(x_val))

plt.plot(x, y, 'or', label='Data points')  
'''
p=2
X = np.zeros((len(x),p+1))
for j in range(p+1):
    for i in range(len(x)):
        X[i,j] = x[i]**j
        
X_new = X.T.dot(X)        
INV = np.linalg.inv(X_new)
beta = INV.dot(X.T.dot(y))
print(beta)

def f(x):
    summer = 0
    for i in range(p+1):
        summer += beta[i]*x**i
    return summer
print(np.sum((f(x)-y)**2)/len(x))
print(np.sum((f(x_val)-y_val)**2)/len(x_val))

plt.plot(x_, f(x_), '-.', label='2nd order')   

p=1
X = np.zeros((len(x),p+1))
for j in range(p+1):
    for i in range(len(x)):
        X[i,j] = x[i]**j
        
X_new = X.T.dot(X)        
INV = np.linalg.inv(X_new)
beta = INV.dot(X.T.dot(y))
print(beta)

def f(x):
    summer = 0
    for i in range(p+1):
        summer += beta[i]*x**i
    return summer
print(np.sum((f(x)-y)**2)/len(x))
print(np.sum((f(x_val)-y_val)**2)/len(x_val))

plt.plot(x_, f(x_), ':', label='1st order')  

plt.plot(x,y, 'or',label='Training points')
plt.plot(x_val,y_val,'vk',label='Validation points')
'''
plt.xlabel('$x$', **label_size)
plt.ylabel('$y$', **label_size) 
plt.legend(loc='lower center', fontsize=size_legend, facecolor='white', framealpha=1)
plt.axis([-0.1, 16.1, -1, 71])
plt.show()

#import tikzplotlib

#tikzplotlib.save("/home/evenmn/Master-thesis/doc/text/pgf/datapointsplot.tex")
