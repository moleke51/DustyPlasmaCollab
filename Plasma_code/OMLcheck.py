import numpy as np
import scipy as sp
from scipy.optimize import bisect
import matplotlib.pyplot as plt
Theta = np.logspace(-7,5,100)
mu = 43
z = 1

def COML(Phi,mu,Theta,z):
    return Phi + (mu*np.sqrt(Theta)/z)*(1+Phi)*np.exp(Phi) - (Theta/z)
def Phi_solve(mu,Theta,z):
    return bisect(COML,-10,10,args = (mu,Theta,z))
def Phi_sol(mu,Theta,z):
    return (Theta - mu*np.sqrt(Theta))/(z+mu*np.sqrt(Theta))
print(Phi_solve(43,0,1))
Phi_a = np.ones(len(Theta))
for i in range(len(Theta)):
    Phi_a[i] *= Phi_sol(mu,Theta[i],z)

plt.plot(Theta,np.zeros(len(Theta)))
plt.plot(Theta,Phi_a)
plt.xscale('log')
plt.show()