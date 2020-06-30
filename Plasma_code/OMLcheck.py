import numpy as np
import scipy as sp
from scipy.optimize import bisect

Theta = 1
mu = 43
z = 1

def COML(Phi,mu,Theta,z):
    return Phi + (mu*np.sqrt(Theta)/z)*(1+Phi)*np.exp(Phi) - (Theta/z)


Phi_a =bisect(COML,-10,0,args = (mu,Theta,z))
print(Phi_a)
print(COML(Phi_a,mu,Theta,z))