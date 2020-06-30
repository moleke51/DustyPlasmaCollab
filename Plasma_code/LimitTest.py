import numpy as np
import scipy as sp
import scipy.special as sps
from scipy.optimize import bisect

def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64 :
        w = sps.lambertw(x)
        if np.imag(w)==0:
            W = np.real(w)
            return(W)
        else:
            return('This value is outside the accepted range')
    elif type(x) is np.ndarray or type(x) is list:
        W = np.zeros(len(x))
        for i in range(0,len(x)):
            
            if np.imag(sps.lambertw(x[i]))==0:
                W[i]= np.real(sps.lambertw(x[i]))
            else:
                print('The value at position ' +str(i)+' is outside the accepted range')
        return(W)
    else:
        return('This is an invalid input')

mu = 43
z = 1

def MOML_limit(Phi,mu,z):
    return ((1+Phi)*np.exp(Phi) - np.sqrt(2*np.pi)/mu*z)
def MOML_limit_solve(mu,z):
    return bisect(MOML_limit,-10,10,args = (mu,z))
#Phi_a = fsolve(MOML_limit, -3, args = (mu,z))
print(MOML_limit_solve(mu,z))
