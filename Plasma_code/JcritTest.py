#============================IMPORT STUFF==========================#
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sps

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

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
def potential_finder(Theta,mu,z,alpha,upsilon): 
    Phi = (Theta/z) - realLambertW((mu*np.sqrt(Theta)/z)*np.exp(Theta/z))
    return np.absolute(Phi) #returned phi is positive

def L_crit(E,Phi):
    return np.sqrt(2*E)*np.sqrt(1 + (Phi/E))

Theta = 1
mu = 43
z = 1 
alpha = 0
upsilon = 0

E = np.linspace(0,5,100)

Phi_oml = potential_finder(Theta,mu,z,alpha,upsilon)*np.ones(len(E))
Phi_ps = 0.5*np.ones(len(E))
Phi_sim = 0.8*np.ones(len(E))#From Willis' paper

L_oml = L_crit(E,Phi_oml)
L_ps = L_crit(E,Phi_ps)
L_sim = L_crit(E,Phi_sim)

plt.plot(E,L_oml,color = "Red", label = "L_oml")
plt.plot(E,L_ps,color = "Blue", label = "L_ps")
plt.plot(E,L_sim,color = "Green", label = "L_sim")
plt.title("Variation of critical angular momentum with ion energy")
plt.ylabel("Critical angular momentum, L_crit")
plt.xlabel("Ion energy")
plt.legend()
plt.grid()
plt.show()