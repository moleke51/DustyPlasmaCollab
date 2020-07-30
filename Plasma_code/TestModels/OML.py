import scipy as sp  
import numpy as np
import scipy.special as sps
from scipy.optimize import fsolve,bisect

def get_name():
    return "OML"

def colour():
    return 'orange'

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
#Define OML equation to solve 
def OML_function(Phi,Theta,mu,z,alpha,upsilon):
    return (np.sqrt(Theta)/mu)*(1 - (z/Theta)*Phi) - np.exp(Phi)

#Solve OML equation for Phi
def potential_finder(Theta,mu,z,alpha,upsilon):
    Phi = bisect(OML_function,-10,10,args = (Theta,mu,z,alpha,upsilon))
    return np.absolute(Phi)

def priority(Theta,alpha,upsilon):
    if Theta >= 1e-4:
        P_t = 1
    else:
        P_t = 0.5
    if alpha > 1.25*Theta**(0.4):
        P_a = 0
    else:
        P_a = 1
    if upsilon > 0:
        P_u = 0
    else:
        P_u = 1
    return (P_t + P_a + P_u)
  
print(potential_finder(1,43,1,0.01,0))