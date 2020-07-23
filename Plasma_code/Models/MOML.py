import scipy as sp  
import numpy as np
import scipy.special as sps
from scipy.optimize import fsolve,bisect

def get_name():
    return "MOML"

def colour():
    return 'red'

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
#Define MOML equation to solve 
def MOML_function(Phi,Theta,mu,z,alpha,upsilon): #gamma = 5/3 for static plasmas
    return (np.sqrt(Theta)/mu)*(1 - (1/Theta)*(Phi - 0.5*(np.log(2*np.pi*(1+(5/3)*Theta))-np.log(mu**2)))) - np.exp(Phi)

#Solve MOML equation for Phi
def potential_finder(Theta,mu,z,alpha,upsilon):
    Phi = bisect(MOML_function,-10,10,args = (Theta,mu,z,alpha,upsilon))
    return np.absolute(Phi)
def planar_presheath(Theta,mu,z,alpha,upsilon):
    return -0.5*np.log((2*np.pi*(1+(5/3)*Theta)/(mu**2))) + np.log(2)

def priority(Theta,alpha,upsilon):
    if Theta >= 1e-4:
        P_t = 1
    else:
        P_t = 0.5
    if alpha >= 50:
        P_a = 1
    else:
        P_a = 0
    if upsilon > 0:
        P_u = 0
    else:
        P_u = 1
    return (P_t + P_a + P_u)             

print(potential_finder(1,43,1,100,0))
print(planar_presheath(1,43,1,100,0))