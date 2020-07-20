import scipy as sp  
import numpy as np
import scipy.special as sps
from scipy.optimize import fsolve,bisect

def get_name():
    return "SMOML"

def colour():
    return 'black'

#SMOML (Shifted Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis
#Define SMOML equation to solve
def SMOML_function(Phi,Theta,mu,z,alpha,upsilon): #gamma = 3 for flowing plasmas
    s_1 = ((np.sqrt(np.pi))*(1+2*(upsilon**2))* sps.erf(upsilon))/(4*upsilon) + 0.5*np.exp(-(upsilon**2))
    s_2 = (np.sqrt(np.pi)* sps.erf(upsilon))/(2*upsilon)
    return (np.sqrt(Theta)/mu)*(s_1 - (s_2/Theta)*(Phi - 0.5*(np.log(2*np.pi*(1+(5/3)*Theta))-np.log(mu**2)))) - np.exp(Phi)


def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 3 for flowing plasmas
    Phi = bisect(SMOML_function,-10,10,args = (Theta,mu,z,alpha,upsilon))
    return np.absolute(Phi) #returned phi is positive

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
        P_u = 1
    else:
        P_u = 0
    return (P_t + P_a + P_u)    

