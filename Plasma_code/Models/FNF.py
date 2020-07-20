import scipy as sp  
import numpy as np
import scipy.special as sps
import SOML as soml
import SMOML as smoml


def get_name():
    return "FNF"

def colour():
    return 'purple'


#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_SOML,phi_TS,alpha_SOML,alpha_TS,alpha):
    x = ((phi_TS - phi_SOML)/(np.log(alpha_TS) - np.log(alpha_SOML)))*np.log((alpha)/(alpha_TS)) + phi_TS
    return x 

def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 3 for flowing plasmas
    alpha_OML = 1.25*(Theta)**0.4 #Assume this is the same as the static case
    alpha_TS = 50
    Phi_SMOML = smoml.potential_finder(Theta,mu,z,alpha,upsilon)
    Phi_SOML = soml.potential_finder(Theta,mu,z,alpha,upsilon)
    Phi = Linear_function(Phi_SOML,Phi_SMOML,alpha_OML,alpha_TS,alpha)
    return Phi #returned phi is positive

def priority(Theta,alpha,upsilon):
    if Theta >= 1e-4:
        P_t = 1
    else:
        P_t = 0.5
    if alpha < 50 and alpha > 1.25*Theta**(0.4):
        P_a = 1
    else:
        P_a = 0
    if upsilon > 0:
        P_u = 1
    else:
        P_u = 0
    return (P_t + P_a + P_u)   