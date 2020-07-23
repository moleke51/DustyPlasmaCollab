import scipy as sp  
import numpy as np
import scipy.special as sps
import OML as oml
import MOML as moml
import matplotlib.pyplot as plt

def get_name():
    return "SNF"

def colour():
    return 'brown'


#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_OML,phi_TS,alpha_OML,alpha_TS,alpha):
    x = ((phi_TS - phi_OML)/(np.log(alpha_TS) - np.log(alpha_OML)))*np.log((alpha)/(alpha_TS)) + phi_TS
    return x 

def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 5/3 for static plasmas
    alpha_OML = 1.25*(Theta)**0.4
    alpha_TS = 50
    Phi_MOML = moml.potential_finder(Theta,mu,z,alpha,upsilon)
    Phi_OML = oml.potential_finder(Theta,mu,z,alpha,upsilon)
    Phi = Linear_function(Phi_OML,Phi_MOML,alpha_OML,alpha_TS,alpha)
    return np.absolute(Phi) #returned phi is positive

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
        P_u = 0
    else:
        P_u = 1
    return (P_t + P_a + P_u)    
