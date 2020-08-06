import scipy as sp  
import numpy as np
import scipy.special as sps
import OML as oml
import MOML as moml
from termcolor import colored

def get_name():
    return "SNF"

def colour():
    return 'brown'

def get_info():
    assumptions_list = ["Model assumptions:\n","Spherical symmetry\n","No collisions\n","No magnetic field\n","No external electric field\n","No electron emission of any kind\n","Quasi-neutrality in bulk plasma\n"]
    validity_list = ["Validity:\n","Static plasma\n","Intermediate sized dust\n","Any ion temperature\n"]
    reference_list = ["References:\n","C. T. N. Willis, “Dust in stationary and flowing plasmas,” Physics PhD Thesis, Imperial College London, March 2012\n","D. M. Thomas, “Theory and simulation of the charging of dust in plasmas,” Physics PhD Thesis, Imperial College London, March 2016\n"]
    string = " ".join(assumptions_list) + " ".join(validity_list) + " ".join(reference_list)
    return print(colored(string,'yellow'))

#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_OML,phi_TS,alpha_OML,alpha_TS,alpha):
    x = ((phi_TS - phi_OML)/(np.log(alpha_TS) - np.log(alpha_OML)))*np.log((alpha)/(alpha_TS)) + phi_TS
    return x 

def potential_finder(dictionarylist): #gamma = 5/3 for static plasmas
    for _vardict in dictionarylist:
        if _vardict.get('Norm_var_name') != None:
            if _vardict.get('Norm_var_name') == 'alpha':
                alpha = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'z':
                z = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'mu':
                mu = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'upsilon':
                upsilon = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'Theta':
                Theta = _vardict.get('Norm_value')
            
    alpha_OML = 1.25*(Theta)**0.4
    alpha_TS = 50
    Phi_MOML = moml.potential_finder(Theta,mu,z,alpha,upsilon)
    Phi_OML = oml.potential_finder(Theta,mu,z,alpha,upsilon)
    Phi = Linear_function(Phi_OML,Phi_MOML,alpha_OML,alpha_TS,alpha)
    return np.absolute(Phi) #returned phi is positive

def priority(dictionarylist):
    for _vardict in dictionarylist:
        if _vardict.get('Norm_var_name') != None:
            if _vardict.get('Norm_var_name') == 'alpha':
                alpha = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'z':
                z = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'mu':
                mu = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'upsilon':
                upsilon = _vardict.get('Norm_value')
            elif _vardict.get('Norm_var_name') == 'Theta':
                Theta = _vardict.get('Norm_value')
            else:
                if _vardict.get('Norm_value') != _vardict.get('default value'):
                    return 0
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
