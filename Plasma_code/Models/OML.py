import scipy as sp  
import numpy as np
import scipy.special as sps
from scipy.optimize import fsolve,bisect


def get_name():
    return "OML"

def colour():
    return 'orange'

def get_info():
    assumptions_list = ["Model assumptions:\n","Spherical symmetry\n","No collisions\n","No magnetic field\n","No external electric field\n","No electron emission of any kind\n","Quasi-neutrality in bulk plasma\n","Conservation of particle energy\n","Conservation of particle angular momentum\n","Limiting trajectory is the grazing incidence\n"]
    validity_list = ["Validity:\n","Static plasma\n","Small dust\n","Any ion temperature\n"]
    reference_list = ["References:\n","C. T. N. Willis, “Dust in stationary and flowing plasmas,” Physics PhD Thesis, Imperial College London, March 2012\n","D. M. Thomas, “Theory and simulation of the charging of dust in plasmas,” Physics PhD Thesis, Imperial College London, March 2016\n","K. R. V. and A. J. E., “The floating potential of spherical probes and dust grains. ii: Orbital motion theory,” Journal of Plasma Physics, vol. 69.6, pp. 485–506, 2002"]
    string = " ".join(assumptions_list) + " ".join(validity_list) + " ".join(reference_list)
    return string

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
#Define OML equation to solve 
def OML_function(Phi,Theta,mu,z,alpha,upsilon):
    return (np.sqrt(Theta)/mu)*(1 - (z/Theta)*Phi) - np.exp(Phi)

#Solve OML equation for Phi
def potential_finder(dictionarylist):
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
            
    Phi = bisect(OML_function,-10,10,args = (Theta,mu,z,alpha,upsilon))
    return np.absolute(Phi)

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
    if alpha > 1.25*Theta**(0.4):
        P_a = 0
    else:
        P_a = 1
    if upsilon > 0:
        P_u = 0
    else:
        P_u = 1
    return (P_t + P_a + P_u)

  