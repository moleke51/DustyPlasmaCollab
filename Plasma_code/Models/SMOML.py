import scipy as sp  
import numpy as np
import scipy.special as sps
from scipy.optimize import fsolve,bisect
from termcolor import colored

def get_name():
    return "SMOML"

def colour():
    return 'black'

def get_info():
    string = ("Base assumptions: Spherical symmetry; no collisions; no magnetic field; no external electric field; no electron emission of any kind; quasi-neutrality in bulk plasma.\n"
              "Model assumptions: Conservation of particle energy; conservation of particle angular momentum; current at sheath edge is the same as current at dust surface; flowing" + "\n" + "Maxwell-Boltzmann velocity distribution for ions.\n"
              "Vality: Flowing plasma; any " + "\u0398" + "; large " + "\u03B1" + " (" + "\u03B1" " greater than or equal to 100).\n"
              "References: C. T. N. Willis, “Dust in stationary and flowing plasmas,” Physics PhD Thesis, Imperial College London, March 2012;\n" +
              "D. M. Thomas, “Theory and simulation of the charging of dust in plasmas,” Physics PhD Thesis, Imperial College London, March 2016.")
    return print(colored(string,'blue'))

#SMOML (Shifted Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis
#Define SMOML equation to solve
def SMOML_function(Phi,Theta,mu,z,alpha,upsilon): #gamma = 3 for flowing plasmas
    s_1 = ((np.sqrt(np.pi))*(1+2*(upsilon**2))* sps.erf(upsilon))/(4*upsilon) + 0.5*np.exp(-(upsilon**2))
    s_2 = (np.sqrt(np.pi)* sps.erf(upsilon))/(2*upsilon)
    return (np.sqrt(Theta)/mu)*(s_1 - (s_2/Theta)*(Phi - 0.5*(np.log(2*np.pi*(1+(5/3)*Theta))-np.log(mu**2)))) - np.exp(Phi)


def potential_finder(dictionarylist): #gamma = 3 for flowing plasmas
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
            
    Phi = bisect(SMOML_function,-10,10,args = (Theta,mu,z,alpha,upsilon))
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
    if alpha >= 50:
        P_a = 1
    else:
        P_a = 0
    if upsilon > 0:
        P_u = 1
    else:
        P_u = 0
    return (P_t + P_a + P_u)    
