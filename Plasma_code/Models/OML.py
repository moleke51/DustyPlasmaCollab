import scipy as sp  
import numpy as np
import scipy.special as sps
from scipy.optimize import fsolve,bisect
from termcolor import colored

def get_name():
    return "OML"

def colour():
    return 'orange'

def get_info():
    string = ("Base assumptions: Spherical symmetry; no collisions; no magnetic field; no external electric field; no electron emission of any kind; quasi-neutrality in bulk plasma.\n"
              "Model assumptions: Conservation of particle energy; conservation of particle angular momentum; limiting trajectory is the grazing incidence.\n"
              "Vality: Static plasma; any " + "\u0398" + "; small " + "\u03B1" + " (" + "\u03B1" " less than or equal to 1.25*" + "\u0398" + "^(0.4)).\n"
              "References: C. T. N. Willis, “Dust in stationary and flowing plasmas,” Physics PhD Thesis, Imperial College London, March 2012;\n" +
              "D. M. Thomas, “Theory and simulation of the charging of dust in plasmas,” Physics PhD Thesis, Imperial College London, March 2016;\n" +
              "K. R. V. and A. J. E., “The floating potential of spherical probes and dust grains. ii: Orbital motion theory,” Journal of Plasma Physics, vol. 69.6, pp. 485–506, 2002.")
    return print(colored(string,'blue'))

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
  
