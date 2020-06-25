import scipy as sp  
import numpy as np
import scipy.special as sps


def get_name():
    return "OML"

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

def priority():
    return 1
  
