import scipy as sp  
import numpy as np
import scipy.special as sps


def get_name():
    return "SMOML"

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

#SMOML (Shifted Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis
def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 5/3
    s_1 = ((np.sqrt(np.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*np.exp(-(upsilon**2))
    s_2 = (np.sqrt(np.pi)*sps.erf(upsilon))/(2*upsilon)
    Phi = (Theta*s_1)/(s_2) - realLambertW(((np.sqrt(2*np.pi*Theta*(1+(5/3)*Theta)))/(s_2))*np.exp((Theta*s_1)/(s_2))) + 0.5*np.log((2*np.pi*(1+(5/3)*Theta))/((z**2)*(mu)**2))
    return Phi #returned phi is positive

def priority():
    return 1    