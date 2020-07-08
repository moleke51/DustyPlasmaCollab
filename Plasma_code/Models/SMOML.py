import scipy as sp  
import numpy as np
import scipy.special as sps


def get_name():
    return "SMOML"

def colour():
    return 'black'

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
def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 3 for flowing plasmas
    s_1 = ((np.sqrt(np.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*np.exp(-(upsilon**2))
    s_2 = (np.sqrt(np.pi)*sps.erf(upsilon))/(2*upsilon)
    Phi = (Theta*s_1)/(z*s_2) - realLambertW(((np.sqrt(2*np.pi*Theta*(1+(3)*Theta)))/(s_2))*np.exp((Theta*s_1)/(z*s_2))) + 0.5*np.log(((z**2)*2*np.pi*(1+(3)*Theta))/((mu)**2))
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