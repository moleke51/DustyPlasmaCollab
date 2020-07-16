import scipy as sp  
import numpy as np
import scipy.special as sps


def get_name():
    return "MOML"

def colour():
    return 'red'

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

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 5/3 for static plasmas
    Phi = -0.5*((np.sqrt(1+(5/3)*Theta) - np.sqrt(8*Theta/np.pi))**2) + 0.5*np.log(((z**2)*2*np.pi*(1+(5/3)*Theta))/((mu)**2))
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
        P_u = 0
    else:
        P_u = 1
    return (P_t + P_a + P_u)             

print(potential_finder(1,42.82,1,100,0))