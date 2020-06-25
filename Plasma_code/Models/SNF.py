import scipy as sp  
import numpy as np
import scipy.special as sps


def get_name():
    return "SNF"

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

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis (absolute values)
def OML_surface_potential_finder_abs(Theta,mu,z): 
    x = np.absolute((Theta/z) - realLambertW((mu*np.sqrt(Theta)/z)*np.exp(Theta/z)))
    return x

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis (absolute values)
def MOML_surface_potential_finder_abs(Theta,mu,z): #gamma = 5/3
    x = np.absolute(Theta - realLambertW((np.sqrt(2*np.pi*Theta*(1+(5/3)*Theta)))*np.exp(Theta)) + 0.5*np.log((2*np.pi*(1+(5/3)*Theta))/((z**2)*(mu)**2)))
    return x

#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_OML,phi_TS,alpha_OML,alpha_TS,alpha):
    x = ((phi_TS - phi_OML)/(np.log(alpha_TS) - np.log(alpha_OML)))*np.log((alpha)/(alpha_TS)) + phi_TS
    return x 

def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 5/3
    alpha_OML = 1.25*(Theta)**0.4
    alpha_TS = 50
    Phi_MOML = MOML_surface_potential_finder_abs(Theta,mu,z) 
    Phi_OML = OML_surface_potential_finder_abs(Theta,mu,z)
    Phi = Linear_function(Phi_OML,Phi_MOML,alpha_OML,alpha_TS,alpha)
    return Phi #returned phi is positive

def priority():
    return 1    