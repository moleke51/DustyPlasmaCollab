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

#SOML (Shifted OML) model for normalised dust surface potential - eqn 2.138 in Thomas' thesis (absolute values)
def SOML_surface_potential_finder_abs(Theta,mu,z,upsilon): 
    s_1 = ((np.sqrt(np.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*np.exp(-(upsilon**2))
    s_2 = (np.sqrt(np.pi)*sps.erf(upsilon))/(2*upsilon)
    x = np.absolute((Theta*s_1)/(s_2) - realLambertW(((mu*z*np.sqrt(Theta))/(s_2))*np.exp((Theta*s_1)/(s_2))))
    return x

#SMOML (Shifted Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis (absolute values)
def SMOML_surface_potential_finder_abs(Theta,mu,z,upsilon): #gamma = 5/3
    s_1 = ((np.sqrt(np.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*np.exp(-(upsilon**2))
    s_2 = (np.sqrt(np.pi)*sps.erf(upsilon))/(2*upsilon)
    x = np.absolute((Theta*s_1)/(s_2) - realLambertW(((np.sqrt(2*np.pi*Theta*(1+(5/3)*Theta)))/(s_2))*np.exp((Theta*s_1)/(s_2))) + 0.5*np.log((2*np.pi*(1+(5/3)*Theta))/((z**2)*(mu)**2)))
    return x

#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_SOML,phi_TS,alpha_SOML,alpha_TS,alpha):
    x = ((phi_TS - phi_SOML)/(np.log(alpha_TS) - np.log(alpha_SOML)))*np.log((alpha)/(alpha_TS)) + phi_TS
    return x 

def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 5/3
    alpha_OML = 1.25*(Theta)**0.4 #Assume this is the same as the static case
    alpha_TS = 50
    Phi_SMOML = SMOML_surface_potential_finder_abs(Theta,mu,z,upsilon) 
    Phi_SOML = SOML_surface_potential_finder_abs(Theta,mu,z,upsilon)
    Phi = Linear_function(Phi_SOML,Phi_SMOML,alpha_OML,alpha_TS,alpha)
    return Phi #returned phi is positive

def priority():
    return 1    