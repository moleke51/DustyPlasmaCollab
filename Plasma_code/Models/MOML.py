import scipy as sp  
import numpy as np
import scipy.special as sps
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,bisect

Theta = np.logspace(-3,3,1001)
mu = 43*np.ones(len(Theta))
z = 1*np.ones(len(Theta))
alpha = 50*np.ones(len(Theta))
upsilon = 0*np.ones(len(Theta))


def get_name():
    return "MOML"

def colour():
    return 'red'

#def realLambertW(x):
#    if type(x) is float or type(x) is int or type(x) is np.float64 :
#        w = sps.lambertw(x)
#        if np.imag(w)==0:
#            W = np.real(w)
#            return(W)
#        else:
#            return('This value is outside the accepted range')
#    elif type(x) is np.ndarray or type(x) is list:
#        W = np.zeros(len(x))
#        for i in range(0,len(x)):
#            
#            if np.imag(sps.lambertw(x[i]))==0:
#                W[i]= np.real(sps.lambertw(x[i]))
#            else:
#                print('The value at position ' +str(i)+' is outside the accepted range')
#        return(W)
#    else:
#        return('This is an invalid input')

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
#def potential_finder(Theta,mu,z,alpha,upsilon): #gamma = 5/3 for static plasmas
#    Phi = (Theta/z) - realLambertW((np.sqrt(2*np.pi*Theta*(1+(5/3)*Theta)))*np.exp(Theta/z)) + 0.5*np.log(((z**2)*2*np.pi*(1+(5/3)*Theta))/((mu)**2))
#    return np.absolute(Phi) #returned phi is positive  
#    
#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
#Define MOML equation to solve 
def MOML_function(Phi,Theta,mu,z,alpha,upsilon):
    return (np.sqrt(Theta)/mu)*(1 - (1/Theta))*(Phi - 0.5*(np.log(2*np.pi*(1+(5/3)*Theta))-np.log(mu**2))) - np.exp(Phi)

#Solve OML equation for Phi
def potential_finder(Theta,mu,z,alpha,upsilon):
    Phi_list = []
    for i in range(len(Theta)):
        Phi = bisect(MOML_function,-10,10,args = (Theta[i],mu[i],z[i],alpha[i],upsilon[i]))
        Phi_list.append(Phi)
    Phi = np.array(Phi_list)
    return Phi

Phi = potential_finder(Theta,mu,z,alpha,upsilon)
plt.plot(Theta,Phi)
plt.xscale("log")
plt.show()

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

print(potential_finder(0.01,42.82,1,100,0))