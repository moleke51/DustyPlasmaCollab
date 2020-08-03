import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve,bisect
import scipy.special as sps

def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64:
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
#Define MOML equation to solve 
def MOML_function(Phi,Theta,mu,z,gamma): #gamma = 5/3 for static plasmas
    return (np.sqrt(Theta)/mu)*(1 - (1/Theta)*(Phi - 0.5*(np.log(2*np.pi*(1+gamma*Theta))-np.log(mu**2)))) - np.exp(Phi)

def FlowingSheathMOML(Theta,mu,gamma,kappa):
    return -0.5*(1+Theta*(gamma+3-2*kappa)) + 0.5*(np.log(2*np.pi*(1+gamma*Theta))-np.log(mu**2))
def FlowingPresheath(Theta,mu,gamma,kappa):
    return -0.5*(1+Theta*(gamma+3-2*kappa))
'''
Theta = np.logspace(-5,0,101)
#Theta = np.logspace(-1,0,101)
mu = 43
z = 1
Phi_1 = np.zeros(len(Theta)) #gamma = 1
Phi_2 = np.zeros(len(Theta)) #gamma = 5/3
Phi_3 = np.zeros(len(Theta)) #gamma = 3


for i in range(len(Theta)):
    Phi_1[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,1))
    Phi_2[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,5/3))
    Phi_3[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,3))


Phi_flow = FlowingSheathMOML(Theta,43*np.ones(len(Theta)),(5/3)*np.ones(len(Theta)))

#plt.plot(Theta,Phi_1,color = "Red", label = "Gamma = 1")
plt.plot(Theta,Phi_2,color = "Black", label = "Gamma = 5/3")
#plt.plot(Theta,Phi_3,color = "Green", label = "Gamma = 3")
plt.plot(Theta,Phi_flow,color = 'red',label = 'Flowing sheath')

plt.xlabel("Theta")
plt.ylabel("Normalised dust potential")
plt.legend()
plt.grid()
plt.xscale("log")
plt.show()

a = 0
while a == 0:
    Theta = float(input('Select theta value: '))
    #print(FlowingSheathMOML(float(input('Select theta value: ')),43,5/3,2))
    print(bisect(MOML_function,-10,10,args = (Theta,43,1,5/3))+FlowingPresheath(Theta,43,5/3,2))

'''
