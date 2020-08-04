import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve,bisect
import scipy.special as sps
pi = np.pi
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
'''
Theta = np.logspace(-5,1,100001)
mu = 43
z = 1
#Phi_1 = np.zeros(len(Theta)) #gamma = 1
Phi = np.zeros(len(Theta)) #gamma = 5/3
#Phi_3 = np.zeros(len(Theta)) #gamma = 3
for i in range(len(Theta)):
    #Phi_1[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,1))
    Phi[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,5/3))
    #Phi_3[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,3))

#Find the points where MOML and FlowingSheath are the closest, within 0.01%
index = []
Phi_intersect = []
Phi_flow = FlowingSheathMOML(Theta,43*np.ones(len(Theta)),(5/3)*np.ones(len(Theta)),2*np.ones(len(Theta)))
for i in range(len(Theta)):
    if (np.absolute(Phi_flow[i] - Phi[i]) < 0.0001):
        index.append(i)
        Phi_intersect.append(Phi_flow[i])

#Find the biggest and smallest points 
Phi_intersect_index1 = np.argmin(Phi_intersect)
Phi_intersect_index2 = np.argmax(Phi_intersect)    

#plt.plot(Theta,Phi_1,color = "Red", label = "Gamma = 1")
plt.plot(Theta,Phi,color = "Black", label = "Gamma = 5/3")
#plt.plot(Theta,Phi_3,color = "Green", label = "Gamma = 3")
plt.plot(Theta,Phi_flow,color = 'red',label = 'Flowing sheath, gamma = 5/3, kappa = 2')
plt.plot(Theta[index[Phi_intersect_index1]],Phi_flow[index[Phi_intersect_index1]],"o",color = "Green")
plt.plot(Theta[index[Phi_intersect_index2]],Phi_flow[index[Phi_intersect_index2]],"o",color = "Green")
plt.xlabel("Theta")
plt.ylabel("Normalised dust potential")
plt.legend()
plt.grid()
plt.xscale("log")
plt.show()


a = 0
while a == 0:
    print(FlowingSheathMOML(float(input('Select theta value: ')),43,5/3,2))
'''
def ThetaIntersect(Theta,mu,gamma,kappa):
    return realLambertW(np.sqrt(2*pi*Theta*(1+gamma*Theta))*np.exp(Theta)) - 0.5*(1+Theta*(gamma+5-2*kappa))
print(fsolve(ThetaIntersect,10,args = (43,5/3,2)))