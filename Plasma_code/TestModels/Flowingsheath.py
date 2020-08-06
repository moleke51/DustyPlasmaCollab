import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve,bisect
import scipy.special as sps
from scipy.integrate import simps
import scipy.optimize as spo
pi = np.pi
u = 1.66e-27

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

def FlowingSheathMOML(Theta,kappa):
    mu = 43
    z = 1 
    gamma = 5/3
    return (-0.5/z)*(1+Theta*(gamma+3-2*kappa)) + 0.5*(np.log(2*np.pi*(1+gamma*Theta))-np.log(mu**2))

def ThinSheathFit(Theta):
    Phi = np.zeros(len(Theta))
    for i in range(len(Theta)):
        if Theta[i] <= 2:
            A = 1 #Mass number
            eta = 0.456
            psi = 0
            C = 3.179
            Phi[i] = eta*np.log(A) + psi*np.log(Theta[i]) + C
        else:
            A = 1 #Mass number
            eta = 0.557
            psi = -0.386-0.024*np.log(A)
            C = 3.399
            Phi[i] = eta*np.log(A) + psi*np.log(Theta[i]) + C

    return -1*Phi

#Parameters
Theta = np.logspace(-5,0,1001)
mu = 43
z = 1

#Calculate MOML values
Phi_moml = np.zeros(len(Theta)) #gamma = 5/3
for i in range(len(Theta)):
    Phi_moml[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,5/3))

#Find best fit kappa
Phi_fit = ThinSheathFit(Theta)
kappa,cov_kappa = spo.curve_fit(FlowingSheathMOML,Theta,Phi_fit)
kappa = kappa[0]
print(f"Best fit: kappa = {kappa} +/- {np.sqrt(cov_kappa[0,0])}")

#Find the flowing sheath values
Phi_flow = FlowingSheathMOML(Theta,(2)*np.ones(len(Theta)))

#Redsiduals
Phi_mfit = np.zeros(len(Theta))
Phi_ffit = np.zeros(len(Theta))
for i in range(len(Theta)):
    Phi_mfit[i] = np.absolute((Phi_fit[i] - Phi_moml[i])/Phi_fit[i])*100
    Phi_ffit[i] = np.absolute((Phi_flow[i] - Phi_fit[i])/Phi_fit[i])*100

Phi_marea = simps(Phi_mfit,dx = 0.001)
Phi_farea = simps(Phi_ffit,dx = 0.001)
print(f"Average percentage difference between SCEPTIC fit and flowing sheath is {str(Phi_farea)}%")
print(f"Average percentage difference between SCEPTIC fit and MOML is {str(Phi_marea)}%")

Phi_fit = ThinSheathFit(Theta)
plt.figure(1)
plt.plot(Theta,Phi_moml,color = "Black", label = "MOML, gamma = 5/3")
plt.plot(Theta,Phi_fit,"--",color = 'Blue',label = 'SCEPTIC Fit')
plt.plot(Theta,Phi_flow,color = 'red',label = 'Flowing sheath, gamma = 5/3, kappa = 2')
#plt.plot(Theta[index[Phi_intersect_index1]],Phi_flow[index[Phi_intersect_index1]],"o",color = "Green")
#plt.plot(Theta[index[Phi_intersect_index2]],Phi_flow[index[Phi_intersect_index2]],"o",color = "Green")
plt.title("Normalised potential as a function of theta")
plt.xlabel("Theta")
plt.ylabel("Normalised dust potential")
plt.legend(loc = 'upper right')
plt.grid()
plt.xscale("log")

plt.figure(2)
plt.plot(Theta,Phi_mfit,'--',color = 'Red',label = "MOML")
plt.plot(Theta,Phi_ffit,color = 'Black',label = "Flowing sheath")
plt.title("Percentage difference between models and SCEPTIC fit as a function of theta")
plt.xlabel("Theta")
plt.ylabel("Percentage difference")
plt.legend()
plt.grid()
plt.xscale("log")


plt.show()


#a = 0
#while a == 0:
#    print(FlowingSheathMOML(float(input('Select theta value: ')),43,5/3,2))
