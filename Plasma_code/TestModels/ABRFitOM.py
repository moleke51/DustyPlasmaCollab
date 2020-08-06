import scipy as sp 
import numpy as np
import matplotlib.pyplot as plt
import HotABR as abr
import scipy.optimize as spo
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

'''
def function(alpha,kappa):
    Theta = 1
    mu = 43
    z = 1
    upsilon = 0
    Phi_list = []
    for i in range(len(alpha)):
        Phi = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,kappa)
        Phi_list.append(Phi)
    Phi = np.array(Phi_list)
    return Phi
'''
#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
def MOML(Theta,gamma): #we are ingnoring absorption radii outside the dust grain
    mu = 43
    z = 1
    x = Theta/z - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*np.exp(Theta/z)) + 0.5*np.log((((z**2)/(mu**2))*2*sp.pi*(1+gamma*Theta)))
    return x

def FlowingSheathMOML(Theta,kappa):
    mu = 43
    gamma = 5/3
    z = 1
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

Theta = np.logspace(-5,0,1001)
mu = 43
z = 1
Phi_fit = Phi_fit = ThinSheathFit(Theta)
#alpha,Phi = np.loadtxt("/Users/doganakpinar/Documents/Physics_Research/DustyPlasmaCollab/Plasma_code/Models/OMData3.txt",skiprows=1,unpack=True)
#Kappa, cov = spo.curve_fit(function,alpha,Phi,0.25)
#print(Kappa,cov)
kappa,cov_kappa = spo.curve_fit(FlowingSheathMOML,Theta,Phi_fit)
print(f"Best fit: kappa = {kappa} +/- {np.sqrt(cov_kappa[0,0])}")
gamma,cov_gamma = spo.curve_fit(MOML,Theta,Phi_fit)
print(f"Best fit: gamma = {gamma} +/- {np.sqrt(cov_gamma[0,0])}")
#plt.plot(alpha,function(alpha,Kappa))
#plt.grid()
#plt.show()