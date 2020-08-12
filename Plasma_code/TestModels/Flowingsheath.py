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
def MOML_function(Phi,Theta,mu,z,alpha): #gamma = 5/3 for static plasmas
    return (np.sqrt(Theta)/mu)*(1 - (1/Theta)*(Phi - 0.5*(np.log(2*np.pi*(1+(5/3)*Theta))-np.log(mu**2)))) - np.exp(Phi)

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

def MOML_ps(Theta,gamma): 
    mu = 43
    z = 1
    return Theta - realLambertW((np.sqrt(2*np.pi*Theta*(1+gamma*Theta)))*np.exp(Theta)) 

#Parameters
Theta = np.logspace(-3,3,1001)
Theta2 = np.logspace(-7,3,1001)
mu = 43
z = 1

#Calculate MOML values
Phi_moml = np.zeros(len(Theta)) #gamma = 5/3
for i in range(len(Theta)):
    Phi_moml[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,5/3))

Phi_moml_ps = MOML_ps(Theta2,5/3)

'''
fig, (ax1, ax2) = plt.subplots(2, figsize = (7,7))
ax1.plot(Theta,Phi_moml, color = "red", label = "MOML")
ax1.plot(Theta, -3.34*np.ones(len(Theta)),'--', color = "red", label = "ABR value")
ax2.plot(Theta2,Phi_moml_ps, color = "red")
ax1.set(title = "(a) Normalised potential as a function of " + r"$\Theta$")
ax2.set(title = "(b) MOML pre-sheath potential drop as a function of " + r"$\Theta$")
ax1.set(ylabel = "Normalised potential, " + r"$\Phi_a$")
ax2.set(ylabel = "Normalised sheath edge potential, " + r"$\Phi_p$")
ax1.grid()
ax2.grid()
ax2.set(xlabel = r"$\Theta$")
ax1.set(xscale = "log")
ax2.set(xscale = "log")
ax1.legend()
plt.show()

'''


#Find best fit kappa
Phi_fit = ThinSheathFit(Theta)
kappa,cov_kappa = spo.curve_fit(FlowingSheathMOML,Theta,Phi_fit)
kappa = kappa[0]
print(f"Best fit: kappa = {kappa} +/- {np.sqrt(cov_kappa[0,0])}")


#Find the flowing sheath values
Phi_flow = FlowingSheathMOML(Theta,(2)*np.ones(len(Theta)))

'''
plt.plot(Theta,Phi_flow, color = "red")
plt.title("Normalised potential as a function of " + r"$\Theta$" + " for " + r"$\kappa$" + " = 2")
plt.ylabel("Normalised potential, " + r"$\Phi_a$")
plt.xlabel(r"$\Theta$")
plt.grid()
plt.xscale("log")
'''

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
fig, (ax1, ax2) = plt.subplots(2, sharex= True)
ax1.plot(Theta,Phi_moml,color = "Black", label = "MOML, " + r"$\gamma$" + "= 5/3")
ax1.plot(Theta,Phi_fit,"--",color = 'Blue',label = 'SCEPTIC Fit')
ax1.plot(Theta,Phi_flow,color = 'red',label = "Flowing sheath, " + r"$\gamma$" + "= 5/3, " + r"$\kappa$" + " = 2")
ax1.plot(Theta, -3.34*np.ones(len(Theta)), ':', color = "Purple", label = "ABR value")
ax1.set(title = "(a)")
ax1.set(ylabel = "Normalised potential, " + r"$\Phi_a$")
ax1.set(xscale = "log")
ax1.legend(loc = 'upper right')
ax1.grid()


ax2.plot(Theta,Phi_mfit,'--',color = 'Black',label = "MOML")
plt.plot(Theta,Phi_ffit,color = 'Red',label = "Flowing sheath")
plt.xlabel(r"$\Theta$")
ax2.set(title = "(b)")
ax2.set(ylabel = "Percentage difference")
ax2.legend()
ax2.grid()
ax2.set(xscale = "log")

plt.show()


#a = 0
#while a == 0:
#    print(FlowingSheathMOML(float(input('Select theta value: ')),43,5/3,2))
