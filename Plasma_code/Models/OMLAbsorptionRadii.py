import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve,bisect
import scipy.optimize as spo
import OML as oml

def OML_AR(Phi,Theta,mu,z,alpha,rho):
    return (rho**2/alpha**2)*Theta - z*(rho/alpha)*np.exp(rho-alpha)*Phi - mu*np.sqrt(Theta)*np.exp(Phi)

alpha = np.logspace(-3,3,101)
Theta = 1*np.ones(len(alpha))
mu = 43*np.ones(len(alpha))
z = 1*np.ones(len(alpha))
rho = 1.08*alpha
upsilon = 0*np.ones(len(alpha))

Phi_list = []
for i in range(len(alpha)):
    Phi = bisect(OML_AR,-10,10,args = (Theta[i],mu[i],z[i],alpha[i],rho[i]))
    Phi_list.append(Phi)
Phi = np.array(Phi_list)

Phi_oml_list = []
for i in range(len(alpha)):
    Phi_oml = oml.potential_finder(Theta[i],mu[i],z[i],alpha[i],upsilon[i])
    Phi_oml_list.append(Phi_oml)
Phi_oml = -1*np.array(Phi_oml_list)

plt.plot(alpha,Phi,color = "Red", label = "New OML")
plt.plot(alpha,Phi_oml,color = "Black", label = "OML")
plt.xscale("log")
plt.title("Normalised dust potential as a function of alpha")
plt.xlabel("alpha")
plt.ylabel("Normalised dust potential")
plt.grid()
plt.legend()
plt.show()