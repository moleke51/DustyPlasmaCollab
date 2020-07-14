import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve

Theta = np.logspace(0,3,101)
#Theta = np.linspace(0,3,101)
gamma = (5/3)*np.ones(len(Theta))

#def pre_sheath(Phi,Theta,gamma):
#    return ((np.exp(-2*Phi) - 1) + 2*Phi/(1+gamma*Theta))
#def grapher(Theta, gamma):
#   Phi_list = []
#    for i in range(len(Theta)):
#        Phi = fsolve(pre_sheath,0.5,args = (Theta[i],gamma[i]))
#        Phi_list.append(Phi)
 #   Phi = np.array(Phi_list)
 #   return Phi

def potential(Theta):
    return (-1*0.5*Theta)/(Theta + 1)

#def potential(Theta):
#    return (-1*0.5/(1+Theta))

Phi = potential(Theta)
plt.plot(Theta,Phi)
plt.xscale('log')
plt.grid()
plt.show()