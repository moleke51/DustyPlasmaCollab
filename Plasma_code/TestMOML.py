import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt 
theta = 1
mu = 43
gamma = 5/3
def norm_ion_current(eta,theta,mu,gamma):
    return mu*np.exp(theta-mu*np.sqrt(theta)*np.exp(-eta))*np.sqrt(1+gamma*theta)
def norm_electron_current(eta):
    return (1/np.sqrt(2*np.pi))*np.exp(-eta)

eta_a = np.logspace(-9,0.3,1000)
theta = theta*np.ones(len(eta_a))
mu = mu*np.ones(len(eta_a))
gamma = gamma*np.ones(len(eta_a))
J_i = norm_ion_current(eta_a,theta,mu,gamma)
J_e = norm_electron_current(eta_a)
plt.plot(eta_a,J_i,label = 'Ion current',color = 'red')
plt.plot(eta_a,J_e,label = 'Electron current',color = 'blue')
plt.grid()
plt.legend()
plt.xscale('log')
plt.show()