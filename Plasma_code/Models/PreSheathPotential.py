import numpy as np
import matplotlib.pyplot as plt 

'''
def PSD(Theta,mu,z,rho):
    gamma = 1
    Phi_hot = np.absolute(0.5*np.log(2*np.pi*(1+gamma*Theta)) + np.log(z/mu))
    rho_hot = (np.sqrt(2)/3)*(2*Phi_hot)**(3/4)
    Phi_cold = np.absolute(0.5*np.log(2*np.pi) + np.log(z/mu))
    rho_cold = (np.sqrt(2)/3)*(2*Phi_cold)**(3/4)

    A = (Phi_hot/(2*Phi_hot - 1))*(rho_cold - rho_hot)
    C = (1/(2*Phi_hot - 1))*(rho_cold - 2*Phi_hot*rho_hot)
    return -A/(rho+C) - 0.5
'''
'''
def PSD(Theta,mu,z):
    Phi_cold = np.absolute(np.log(2*np.pi) + np.log(z/mu))
    C = -1*(np.sqrt(2)/3)*Phi_cold**(3/4)
    A = (np.sqrt(2)/6)*((Phi_cold)**(3/4))*(1-(Phi_cold)**(1/4))
    gamma = 5/3
    Phi_hot = np.absolute(0.5*np.log(2*np.pi*(1+gamma*Theta)) + np.log(z/mu))
    rho = (np.sqrt(2)/3)*(2*Phi_hot)**(3/4)
    return -A/(rho+C) 

'''

def PSD(Theta,mu,z):
    gamma = 1
    Phi_hot = np.absolute(0.5*np.log(2*np.pi*(1+gamma*Theta)) + np.log(z/mu))
    rho_hot = (np.sqrt(2)/3)*(2*Phi_hot)**(3/4)
    Phi_cold = np.absolute(0.5*np.log(2*np.pi) + np.log(z/mu))
    rho_cold = (np.sqrt(2)/3)*(2*Phi_cold)**(3/4)
    A = 0.5*rho_cold
    return -A/rho_hot 

Theta = np.linspace(0,10,101)
mu = 43
z = 1
#rho = np.linspace(0,10,101)

'''
gamma = 5/3
Phi_hot = np.absolute(0.5*np.log(2*np.pi*(1+gamma*Theta)) + np.log(z/mu))
rho_hot = (np.sqrt(2)/3)*(2*Phi_hot)**(3/4)
Phi_range = np.linspace(-3,3,101)
rho_hot = (np.sqrt(2)/3)*(2*Phi_hot)**(3/4)*np.ones(len(Phi_range))
'''

Phi = PSD(Theta,mu,z)
plt.plot(Theta,Phi)
#plt.plot(rho_hot,Phi_range,'--')
plt.show()