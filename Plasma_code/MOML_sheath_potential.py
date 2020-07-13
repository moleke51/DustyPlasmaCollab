import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve

def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64 :
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

def sheath_pot(Theta,gamma): #gamma = 5/3
    Phi = Theta - realLambertW((np.sqrt(2*np.pi*Theta*(1+(gamma)*Theta)))*np.exp(Theta)) 
    return Phi #returned phi is positive     

def presheath_drop(Theta,gamma):
    return -1*0.5 + (-1*gamma + 3/2)*Theta

def RI(Theta,gamma = 5/3):
    return (sheath_pot(Theta,gamma) - presheath_drop(Theta,gamma))

Theta = np.logspace(-9,2,1000)
gamma = 5/3*np.ones(len(Theta))
Phi_s = sheath_pot(Theta,gamma)
Phi_ps = presheath_drop(Theta,gamma)

Theta1 = fsolve(RI,0)
Theta2 = fsolve(RI,1e1)
print(Theta1, Theta2)

plt.plot(Theta,Phi_s,color = 'Black',label = "MOML prediction")
plt.plot(Theta,Phi_ps,'--',color = 'Red', label = "Derived prediction")
plt.plot(Theta1,presheath_drop(Theta1,5/3),'x',color = "Blue", label = f"Theta1 = {Theta1}")
plt.plot(Theta2,presheath_drop(Theta2,5/3),'x',color = "Blue", label = f"Theta2 = {Theta2}")
plt.grid()
plt.title('The potential drop across the presheath')
plt.xlabel('$\Theta$')
plt.ylabel('$\Phi$')
plt.xscale('log')
plt.legend()
plt.show()
