import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt 

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
    return np.absolute(Phi) #returned phi is positive     


theta = np.logspace(-9,2,1000)
gamma = 5/3*np.ones(len(theta))
Phi_s = sheath_pot(theta,gamma)
plt.plot(theta,Phi_s,color = 'blue')
plt.grid()
plt.title('The potential drop from the far field plasma to the sheath edge')
plt.xlabel('$\Theta$')
plt.ylabel('|$\Phi_s$|')
plt.xscale('log')
plt.show()
