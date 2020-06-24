import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.special import erf

#This function gives the real value of the Lambert W0 function provided x > -1/e
def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64 :
        w = sps.lambertw(x)
        if sp.imag(w)==0:
            W = sp.real(w)
            return(W)
        else:
            return('This value is outside the accepted range')
    elif type(x) is np.ndarray or type(x) is list:
        W = sp.zeros(len(x))
        for i in range(0,len(x)):
            
            if sp.imag(sps.lambertw(x[i]))==0:
                W[i]= sp.real(sps.lambertw(x[i]))
            else:
                print('The value at position ' +str(i)+' is outside the accepted range')
        return(W)
    else:
        return('This is an invalid input')

def dh_drho(h,rho):

    #h = [Phi, Phi_d1] 

    #h_d = [Phi_d1,Phi_d2]

    return [h[1], sp.exp(h[0]) - (nu*sp.exp(h[0]))/(rho*sp.sqrt(-h[0])) - (1/rho)*h[1]]

mu = 43
alpha = 1
rho_b = 50
z = 1
rhos = sp.linspace(rho_b,alpha,100)
Phi_a = -13

nu = ((alpha)**2)*((mu)/(sp.sqrt(4*sp.pi*z)))

Phi_b = 0.5*realLambertW((-2*(nu**2)*sp.exp(2*Phi_a))/((rho_b)**4))

h0 = [Phi_b,(-4*Phi_b)/((1+2*Phi_b)*rho_b)]

Phis = integrate.odeint(dh_drho,h0,rhos)
print(Phis)
x = Phis[:,0]*0 +1
print(x)
print(Phis[-1][0])
#plt.plot(rhos,Phis[:,1],'x',label = 'Phi')
plt.plot(rhos,Phis[:,0],'+',label = 'Phi')
plt.plot(x,Phis[:,0])
plt.grid()
plt.legend()
plt.show()



