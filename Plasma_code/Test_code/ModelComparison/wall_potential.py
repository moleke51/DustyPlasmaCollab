import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt 

e = 1.6 * 10**(-19)
k_B = 1.38 *10**(-23)

def wall_potential(mu,T_e):
    phi = (1/2) * (sp.log(2*sp.pi) - 2*sp.log(mu)) - 0.5
    return(sp.absolute(phi))
def ABR_limit(mu):
    return(-(0.4189 - sp.log(mu)))
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

def MOML_mu(Theta,mu,z,gamma): #we are ingnoring absorption radii outside the dust grain

    x = Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2))

    return(-x )
T_e = 1*e/k_B

mu = sp.linspace(43,1000,10000)



phi= wall_potential(mu,T_e)
plt.plot(mu,phi,color = 'black' ,label = 'Wall potential')

lims = ABR_limit(mu)
plt.plot(mu,lims,color = 'black', label = 'ABR limits') 

Theta = sp.logspace(-6,0,7)
#colour = ['purple','blue','green']
for i in range(0,len(Theta)):
    moml = MOML_mu(Theta[i],mu,1,5/3)
    plt.plot(mu, moml ,label = 'MOML theta = ' + str(Theta[i]))
plt.xlabel('$\mu$')
plt.ylabel('$\Phi$ ')
plt.grid()
plt.legend()
plt.show()