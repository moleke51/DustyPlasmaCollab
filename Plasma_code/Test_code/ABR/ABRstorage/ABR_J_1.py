import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.special import erf
from scipy.optimize import fsolve

m_i = 1
gamma = 10000
F = 1
mu = 43
alpha_lim = 0.00000000001
z = 1
m_e = 9.11e-31 #kg



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

def plasma_criterion(phi_b):

    x = (4*(phi_b)**(3/2)*((2*phi_b)-3)*((2*phi_b)+1))/((2*phi_b)-1)**3

    return x

##Define the differential equation ABR 9. 
def JDE(H,rho, J, F):
    return( F*H[1], (-2/rho * H[1]) + (J/(sp.sqrt(H[0]) * rho**2)) -sp.exp(-H[0]))

#Define equation ABR 13 for the boundary potential Phi_b.   
def boundary_func(Phi_b, J, gamma):
    
    sol = (4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3) - J/gamma 
    
    return(sol)

def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*sp.pi)**0.5))*sp.exp(-Phi)
    return(j)

#Define a function to solve equation ABR (9).
# Guess Phi_a 
def ABR_potential_finder(mu,alpha_lim,z,J,iter,gamma = 10000, resolution = 0.0001):

    m_e = 9.11e-31 #kg
    Nu = mu/sp.sqrt((4*sp.pi))
    F=1
    Phi_b_initial_guess = 0.25
    Phi_b_solution = fsolve(boundary_func, Phi_b_initial_guess, args = (J, gamma))
    Phi_b = Phi_b_solution[0]
    rho_b = (sp.sqrt(J) * sp.exp(Phi_b/2))/(Phi_b**(1/4))
    #print(rho_b)
    dPhi_drho_b = (2*rho_b/J) * (Phi_b**(3/2))/(Phi_b - 1/2) * sp.exp(-Phi_b)
    H0 = [Phi_b,dPhi_drho_b]
    rhos = sp.linspace(rho_b,alpha_lim,10000000)
    Phis = integrate.odeint(JDE,H0,rhos,args=(J,F))
    Phi_sol = Phis[:,0]
    Phi_a = sp.log(((rhos**2)*Nu)/J)
    for i in range (0,len(Phi_a)):
        if sp.absolute((Phi_a[i] - Phi_sol[i])/Phi_a[i]) < resolution[iter]:
            if i != 0:
                turning_point = sp.absolute(((norm_J_current(rhos[i],Phi_sol[i],mu) - J)//1) - ((norm_J_current(rhos[i-1],Phi_sol[i-1],mu) - J)//1))
                if turning_point:
                    alpha = rhos[i]
                    Phi = Phi_sol[i]
                    print(str.format('{0:.4f}', alpha), str.format('{0:.4f}', Phi) )
                    return(alpha,Phi)


J_group = sp.logspace(-7,7,100)
alpha = sp.zeros(len(J_group))
Phi_a = sp.zeros(len(J_group))
res = sp.zeros(len(J_group))
for i in range(0, len(J_group)):
    if J_group[i] < 1e-5 :
        res[i] = 1
    elif J_group[i] < 1e-3:
        res[i] = 0.1
    else:
        res[i] = 0.0001
for i in range(0,len(J_group)):
    alpha[i],Phi_a[i] = ABR_potential_finder(mu,alpha_lim,z,J_group[i], i ,gamma,res)

with open("Hdata.txt","w") as f:
    for (distance,potential) in zip(alpha,Phi_a):
        f.write(f"{distance},{potential}\n")

plt.plot(alpha, Phi_a, label = 'Hydrogen plasma')
plt.title('ABR')
plt.ylabel('Normalised surface potential')
plt.xlabel('Normalised dust grain radius')
plt.xscale('log')
plt.grid()
plt.legend()
plt.show()



