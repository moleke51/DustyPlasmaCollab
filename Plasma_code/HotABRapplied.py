#============================IMPORT STUFF==========================#

import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve,bisect

def get_name():
    return "HotABR"

def colour():
    return 'blue'
    

#Define the normalised current J from equation ABR 12
def norm_J_current(alpha,Psi,tau,mu,z):
    j = (alpha**2)*(mu/((4*np.pi)**0.5))*np.exp((tau-Psi)/z)
    return(j)
#Define equation ABR 13 for the boundary potential tau + z*Phi_b.
def Eq_Psi_b(Psi_b, J, z=1,gamma = 10000):
    return 4 * (Psi_b**(3/2))*(2*Psi_b-3*z)*(2*Psi_b+z) - J/(gamma) * ((2*Psi_b-z)**3)
#Calculate the boundary conditions.
def get_boundary(J,tau,z,gamma):
    Psi_b_initial_guess = 0.25
    Psi_b_solution = fsolve(Eq_Psi_b, Psi_b_initial_guess, args = (J,z,gamma))
    Psi_b = Psi_b_solution[0]
    
    #Calculate rho at the boundary. 
    rho_b = (np.sqrt(J) * np.exp((Psi_b-tau)/(2*z)))/((Psi_b)**(1/4))
    #Calculate the first derivative of Phi with renpect to rho evaluated at the boundary.
    dPsi_drho_b = (((2*rho_b/J) * (Psi_b**(3/2))/(Psi_b - 0.5*z)) * np.exp((tau-Psi_b)/z))*z
    return(rho_b,Psi_b,dPsi_drho_b)
#Runge-Kutta 4th order ODE solver
def ABR_f(x,y,w):
    return w
def ABR_g(x,y,w,J,tau,z):
    return -(2/x)*w + (J*z/(x**2))*(y**(-0.5)) - z*np.exp(tau/z)*np.exp(-y/z)
    
def ABR_RK(x0,y0,w0,X,J,tau,z,N):
    h = (X-x0)/N
    x = X - N*h
    y = y0
    w = w0
    for n in range(N):
        k1 = ABR_f(x,y,w)
        l1 = ABR_g(x,y,w,J,tau,z)

        k2 = ABR_f(x+h/2,y+k1*h/2,w+l1*h/2)
        l2 = ABR_g(x+h/2,y+k1*h/2,w+l1*h/2,J,tau,z)

        k3 = ABR_f(x+h/2,y+k2*h/2,w+l2*h/2)
        l3 = ABR_g(x+h/2,y+k2*h/2,w+l2*h/2,J,tau,z)

        k4 = ABR_f(x+h,y+k3*h,w+l3*h)
        l4 = ABR_g(x+h,y+k3*h,w+l3*h,J,tau,z)

        x += h
        y += 1/6 * (k1 + 2*k2 + 2*k3 + k4)*h
        w += 1/6 * (l1 + 2*l2 + 2*l3 + l4)*h
                
    return(y0,y)

def Psi_J(J,alpha,mu,tau,z,gamma = 10000):
    x0,y0,w0 = get_boundary(J,tau,z,gamma) 
    Psi_b,Psi_alpha = ABR_RK(x0,y0,w0,alpha,J,tau,z,N=100000)
    return Psi_alpha
def delta_J(J,alpha,mu,tau,z,gamma = 10000):
    return J - norm_J_current(alpha,Psi_J(J,alpha,mu,tau,z,gamma = 10000),tau,mu,z)
def delta_Psi(Psi_a,alpha,mu,tau,z,gamma = 10000):
    return Psi_a - Psi_J(norm_J_current(alpha,Psi_a,tau,mu,z),alpha,mu,tau,z,gamma = 10000)
def retrive_Psi_a(J,mu,alpha,z,tau):
    return ((2*np.log(alpha) + np.log(mu) - np.log(J) - 0.5*np.log(4*np.pi))*(z))+tau

def potential_finder(Theta,mu,z,alpha,upsilon,kappa = 0.5,gamma=10000):
    tau = kappa*Theta
    #Guess Phi_a (Its likely to be between 0 and 10)
    if Theta == 0:
        if alpha != 0 :
            if alpha < 1e-5 and alpha > 0: #alpha,Psi,tau,mu,z
                Jsol = fsolve(delta_J,norm_J_current(alpha,0,tau,mu,z),args = (alpha,mu,tau,z,gamma))
            elif alpha > 1e12:
                Jsol = fsolve(delta_J,norm_J_current(alpha,-0.5*np.log(2*np.pi)+0.5+np.log(z*mu),tau,mu,z),args = (alpha,mu,tau,z,gamma))
            else:
                Jsol = bisect(delta_J,norm_J_current(alpha,0,tau,mu,z),norm_J_current(alpha,-0.5*np.log(2*np.pi)+0.5+np.log(z*mu),tau,mu,z),args = (alpha,mu,tau,z,gamma))
            return retrive_Psi_a(Jsol,mu,alpha,z,tau)
        else:
            return 0 #As argued by Kennedy and Allen
    else:
        #Jsol = fsolve(delta_J,norm_J_current(alpha,Phi_guess,mu),args = (alpha,mu,tau,z,gamma))[0]
        #Jsol = bisect(delta_J,norm_J_current(alpha,2.0,mu),norm_J_current(alpha,5,mu),args = (alpha,mu,tau,z,gamma))
        #return retrive_Phi_a(Jsol,mu,alpha)
        return((fsolve(delta_Psi,3,args= (alpha,mu,tau,z,gamma))[0]-tau)/z)
def priority(Theta,alpha,upsilon):
    if Theta > 1e-4:
        P_t = 0
    else:
        P_t = 1
    P_a = 1
    if upsilon > 0:
        P_u = 0
    else:
        P_u = 1
    return (P_t + P_a + P_u)

Theta = np.logspace(-5,0,6)
Phi = np.zeros(len(Theta))
for i in range(len(Theta)):
   Phi[i] = potential_finder(Theta[i],43,1,100,0,2)
   print(Theta[i],Phi[i])
plt.plot(Theta,Phi,label = 'Kappa = 2')
plt.xlabel("Theta")
plt.ylabel("Normalised potential")
plt.grid()
plt.xscale('log')
plt.legend()
plt.show()
'''
Philip = np.linspace(-100,100,10000)
Jay = norm_J_current(np.ones(len(Philip)),Philip,43*np.ones(len(Philip)))

Logan = 'log'
plt.plot(Philip,Jay)
plt.grid()
plt.xscale(Logan)
plt.yscale(Logan)
plt.show()
'''
