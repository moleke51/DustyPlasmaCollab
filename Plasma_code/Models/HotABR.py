#============================IMPORT STUFF==========================#

import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve,bisect
import sys
'''
sys.path.append('.')
import Model as mdl


def get_Norm_Potential(Theta,mu,z,alpha,upsilon,variable_counter,previous_model = None,previous_phi = None):
    modellist = mdl.modelpicker('Models/',Theta,mu,z,alpha,upsilon)
    priority = 0
    for model in modellist:
        __import__(model.get_name())
        if model.priority() > priority:
            
            priority = model.priority()
            modelindex = modellist.index(model)
    
    if variable_counter == 0:
        print(modellist[modelindex])
        return modellist[modelindex].potential_finder()
    else:
        if modellist[modelindex].get_name() == 'ABR' and previous_model == 'ABR':
            return(previous_phi, previous_model)
        else:
            return(modellist[modelindex].potential_finder(), modellist[modelindex].get_name(),modellist[modelindex].get_colour())
'''

def get_name():
    return "HotABR"

def colour():
    return 'blue'
    

#Define the normalised current J from equation ABR 12
def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*np.pi)**0.5))*np.exp(-Phi)
    return(j)
#Define equation ABR 13 for the boundary potential tau + z*Phi_b.
def Eq_Phi_b(Phi_b, J, z=1,gamma = 10000):
    return 4 * (Phi_b**(3/2))*(2*Phi_b-3*z)*(2*Phi_b+z) - J/(gamma) * ((2*Phi_b-z)**3)
#Calculate the boundary conditions.
def get_boundary(J,tau,z,gamma):
    P_b_initial_guess = 0.25
    P_b_solution = fsolve(Eq_Phi_b, P_b_initial_guess, args = (J,z,gamma))
    P_b = P_b_solution[0]
    Phi_b = (P_b-tau)/z
    #Calculate rho at the boundary. 
    rho_b = (np.sqrt(J) * np.exp(Phi_b/2))/((P_b)**(1/4))
    #Calculate the first derivative of Phi with renpect to rho evaluated at the boundary.
    dPhi_drho_b = (((2*rho_b/J) * (P_b**(3/2))/(P_b - 0.5*z)) * np.exp(-Phi_b))
    return(rho_b,Phi_b,dPhi_drho_b)
#Runge-Kutta 4th order ODE solver
def ABR_f(x,y,w):
    return w
def ABR_g(x,y,w,J,tau,z):
    return -(2/x)*w + (J/(x**2))*(tau + z*(y**(-0.5))) - np.exp(-y)
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

def Phi_J(J,alpha,mu,tau,z,gamma = 10000):
    x0,y0,w0 = get_boundary(J,tau,z,gamma) 
    Phi_b,Phi_alpha = ABR_RK(x0,y0,w0,alpha,J,tau,z,N=100000)
    return Phi_alpha

def delta_J(J,alpha,mu,tau,z,gamma = 10000):
    return J - norm_J_current(alpha,Phi_J(J,alpha,mu,tau,z,gamma = 10000),mu)
def delta_Phi(Phi_a,alpha,mu,tau,z,gamma = 10000):
    return Phi_a - Phi_J(norm_J_current(alpha,Phi_a,mu),alpha,mu,tau,z,gamma = 10000)
def retrive_Phi_a(J,mu,alpha):
    return 2*np.log(alpha) + np.log(mu) - np.log(J) - 0.5*np.log(4*np.pi)

def potential_finder(Theta,mu,z,alpha,upsilon,kappa = 0.5,gamma=10000):
    tau = kappa*Theta
    #Guess Phi_a (Its likely to be between 0 and 10)
    if Theta == 0:
        if alpha != 0 :
            if alpha < 1e-5 and alpha > 0:
                Jsol = fsolve(delta_J,norm_J_current(alpha,0,mu),args = (alpha,mu,tau,z,gamma))
            elif alpha > 1e12:
                Jsol = fsolve(delta_J,norm_J_current(alpha,-0.5*np.log(2*np.pi)+0.5+np.log(z*mu),mu),args = (alpha,mu,tau,z,gamma))
            else:
                Jsol = bisect(delta_J,norm_J_current(alpha,0,mu),norm_J_current(alpha,-0.5*np.log(2*np.pi)+0.5+np.log(z*mu),mu),args = (alpha,mu,tau,z,gamma))
            return retrive_Phi_a(Jsol,mu,alpha)
        else:
            return 0 #As argued by Kennedy and Allen
    else:
        Phi_guess = potential_finder(0,mu,z,alpha,upsilon)
        print(Phi_guess)
        #Jsol = fsolve(delta_J,norm_J_current(alpha,Phi_guess,mu),args = (alpha,mu,tau,z,gamma))[0]
        #Jsol = bisect(delta_J,norm_J_current(alpha,2.0,mu),norm_J_current(alpha,5,mu),args = (alpha,mu,tau,z,gamma))
        #return retrive_Phi_a(Jsol,mu,alpha)
        return(fsolve(delta_Phi,Phi_guess,args= (alpha,mu,tau,z,gamma))[0])
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

print(potential_finder(0.01,43,1,100,0,2))


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
