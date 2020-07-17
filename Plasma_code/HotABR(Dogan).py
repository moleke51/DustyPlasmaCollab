import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve,bisect

Theta = 0
gamma = 5/3
z = 1
upsilon = 0
alpha = 1.08*(50)
mu = 43

#Define the normalised current J from equation ABR 12
def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*np.pi)**0.5))*np.exp(-Phi)
    return(j)

def ABR_f(x,y,z):
    return z
def ABR_g(x,y,z,A):
    return -(2/x)*z + (A/(x**2)) - np.exp(-y)
def ABR_RK(x0,y0,z0,X,A,N):
    h = (X-x0)/N
    x = X - N*h
    y = y0
    z = z0

    for n in range(N):
        k1 = ABR_f(x,y,z)
        l1 = ABR_g(x,y,z,A)

        k2 = ABR_f(x+h/2,y+k1*h/2,z+l1*h/2)
        l2 = ABR_g(x+h/2,y+k1*h/2,z+l1*h/2,A)

        k3 = ABR_f(x+h/2,y+k2*h/2,z+l2*h/2)
        l3 = ABR_g(x+h/2,y+k2*h/2,z+l2*h/2,A)

        k4 = ABR_f(x+h,y+k3*h,z+l3*h)
        l4 = ABR_g(x+h,y+k3*h,z+l3*h,A)

        x += h
        y += 1/6 * (k1 + 2*k2 + 2*k3 + k4)*h
        z += 1/6 * (l1 + 2*l2 + 2*l3 + l4)*h
                
    return(y0,y)

def Phi_J(Theta,mu,z,alpha,gamma,J):
    #x0,y0,z0 = get_boundary(Theta,z,gamma,J) 
    x0 = np.sqrt(J/z)*(2/(1+gamma*Theta))**(1/4)
    y0 = 0
    z0 = (2*np.sqrt(z)/J)*((1+gamma*Theta)/2)**(1/4)
    Phi_b,Phi_alpha = ABR_RK(x0,y0,z0,alpha,(J/z)*np.sqrt(2/(1+gamma*Theta)),N=10000)
    print(Phi_alpha)
    return Phi_alpha

def delta_J(J,alpha,mu,z):
    return J - norm_J_current(alpha,Phi_J(Theta,mu,z,alpha,gamma,J),mu)
def retrive_Phi_a(J,mu,alpha):
    return 2*np.log(alpha) + np.log(mu) - np.log(J) - 0.5*np.log(4*np.pi)

def potential_finder(Theta,mu,z,alpha,upsilon,gamma):
    #Guess Phi_a (Its likely to be between 0 and 10)
    Jsol = bisect(delta_J,norm_J_current(alpha,0,mu),norm_J_current(alpha,10,mu),args = (alpha,mu,z))
    #Jsol = fsolve(delta_J,norm_J_current(alpha,0.5,mu),args = (alpha,mu,z))
    return retrive_Phi_a(Jsol,mu,alpha)

Phi = potential_finder(Theta,mu,z,alpha,upsilon,gamma)
print(Phi)