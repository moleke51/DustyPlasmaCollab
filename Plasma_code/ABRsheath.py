#============================IMPORT STUFF==========================#
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve,bisect


def get_name():
    return "ABR"

def colour():
    return 'blue'
    
#Define the normalised current J from equation ABR 12
def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*np.pi)**0.5))*np.exp(-Phi)
    return(j)
#Define equation ABR 13 for the boundary potential Phi_b.
def Eq_Phi_b(Phi_b, J, z=1,gamma = 10000):
    return 4 * (Phi_b**(3/2))*(2*Phi_b-3)*(2*Phi_b+1) - J/(gamma*np.sqrt(z)) * ((2*Phi_b-1)**3)
#Calculate the boundary conditions.
def get_boundary(J,z,gamma):
    Phi_b_initial_guess = 0.25
    Phi_b_solution = fsolve(Eq_Phi_b, Phi_b_initial_guess, args = (J,z,gamma))
    Phi_b = Phi_b_solution[0]
    #Calculate rho at the boundary. 
    rho_b = (np.sqrt(J) * np.exp(Phi_b/2))/((Phi_b*z)**(1/4))
    #Calculate the first derivative of Phi with renpect to rho evaluated at the boundary.
    dPhi_drho_b = ((2*rho_b/J) * (Phi_b**(3/2))/(Phi_b - 1/2) * np.exp(-Phi_b))*(z**0.5)
    return(rho_b,Phi_b,dPhi_drho_b)
#Runge-Kutta 4th order ODE solver
def ABR_f(x,y,z):
    return z
def ABR_g(x,y,z,A):
    return -(2/x)*z + (A/(x**2))*(y**(-0.5)) - np.exp(-y)
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

def Phi_J(J,alpha,mu,z,gamma = 10000):
    x0,y0,z0 = get_boundary(J,z,gamma) 
    Phi_b,Phi_alpha = ABR_RK(x0,y0,z0,alpha,J/np.sqrt(z),N=10000)
    return Phi_alpha

def delta_J(J,alpha,mu,z,gamma = 10000):
    return J - norm_J_current(alpha,Phi_J(J,alpha,mu,z,gamma = 10000),mu)
def retrive_Phi_a(J,mu,alpha):
    return 2*np.log(alpha) + np.log(mu) - np.log(J) - 0.5*np.log(4*np.pi)

def potential_finder(Theta,mu,z,alpha,upsilon,gamma=10000):
    #Guess Phi_a (Its likely to be between 0 and 10)
    if alpha != 0 :
        if alpha < 1e-5 and alpha > 0:
            Jsol = fsolve(delta_J,norm_J_current(alpha,0,mu),args = (alpha,mu,z,gamma))
        elif alpha > 1e12:
            Jsol = fsolve(delta_J,norm_J_current(alpha,-0.5*np.log(2*np.pi)+0.5+np.log(z*mu),mu),args = (alpha,mu,z,gamma))
        else:
            Jsol = bisect(delta_J,norm_J_current(alpha,0,mu),norm_J_current(alpha,-0.5*np.log(2*np.pi)+0.5+np.log(z*mu),mu),args = (alpha,mu,z,gamma))
        return retrive_Phi_a(Jsol,mu,alpha)
    else:
        return 0 #As argued by Kennedy and Allen
def presheath_pot(mu,z,alpha):
    Phi_a = potential_finder(0,mu,z,alpha,0)
    boundary = get_boundary(norm_J_current(alpha,Phi_a,mu),z,gamma = 10000)
    return(Phi_a - boundary[1],boundary[1],Phi_a)
def priority(Theta,alpha,upsilon):
    if Theta > 1e-4:
        P_t = 0
    else:
        P_t = 0
    P_a = 0
    if upsilon > 0:
        P_u = 0
    else:
        P_u = 0
    return (P_t + P_a + P_u)


Theta = 0
alpha = 1e14
mu = 43
z = 1
upsilon = 0
Phi = potential_finder(Theta,mu,z,alpha,upsilon)
print(Phi)
'''

alpha = np.logspace(10,20,11)
mu = 43*np.ones(len(alpha))
z = np.ones(len(alpha))
Phi_s_list = []
Phi_b_list = []
Phi_a_list = []
for i in range(len(alpha)):
    print(i)
    Phi_s,Phi_b,Phi_a = presheath_pot(mu[i],z[i],alpha[i])
    Phi_s_list.append(Phi_s)
    Phi_b_list.append(Phi_b)
    Phi_a_list.append(Phi_a)
Phi_s_arr = np.array(Phi_s_list)
Phi_b_arr = np.array(Phi_b_list)
Phi_a_arr = np.array(Phi_a_list)

plt.plot(alpha,Phi_s_arr,label = 'Phi_a - Phi_b',color = 'purple')
plt.plot(alpha,Phi_b_arr,label = 'Phi_b',color = 'blue')
plt.plot(alpha,Phi_a_arr,label = 'Phi_a',color = 'red')
plt.grid()
plt.title('The potential drop from the far field plasma to the sheath edge (across the presheath)')
plt.xlabel('alpha')
plt.ylabel('|$\Phi_s$|')
plt.xscale('log')
plt.show()

'''