import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve,bisect


def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*np.pi)**0.5))*np.exp(-Phi)
    return(j)
#Define equation ABR 13 for the boundary potential tau + z*Phi_b.
def Eq_Phi_b(Phi_b, J, z=1,gamma = 10000):
    return 4 * (Phi_b**(3/2))*(2*Phi_b-3*z)*(2*Phi_b+z) - J/(gamma) * ((2*Phi_b-z)**3)
#Calculate the boundary conditions.
def get_boundary(J,tau,z,gamma = 10000):
    P_b_initial_guess = 0.25
    P_b_solution = fsolve(Eq_Phi_b, P_b_initial_guess, args = (J,z,gamma))
    P_b = P_b_solution[0]
    Phi_b = (P_b-tau)/z
    if Phi_b < 0:
        Phi_b = 0
        P_b = tau
    #Calculate rho at the boundary. 
    rho_b = (np.sqrt(J) * np.exp(Phi_b/2))/((P_b)**(1/4))
    #Calculate the first derivative of Phi with renpect to rho evaluated at the boundary.
    dPhi_drho_b = (((2*rho_b/J) * (P_b**(3/2))/(P_b - 0.5*z)) * np.exp(-Phi_b))*(z**0.5)
    return(rho_b,Phi_b,dPhi_drho_b)
#Runge-Kutta 4th order ODE solver
def ABR_f(x,y,w):
    return w
def ABR_g(x,y,w,J,tau,z):
    if y != 0:
        return -(2/x)*w + (J/(x**2))*(tau + z*(y**(-0.5))) - np.exp(-y)
    else:
        return -(2/x)*w + (J/(x**2))*(tau) - 1
def ABR_RK(x0,y0,w0,X,J,tau,z,N):
    h = (X-x0)/N
    x = X - N*h
    y = y0
    w = w0
    x_arr = [x]
    y_arr = [y]
    w_arr = [w]

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
        x_arr.append(x) 
        y_arr.append(y)
        w_arr.append(w)    
    return(y0,y,x_arr,y_arr,w_arr)

J = 5e3
tau = 0.02
z = 1
alpha = 10
x0,y0,w0 = get_boundary(J,tau,z) 
Phi_b,Phi_alpha,X,Y,W = ABR_RK(x0,0,w0,alpha,J,tau,z,N=100000)
print(X[0],Y[0])
'''
for i in range(len(X)):
    print(X[i],Y[i])
'''
plt.figure(1)
plt.plot(X,Y,color = 'Blue',label = 'Function')
plt.xlabel('Distance')
plt.ylabel('Potential')
plt.grid()
plt.legend()
plt.figure(2)
plt.plot(X,W,color = 'Red',label = 'Derivative')
plt.xlabel('Distance')
plt.ylabel('Derivative of potential')
plt.grid()
plt.legend()
plt.show()