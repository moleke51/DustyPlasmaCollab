#============================IMPORT STUFF============================

import scipy as sp
import matplotlib.pyplot as plt
import periodictable as pt
import numpy as np
from scipy.special import lambertw
from scipy.optimize import fsolve


#============================FUNCTIONS============================

J = sp.logspace(-7,13,1000)
#alpha = sp.zeros(len(J))
#Phi_a = sp.zeros(len(J))
#res = sp.zeros(len(J))
gamma = 10000
z=1


def boundary_func(Phi_b, J, z=1,gamma = 10000):
    
    sol = sp.sqrt(z)*((4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3)) - J/gamma 
    
    return(sol)

Phi_b_initial_guess = sp.zeros(len(J))
for i in range(len(J)):
    #if J[i]/gamma > 1e3:
        #Phi_b_initial_guess[i] = 0.45
    if J[i]/gamma > 1e2:
        Phi_b_initial_guess[i] = 0.5
    if J[i]/gamma > 1e1:
        Phi_b_initial_guess[i] = 0.30
    else:
        Phi_b_initial_guess[i] = 0.25
    
Phi_b = sp.zeros(len(J))
for i in range(len(J)):
    Phi_b_solution = fsolve(boundary_func, Phi_b_initial_guess[i], args = (J[i],z,gamma))
    Phi_b[i] = Phi_b_solution[0]
#print(Phi_b)
plt.plot(Phi_b,J/gamma)
plt.title('Phi_b against J')
plt.xlabel('Phi_b')
plt.ylabel('J')
#plt.legend()
plt.yscale('log')
plt.grid()
plt.show()