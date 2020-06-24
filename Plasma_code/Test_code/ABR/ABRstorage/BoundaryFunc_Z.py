#============================IMPORT STUFF============================

import scipy as sp
import matplotlib.pyplot as plt
import periodictable as pt
import numpy as np
from scipy.special import lambertw
from scipy.optimize import fsolve


#============================FUNCTIONS============================

J = sp.logspace(-7,10,100)
#alpha = sp.zeros(len(J))
#Phi_a = sp.zeros(len(J))
#res = sp.zeros(len(J))
gamma = 10000



def boundary_func(Phi_b, J, z,gamma = 10000):
    
    sol = sp.sqrt(z)*((4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3)) - J/gamma 
    
    return(sol)
Phi_b_initial_guess = 0.25
Phi_b = sp.zeros(len(J))
for i in range(len(J)):
    Phi_b_solution = fsolve(boundary_func, Phi_b_initial_guess, args = (J[i], gamma))
    Phi_b[i] = Phi_b_solution[0]
print(Phi_b)
plt.plot(J, Phi_b)
plt.title('Phi_b against J')
plt.ylabel('Phi_b')
plt.xlabel('J')
plt.legend()
plt.xscale('log')
plt.grid()
plt.show()