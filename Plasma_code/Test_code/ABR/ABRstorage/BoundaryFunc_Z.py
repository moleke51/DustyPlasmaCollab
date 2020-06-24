#============================IMPORT STUFF==========================#

import scipy as sp
import matplotlib.pyplot as plt
import periodictable as pt
import numpy as np
from scipy.optimize import fsolve


#============================FUNCTIONS============================#

J = sp.logspace(-6,14,1000)
gamma = 10000
z=1

def boundary(Phi_b, J, z=1,gamma = 10000):
    return 4 * (Phi_b**(3/2))*(2*Phi_b-3)*(2*Phi_b+1) - J/(gamma*np.sqrt(z)) * ((2*Phi_b-1)**3)


Phi_b_initial_guess = 0.25
Phi_b = sp.zeros(len(J))
for i in range(len(J)):
    Phi_b_solution = fsolve(boundary, Phi_b_initial_guess, args = (J[i],z,gamma))
    Phi_b[i] = Phi_b_solution[0]

plt.plot(Phi_b,J/gamma)
plt.title('J/$\gamma$ against Phi_b')
plt.xlabel('Phi_b')
plt.ylabel('J/$\gamma$')
plt.yscale('log')
plt.grid()
plt.show()