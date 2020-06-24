#Equations and information from "The floating potential of spherical probes and dust grains. 
# Part 1. Radial motion theory" by R. V. Kennedy and J. E. Allen shall be referenced with ABR.

import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import fsolve

#Define the differential equation ABR 9. 
def JDE(H,rho, J, F):
        return( F*H[1], (-2/rho * H[1]) + (J/(sp.sqrt(H[0]) * rho**2)) -sp.exp(-H[0]))

#Define equation ABR 13 for the boundary potential Phi_b.   
def boundary_func(Phi_b, J, gamma):
    sol = (4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3) - J/gamma 
    return(sol)

#Define a function to solve equation ABR 9.
# Guess Phi_a 
def ABR_potential_finder(mu,alpha,z,Phi_a, gamma = 10000):
    #10000 seems to be a reasonable value for gamma as it must be significantly greater than unity.
    #Find J using the guess for Phi_a, the normalised dust radius and equation ABR 12.
    J = alpha**2 * mu * 1/(sp.sqrt(4*sp.pi)) * sp.exp(-Phi_a)
    #For some unknown reason odeint requires two parameters so F is used, do not change from F = 1.
    F=1
    #Phi_b is always between 0 and 0.5 so an initial guess of 0.25 is always valid.
    Phi_b_initial_guess = 0.25
    #Solve for Phi_b using equation ABR 13.
    Phi_b_solution = fsolve(boundary_func, Phi_b_initial_guess, args = (J, gamma))
    Phi_b = Phi_b_solution[0]
    #Solve for rho_b using equation ABR 10
    rho_b = (sp.sqrt(J) * sp.exp(Phi_b/2))/(Phi_b**(1/4))
    #Find dPhi/drho at the boundary with equation ABR 11.
    dPhi_drho_b = (2*rho_b/J) * (Phi_b**(3/2))/(Phi_b - 1/2) * sp.exp(-Phi_b)
    #Collect the boundary conditions.
    H0 = [Phi_b,dPhi_drho_b]
    rhos = sp.linspace(rho_b,alpha,1000)
    #Solve equation ABR 9 for the surface potential Phi_a
    Phis = integrate.odeint(JDE,H0,rhos,args=(J,F))
    Phi_sol = Phis[-1][0]
    #If the solution is equal to the initial guess return the negative of that value, 
    #as the function is defined for positive potential but the value is always negative.
    if sp.absolute((Phi_a - Phi_sol)/Phi_a) < 0.0001:
        return(-Phi_sol)
    #Otherwise make a new guess at Phi_a equal to the mean of the solution and the original guess.
    else:
        Phi_a = (Phi_a + Phi_sol)/2
        return(ABR_potential_finder(mu,alpha,z,Phi_a, gamma))
        
mu = 271 #mu values of above 63, do not work with alpha = 1
alpha = 1 #For mu values of 43, the function doesn't work with alphas greater than 1 
z = 1 #z doesn't affect anything yet
Phi_a = 3 #Reasonable guesses for Phi_a are between 1 and 10
Phi_sol = ABR_potential_finder(mu,alpha,z,Phi_a)
print(Phi_sol)


