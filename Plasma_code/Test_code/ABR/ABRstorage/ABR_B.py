#Equations and information from "The floating potential of spherical probes and dust grains. 
# Part 1. Radial motion theory" by R. V. Kennedy and J. E. Allen shall be referenced with ABR.

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.special as sps
from scipy.optimize import fsolve

def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64:
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

#Define the differential equation ABR 9. 
def JDE(H,rho, J, F):
    return( F*H[1], (-2/rho * H[1]) + (J/(sp.sqrt(H[0]) * rho**2)) -sp.exp(-H[0]))


#B.12
def BDE(H,rho):
    return(H[1], (4*H[0]*(5 - 4*H[0]*(1-H[0]*(1 - 24*(rho**(-2)) - 324*(rho**(-4))   ))))/ (((1-2*H[0]*(1 - 18*(rho**(-2))))**3) * (rho**2)))


#Define equation ABR 13 for the boundary potential Phi_b.   
def boundary_func(Phi_b, J, gamma):
    sol = (4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3) - J/gamma 
    return(sol)

#Define a function to solve equation ABR 9.
# Guess Phi_a 
def ABR_potential_finder(mu,alpha,z,Phi_a, gamma = 10000):
    
    nu = alpha**2 * mu /sp.sqrt(4*sp.pi)
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
    #rho_b = (sp.sqrt(J) * sp.exp(Phi_b/2))/(Phi_b**(1/4))
    
    #B.8 
    rho_b = sp.sqrt(sp.exp(Phi_b)/sp.sqrt(Phi_b)*(-12*Phi_b*sp.sqrt(Phi_b) + nu*sp.exp(-Phi_a)))
    #Find dPhi/drho at the boundary with equation ABR 11.
    #dPhi_drho_b = (2*rho_b/J) * (Phi_b**(3/2))/(Phi_b - 1/2) * sp.exp(-Phi_b)

    #B.11
    dPhi_drho_b = (-4*Phi_b*sp.exp(-Phi_b)*rho_b) / ( (sp.exp(-Phi_b) * (rho_b**2 -2*Phi_b*(rho_b**2 -6))) +24*Phi_b )

    #Collect the boundary conditions.
    H0 = [Phi_b,dPhi_drho_b]
    rhos = sp.linspace(rho_b,alpha,1000)
    #Solve equation ABR 9 for the surface potential Phi_a
    #Phis = integrate.odeint(JDE,H0,rhos,args=(J,F))
    Phis = integrate.odeint(BDE,H0,rhos)
    Phi_sol = Phis[-1][0]
    #If the solution is equal to the initial guess return the negative of that value, 
    #as the function is defined for positive potential but the value is always negative.
    if sp.absolute((Phi_a - Phi_sol)/Phi_a) < 0.0001:
        return(-Phi_sol)
    #Otherwise make a new guess at Phi_a equal to the mean of the solution and the original guess.
    else:
        Phi_a = (Phi_a + Phi_sol)/2
        return(ABR_potential_finder(mu,alpha,z,Phi_a, gamma))
        
mu = 43 #mu values of above 63, do not work with alpha = 1
alpha = 1 #For mu values of 43, the function doesn't work with alphas greater than 1 
z = 1 #z doesn't affect anything yet
Phi_a = 3 #Reasonable guesses for Phi_a are between 1 and 10
Phi_sol = ABR_potential_finder(mu,alpha,z,Phi_a)
print(Phi_sol)


