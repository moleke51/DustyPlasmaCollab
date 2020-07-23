#============================IMPORT STUFF==========================#
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve,bisect
import OML as oml   
import MOML as moml
import SNF as snf
import HotABR as abr

alpha = np.logspace(-3,3,50).tolist()
Theta = 1
mu = 43
z = 1
upsilon = 0

Phi_ABR_0 = np.zeros(len(alpha))
Phi_ABR_1 = np.zeros(len(alpha))
Phi_ABR_2 = np.zeros(len(alpha))
Phi_ABR_3 = np.zeros(len(alpha))
Phi_ABR = np.zeros(len(alpha))
for i in range(len(alpha)):
    Phi_ABR_0[i] = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,0.25)
    Phi_ABR_1[i] = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,0.5)
    Phi_ABR_2[i] = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,1.5)
    Phi_ABR_3[i] = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,2)
    Phi_ABR[i] = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,0)
    print(i)
plt.plot(alpha,Phi_ABR, color = "Purple", label = "ABR")
plt.plot(alpha,Phi_ABR_0,color = "Black", label = "Kappa = 0.25")
plt.plot(alpha,Phi_ABR_1,color = "Red", label = "Kappa = 0.5")
plt.plot(alpha,Phi_ABR_2,color = "Darkblue", label = "Kappa = 1.5")
plt.plot(alpha,Phi_ABR_3,color = "Darkorange", label = "Kappa = 2")
plt.title(f"Theta = {Theta}, mu = {mu}, z = {z}, upsilon = {upsilon}")
plt.xlabel("alpha")
plt.ylabel("Normalised potential")
plt.grid()
plt.xscale("log")
plt.legend()
plt.show()

