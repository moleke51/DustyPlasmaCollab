import scipy as sp 
import numpy as np
import matplotlib.pyplot as plt
import HotABR as abr
import scipy.optimize as spo

def function(alpha,kappa):
    Theta = 1
    mu = 43
    z = 1
    upsilon = 0
    Phi_list = []
    for i in range(len(alpha)):
        Phi = abr.potential_finder(Theta,mu,z,alpha[i],upsilon,kappa)
        Phi_list.append(Phi)
    Phi = np.array(Phi_list)
    return Phi

<<<<<<< HEAD
alpha,Phi = np.loadtxt("/Users/georgedoran/Documents/'George Doran'/Imperial/Physics/'Dust in plasma'/DustyPlasmaCollab/Plasma_code/Models/OMData1.py",skiprows=1,unpack=True)
=======
>>>>>>> 978249a019bb7387dfc0b9a27a1132b2b3002534

alpha,Phi = np.loadtxt("/Users/doganakpinar/Documents/Physics_Research/DustyPlasmaCollab/Plasma_code/Models/OMData3.txt",skiprows=1,unpack=True)
Kappa, cov = spo.curve_fit(function,alpha,Phi,0.25)
print(Kappa,cov)
plt.plot(alpha,function(alpha,Kappa))
plt.grid()
plt.show()