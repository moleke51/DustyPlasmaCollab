import scipy as sp 
import numpy as np
import matplotlib.pyplot as plt
import HotABR as abr
import scipy.optimize as spl

def function(alpha,kappa):
    Theta = 0.01
    mu = 43
    z = 1
    upsilon = 0
    Phi = abr.potential_finder(Theta,mu,z,alpha,upsilon,kappa)
    return Phi

alpha,Phi = np.loadtxt("/Users/doganakpinar/Documents/Physics_Research/DustyPlasmaCollab/Plasma_code/Models/OMData1.txt",skiprows=1,unpack=True)

Phi_fit, Phi_cov = spl.curve_fit(function,alpha,Phi)
print(Phi_fit)
#plt.plot(alpha,Phi_fit(alpha))
#plt.grid()
#plt.show()