import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.optimize import fsolve
import Model as mdl

################################## Constants ##################################
#We define some of the universal constants
e = 1.60e-19 #[C] The charge on an electron
epsilon_0 = 8.85e-12 #[F][m^-1] The permitivity of free space
k_B = 1.38e-23 #[m^2][kg][s^-2][K^-1] The Boltzmann constant
m_e = 9.11e-31 #[kg] The mass of an electron
u = 1.66e-27 #[kg] The mass of a nucleon 

def get_Norm_Potential(Theta,mu,z,alpha,upsilon):
    modellist = mdl.modelpicker('Plasma_code/Models/',Theta,mu,z,alpha,upsilon)
    priority = 0
    for model in modellist:
        __import__(model.get_name())
        if model.priority() > priority:
            priority = model.priority()
            modelindex = modellist.index(model)

    return (modellist[modelindex].potential_finder())

#Physical inputs
T_i = np.linspace(1,1000,1001).tolist() #Ion temperature
T_e = 1 #Electron temperature
m_i = 1.67e-27 #Mass of ion --> Hydrogen in this case
n_0 = 1e10 #Electron density at infinity
a = 1 #Dust grain radius
v = 1 #velocity
#lambda_d = np.sqrt((epsilon_0*k_B*T_e)/(n_0*e**2))
z = 1

#Theta = np.linspace(0,1000,1001).tolist()
#mu = 43
#z = 1
#alpha = 1
#upsilon = 0

input_list = [T_i,T_e,n_0,z,m_i,a,v] #len = 6
title_list = ["ion temperature","electron temperature","electron density","relative ion charge","ion mass","dust radius","plasma flow speed"]
label_list = [r"$T_{i} \ ($K$)$",r"$T_{e} \ ($K$)$",r"$n_{0} / ($m^{-3}$)$",r"$z$",r"$m_{i} / ($u$)$",r"$a / ($m$)$",r"$v / ($ms^{-1}$)$"]

for i in range(len(input_list)):
    if type(input_list[i]) == list:
        variable_index = i
for i in range(len(input_list)):
    if i != variable_index:
        input_list[i] = input_list[i]*np.ones(len(input_list[variable_index]))
    else:
        input_list[i] = np.array(input_list[i])

def grapher(input_list,variable_index):

    Theta = np.array(input_list[0]/input_list[1])
    mu = np.array(np.sqrt(input_list[4]/m_e))
    z = input_list[3]
    alpha = np.array(input_list[5]/np.sqrt((epsilon_0*k_B*input_list[1])/(input_list[2]*e**2)))
    upsilon = np.array(input_list[6]/np.sqrt(2*k_B*input_list[0]/input_list[4]))

    Phi_list = []
    for i in range(len(input_list[variable_index])):
        Phi = get_Norm_Potential(Theta[i],mu[i],z[i],alpha[i],upsilon[i])
        Phi_list.append(Phi)
    Phi = np.array(Phi_list)
    #un-normalise phi
    Phi = (Phi*k_B*input_list[1])/e
    plt.plot(input_list[variable_index], Phi)
    plt.title(f"Variation of dust surface potential with {title_list[variable_index]}")
    plt.ylabel("Surface potential, " + r"$\phi$ (V)")
    plt.xlabel(f"{label_list[variable_index]}")
    plt.grid()
    plt.show()

grapher(input_list,variable_index)