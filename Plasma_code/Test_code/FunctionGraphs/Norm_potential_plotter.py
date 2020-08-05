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

#Normalised inputs
Theta = np.linspace(0,1000,1001).tolist()
mu = 43
z = 1
alpha = 1
upsilon = 0

input_list = [Theta,mu,z,alpha,upsilon] #len = 4
title_list = ["theta","mu","relative ion charge","normalised dust radius","normalised plasma flow speed"]
label_list = [r"$\Theta$",r"$\mu$",r"$\z$",r"$\alpha$",r"$\upsilon$"]

for i in range(len(input_list)):
    if type(input_list[i]) == list:
        variable_index = i
for i in range(len(input_list)):
    if i != variable_index:
        input_list[i] = input_list[i]*np.ones(len(input_list[variable_index]))
    else:
        input_list[i] = np.array(input_list[i])

def grapher(input_list,variable_index):

    Phi_list = []
    for i in range(len(input_list[variable_index])):
        Phi = get_Norm_Potential(input_list[0][i],input_list[1][i],input_list[2][i],input_list[3][i],input_list[4][i])
        Phi_list.append(Phi)
    Phi = np.array(Phi_list)
    plt.plot(input_list[variable_index], Phi)
    plt.title(f"Variation of normalised dust surface potential with {title_list[variable_index]}")
    plt.ylabel("Normalised surface potential, " + r"$\Phi$")
    plt.xlabel(f"{label_list[variable_index]}")
    plt.grid()
    plt.show()

grapher(input_list,variable_index)
