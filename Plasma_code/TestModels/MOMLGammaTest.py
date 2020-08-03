import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve,bisect
import scipy.special as sps

def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64:
        w = sps.lambertw(x)
        if np.imag(w)==0:
            W = np.real(w)
            return(W)
        else:
            return('This value is outside the accepted range')
    elif type(x) is np.ndarray or type(x) is list:
        W = np.zeros(len(x))
        for i in range(0,len(x)):
            
            if np.imag(sps.lambertw(x[i]))==0:
                W[i]= np.real(sps.lambertw(x[i]))
            else:
                print('The value at position ' +str(i)+' is outside the accepted range')
        return(W)
    else:
        return('This is an invalid input')

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
#Define MOML equation to solve 
def MOML_function(Phi,Theta,mu,z,gamma): #gamma = 5/3 for static plasmas
    return (np.sqrt(Theta)/mu)*(1 - (1/Theta)*(Phi - 0.5*(np.log(2*np.pi*(1+gamma*Theta))-np.log(mu**2)))) - np.exp(Phi)

def planar_presheath(Theta,mu,z,gamma):
    return 0.5*np.log((2*np.pi*(1+gamma*Theta)/(mu**2))) - 0.5

Theta = np.logspace(-3,3,101)
#Theta = np.logspace(-1,0,101)
mu = 43
z = 1
Phi_1 = np.zeros(len(Theta)) #gamma = 1
Phi_2 = np.zeros(len(Theta)) #gamma = 5/3
Phi_3 = np.zeros(len(Theta)) #gamma = 3
Phi_pp1 = np.zeros(len(Theta)) #gamma = 1
Phi_pp2 = np.zeros(len(Theta)) #gamma = 5/3
Phi_pp3 = np.zeros(len(Theta)) #gamma = 3

for i in range(len(Theta)):
    Phi_1[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,1))
    Phi_2[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,5/3))
    Phi_3[i] = bisect(MOML_function,-10,10,args = (Theta[i],mu,z,3))

Phi_pp1 = 1*planar_presheath(Theta,mu,z,1)
Phi_pp2 = 1*planar_presheath(Theta,mu,z,5/3)
Phi_pp3 = 1*planar_presheath(Theta,mu,z,3)

#Cold ion wall limit
Phi_cw = planar_presheath(0,mu,z,1)*np.ones(len(Theta))

#Finding minimum of MOML for gamma = 1
min_point_1 = np.argmin(Phi_1)
print(f"For MOML prediction: Gamma = 1, Theta_min = {Theta[min_point_1]} and Phi_min = {Phi_1[min_point_1]}")

#Find point that planar pre-sheath diverges from the cold ion value
index_list_1 = []
for i in range(len(Phi_1)):
    if (np.absolute(Phi_pp1[i] - Phi_cw[i]) < 0.01):
        index_list_1.append(i)

max_point_1 = np.max(index_list_1)
print(f"For PP prediction: Gamma = 1, Theta_min = {Theta[max_point_1]} and Phi_min = {Phi_pp1[max_point_1]}")


'''
#Finding minimum of MOML for gamma = 5/3
min_index_2 = np.argmin(Phi_2)
print(f"Gamma = 5/3, Theta_min = {Theta[min_index_2]} and Phi_min = {Phi_1[min_index_2]}")
Theta_min_2 = Theta[min_index_2]
Phi_min_2 = Phi_2[min_index_2]

#Finding minimum of MOML for gamma = 3
min_index_3 = np.argmin(Phi_3)
print(f"Gamma = 3, Theta_min = {Theta[min_index_3]} and Phi_min = {Phi_1[min_index_3]}")
Theta_min_3 = Theta[min_index_3]
Phi_min_3 = Phi_3[min_index_3]

Theta_min = np.array([Theta_min_1,Theta_min_2,Theta_min_3])
Phi_min = np.array([Phi_min_1,Phi_min_2,Phi_min_3])
'''

plt.plot(Theta,Phi_1,color = "Red", label = "Gamma = 1")
#plt.plot(Theta,Phi_2,color = "Blue", label = "Gamma = 5/3")
#plt.plot(Theta,Phi_3,color = "Green", label = "Gamma = 3")
plt.plot(Theta,Phi_pp1,':',color = "DarkRed", label = "PP, gamma = 1")
#plt.plot(Theta,Phi_pp2,'-.',color = "DarkBlue", label = "PP, gamma = 5/3")
#plt.plot(Theta,Phi_pp3,'--',color = "DarkGreen", label = "PP, gamma = 3")
plt.plot(Theta[min_point_1],Phi_1[min_point_1],'o',color = "Black", label = "Minimum of MOML")
plt.plot(Theta[max_point_1],Phi_pp1[max_point_1],'o',color = "Black", label = "Max chosen value for PP")
plt.plot(Theta,Phi_cw,'--',color = "DarkBlue", label = "Cold ion wall limit")
plt.xlabel("Theta")
plt.ylabel("Normalised dust potential")
plt.legend()
plt.grid()
plt.xscale("log")
plt.show()