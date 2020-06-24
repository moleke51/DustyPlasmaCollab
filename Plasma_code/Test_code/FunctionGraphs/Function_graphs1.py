#============================IMPORT STUFF============================

import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
import numpy as np
from scipy.special import lambertw


#============================FUNCTIONS============================

#This function gives the real value of the Lambert W0 function provided x > -1/e
def realLambertW(x):
    if type(x) is float or type(x) is int :
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

#limiting values for potential, general form for OML and TS - eqn 4.2 in Willis' thesis
def Limit_potential_finder(nu,mu,C,Theta,m_a):

    x = nu*sp.log(m_a) + mu*sp.log(Theta) + C

    return x

#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_OML,phi_TS,alpha_OML,alpha_TS,alpha):

    x = ((phi_TS - phi_OML)/(sp.log(alpha_TS) - sp.log(alpha_OML)))*sp.log((alpha)/(alpha_TS)) + phi_TS

    return x 

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
def OML_surface_potential_finder(Theta,mu,Z): #we are ingnoring absorption radii outside the dust grain

    x = sp.absolute((Theta/Z) - realLambertW((mu*sp.sqrt(Theta)/Z)*sp.exp(Theta/Z)))

    return x

def ABR_surface_potential_graph(phi):

    x = (4*(phi)**(3/2)*((2*phi)-3)*((2*phi)+1))/((2*phi)-1)**3

    return x

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
def MOML_surface_potential_finder(Theta,mu,Z,gamma): #we are ingnoring absorption radii outside the dust grain

    x = sp.absolute(Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((Z**2)*(mu)**2)))

    return x


#This allows the chemical symbols or element names to be entered in either upper or lower case.
def speciesinput():
    word = input("Enter the plasma ion species; ")
    if len(word)==1 : 
        word = word.upper()
    elif len(word)==2 :
        word = word.lower().capitalize()
    elif word.lower() == 'override':
        word = 'override'
    else:
        word = word.lower()
    return(word)
#============================MAIN CODE=============================

#CONSTANTS

u = 1.66e-27 #kg
m_e = 9.11e-31 #kg
k = 1.38e-23 #m2kgs-2K-1
e = 1.60e-19 #C
epsilon_0 = 8.85e-12 #Fm-1
v_threshold = 5 #m/s 

#ASSUMPTIONS
#SPHRERICAL SYMMETRY
#NO COLLISIONS
#NO MAGNETIC FIELD (B = 0)
#NO THERMEONIC EMISSION - NEGATIVE CHARGED DUST
#ION CURRENT + ELECTRON CURRENT = 0 (STEADY STATE)
#QUASI-NEUTRALITY


elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
elementList.append('override')


#ASK USER FOR INPUT SPECICES 


species = speciesinput()
while (species in elementList) == False:
    print('This species does not exist')
    species = speciesinput()
        
    
else:
    if species == 'override':
        mu = float(input('Enter the mu value; ')) #Mu value

    else:
        proton_number = 'pt.' + species + '.number'
        mass = 'pt.' + species + '.mass'
        m_a = eval(mass)
        m_i = (m_a)*u #kg
        print(species, m_i, ' kg')
        mu = sp.sqrt(m_i/m_e) #Mu value

    mu = sp.sqrt(m_i/m_e)


#ASK FOR CHARGE ON ION
if species.lower()=='override':
    Z = float(input('Enter the charge on the plasma ions; '))
else:
    Z_max = eval(proton_number)
    Z = float(input("Enter the charge on the plasma ions; ")) # Z is the ion charge

    while (Z > Z_max or Z <= 0) == True:
        if (Z > Z_max) == True:
            print('The maximum charge for', species, 'is', Z_max, '+')
            Z = float(input("Enter the charge on the plasma ions; "))
        elif (Z < 0) == True:
            print('The charge value must be greater than 0')
            Z = float(input("Enter the charge on the plasma ions; "))
        else:
            print('Species must be an ion')
            Z = float(input("Enter the charge on the plasma ions; "))
    else:
        print(Z,'e')


#STATIC OR FLOWING

v_flow = sp.absolute(float(input('Set plasma flow speed (m/s); ')))

if v_flow == 0:
    print('This is a static plasma')
elif v_flow <= v_threshold:
    print('We will assume a static plasma for this model as the speed is low')
else:
    print('This is a flowing plasma')

print(v_flow, ' m/s')


#ASK ION TEMPERATURE
T_i = float(input('Set ion temperature (Kelvin); '))
while (T_i <= 0) == True:
    print('Violates 3rd law of Thermodynamics')
    T_i = float(input('Set ion temperature (Kelvin); '))
else:
    print(T_i)


#ASK ELECTRON TEMPERATURE

T_e = float(input('Set electron temperature (Kelvin); '))
while (T_e <= 0) == True:
    print('Violates 3rd law of Thermodynamics')
    T_e = float(input('Set ion temperature (Kelvin); '))
else:
    print(T_e)

#THETA VALUE

#Theta = T_i/T_e
#print(Theta)

Theta = float(input('Enter a theta value, where theta is T_i/T_e'))
while (Theta < 0) == True:
    print('Violates 3rd law of Thermodynamics, we cannot have a negative temperature ratio')
    Theta = float(input('Enter a theta value, where theta is T_i/T_e'))
else:
    print(Theta)

#BOUNDARY CONDITIONS

alpha_OML = 1.25*(Theta)**0.4
alpha_TS = 50
#OML used in this range
alpha_1 = sp.logspace(0,sp.log10(alpha_OML),10)
Theta_OML = sp.array([Theta] * len(alpha_1))
#Linear fit model used in this range
alpha_2 = sp.logspace(sp.log10(alpha_OML), sp.log10(50),10)
#MOML used in this range
alpha_3 = sp.logspace(sp.log10(50),2,10)
Theta_MOML = sp.array([Theta] * len(alpha_3))

#LINEAR FIT

#if Theta <= 2:
    #phi_OML = Limit_potential_finder(0.405,0.253 + 0.021*sp.log(m_a),2.454,Theta,m_a)
    #phi_TS = Limit_potential_finder(0.456,0,3.179,Theta,m_a)
    #print(phi_OML)
    #print(phi_TS)
#else:
    #phi_OML = Limit_potential_finder(0.401,-0.122 + 0.029*sp.log(m_a),2.698,Theta,m_a)
    #phi_TS = Limit_potential_finder(0.557,-0.386 - 0.024*sp.log(m_a),3.399,Theta,m_a)
    #print(phi_OML)
    #print(phi_TS)

Phi_MOML = MOML_surface_potential_finder(Theta_MOML,mu,Z,5/3)
phi_MOML = Phi_MOML[1]
Phi_OML = OML_surface_potential_finder(Theta_OML,mu,Z)
phi_OML = Phi_OML[1]
Phi_func_of_alph = Linear_function(phi_OML,phi_MOML,alpha_OML,alpha_TS,alpha_2)


plt.plot(alpha_1,Phi_OML, label='OML model', color = 'Red' , )
plt.plot(alpha_2, Phi_func_of_alph , label='Linear fit model' , color = 'Green' )
plt.plot(alpha_3, Phi_MOML, label = 'MOML model', color = 'Blue')
plt.title('Dust surface potential variation with normalised distance')
plt.ylabel('Normalised surface potential')
plt.xlabel('Normalsied distance')
plt.xscale('log')
plt.grid()
plt.legend()
plt.show()