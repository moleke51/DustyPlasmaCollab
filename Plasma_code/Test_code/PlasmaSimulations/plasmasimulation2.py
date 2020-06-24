#============================IMPORT STUFF============================


import scipy as sp
import matplotlib.pyplot as plt
import periodictable as pt
import numpy as np
from scipy.special import lambertw


#============================FUNCTIONS============================

def OML_surface_potential_finder(Theta,mu,Z): #we are ingnoring absorption radii outside the dust grain

    x = (Theta/Z) - np.absolute(lambertw((mu*sp.sqrt(Theta)/Z)*sp.exp(Theta/Z)))

    return x



#============================MAIN CODE============================

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


#ASK USER FOR INPUT SPECICES 

species = input("Enter the ion species; ")
#This allows the chemical symbols or element names to be entered in either upper or lower case.
if len(species)==1 : 
    species = species.upper()
elif len(species)==2 :
    species = species.lower().capitalize()
    print(species)

else:
    species = species.lower()
    print(species)
elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())

while (species in elementList) == False:
    print('This species does not exist')
    species = input("Enter the ion species; ")

else:
    proton_number = 'pt.' + species + '.number'
    mass = 'pt.' + species + '.mass'
    m_a = eval(mass)
    m_i = (m_a)*u #kg
    print(species, m_i, ' kg')
#MU VALUE

mu = sp.sqrt(m_i/m_e)
#print(mu)

#ASK FOR CHARGE ON ION

Z_max = eval(proton_number)
Z = float(input("Set charge on ion; ")) # Z is the relative ion charge

while (Z > Z_max or Z <= 0) == True:
    if (Z > Z_max) == True:
        print('The maximum charge for', species, 'is', Z_max, '+')
        Z = float(input("Set charge on ion; "))
    elif (Z < 0) == True:
        print('The charge value must be greater than 0')
        Z = float(input("Set charge on ion; "))
    else:
        print('Species must be an ion')
        Z = float(input("Set charge on ion; "))
else:
    print(Z,'C')


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

Theta = T_i/T_e
#print(Theta)



