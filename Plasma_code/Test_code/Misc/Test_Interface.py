def eval_input(x):
    x = x.replace('^','**')
    if '**' in x:
        x = x.split('**')
        a = x[0].split('*')
        b = x[1].split('*')
        A = 1
        for i in range(len(a)-1):
            A *= float(a[i])
        for i in range(1,len(b)):
            A *= float(b[i])
        B = float(a[-1])**float(b[0])
        X = A*B
    else:
        x = x.split('*')
        X = 1
        for i in range(len(x)):
            X *= float(x[i])
    return X

import numpy as np
import periodictable as pt

e = 1.60e-19 #[C] The charge on an electron
epsilon_0 = 8.85e-12 #[F][m^-1] The permitivity of free space
k_B = 1.38e-23 #[m^2][kg][s^-2][K^-1] The Boltzmann constant
m_e = 9.11e-31 #[kg] The mass of an electron
u = 1.66e-27 #[kg] The mass of a nucleon 

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

#Define the periodic table
elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
elementList.append('override')
override = 'override'
print(elementList)
print(pt.deuterium.mass)
while True:
    species = speciesinput()
    while (species in elementList) == False:
        print('This species does not exist')
        species = speciesinput()
    i = np.floor(elementList.index(species)/2)
    proton_number = elements[i].number
    m_i = (elements[i].mass) #[kg]
    mu = np.sqrt(m_i/m_e) #Mu value
    print(m_i,proton_number)
    