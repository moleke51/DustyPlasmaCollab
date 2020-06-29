#==============================ASSUMPTIONS==============================#
#Spherical symmetry
#No collisions
#No magnetic field (B = 0)
#No thermionic emission ==> negatively charged dust
#Steady state ==> Ion current + electron current = 0 at the dust particle 
#Quasi-neutrality
#============================DEFINE VARIABLES===========================#
# e: the magnitude of the charge on an electron
# epsilon_0: permittivity of free space
# k_B: Boltzmann's constant
# u: atomic mass unit
# a: dust radius 
# lambda_D = [(epsilon_0 * k_B * T_e)/(n_e * e^2)]^1/2: Debye length
# alpha = [a/lambda_D]: normalized dust radius
# rho = [r/lambda_D]: normalised distance from dust centre
# m_a: atom mass
# m_i: ion mass
# m_e: electron mass
# mu = [m_i/m_e]^1/2
# n_e: electron density
# T_i: ion temperature
# T_e: electron temperature
# Theta = T_i/T_e 
# gamma: adiabatic constant
# v: flow speed
# upsilon = v / [2*(k_B * T_i)/m_i]^1/2: normalised flow speed
# z: relative ion charge
#============================IMPORT PACKAGES============================#
import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.optimize import fsolve
import Model as mdl
#==============================KEY WORDS================================#
masteroverride = 'normalise' 
override = 'override'
#==============================CONSTANTS================================#
#We define some of the universal constants
e = 1.60e-19 #[C] The charge on an electron
epsilon_0 = 8.85e-12 #[F][m^-1] The permitivity of free space
k_B = 1.38e-23 #[m^2][kg][s^-2][K^-1] The Boltzmann constant
m_e = 9.11e-31 #[kg] The mass of an electron
u = 1.66e-27 #[kg] The mass of a nucleon 
#==============================FUNCTIONS===============================#
#Define the Debye-Huckle potential.
def DH_potential_to_charge(dustradius,Phi_a,lambda_d):
    x = 4*sp.pi*(epsilon_0)*dustradius*Phi_a*sp.exp((dustradius)/(lambda_d))
    return x

#It should be noted that the Debye-Huckle potential reduces into the point charge potential when a<<lambda_d
def Spherical_potential_to_charge(dustradius,Phi_a):
    x = 4*sp.pi*(epsilon_0)*dustradius*Phi_a
    return x

#This allows the chemical symbols or element names to be entered in either upper or lower case.
def speciesinput():
    word = input("Enter the plasma ion species; ")
    if len(word)==1 : 
        word = word.upper()
    elif len(word)==2 :
        word = word.lower().capitalize()
    else:
        word = word.lower()
    return(word)
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


def is_valid(name,requirements,var_counter,units = None):
    units_message = '; '
    if units != None:
        units_message = f' with units of {units}; '
    user_input = input(f'Enter the {name} value{units_message}')
    prefixes = {'Y': '*10**24', 'Z': '*10**21', 'E': '*10**18', 'P': '*10**15', 'T': '*10**12', 'G': '*10**9', 'M': '*10**6', 'k': '*10**3', 'm': '*10**-3', 'u': '*10**-6', 'n': '*10**-9', 'p': '*10**-12', 'f': '*10**-15', 'a': '*10**-18', 'z': '*10**-21', 'y': '*10**-24'}
    check = True
    message = ''
    if user_input.lower() == 'variable':
        if user_input.lower() == 'variable' and var_counter == 0:
            maximum = None
            minimum = None
            while maximum == None:
                maximum, a = is_valid('maximum '+name,requirements,1,units)
            requirements.append(f'<{maximum}')
            while minimum == None:
                minimum, a = is_valid('minimum '+name,requirements,1,units)
            return(np.linspace(minimum,maximum,1000).tolist(),var_counter+1)
        elif var_counter == 1: 
            print('The variable has already been selected')
        else:
            print('This parameter cannot be varied')
    if 'is_num' in requirements:
        try:
            str_number = eval_input((user_input.replace('eV','ev').replace('ev','*11594')).translate(str.maketrans(prefixes)))
            number = float(str_number)
            if 'is_int' in requirements:
                if '.' in user_input:
                    check = False
                    if message == '':
                        message = f'{name} must be an integer'
                    else:
                        message += f', {name} must be an integer'
            if '>0' in requirements:
                if number <= 0:
                    check = False
                    if message == '':
                        message = f'{name} must be greater than zero'
                    else:
                        message += f', {name} must be greater than zero'
            if '>=0' in requirements: 
                if number < 0:
                    check = False
                    if message == '':
                        message = f'{name} must be greater than or equal to zero'
                    else:
                        message += f', {name} must be greater than or equal to zero' 
            for req in requirements:
                if '<=' in req:
                    max_num = float(req.strip('<='))
                    if number > max_num:
                        check = False
                        if message == '':
                            message = f'{name} must be smaller than {max_num}'
                        else:
                            message += f', {name} must be smaller than {max_num}'
            for req in requirements:
                if '<' in req and '<=' not in req:
                    max_num = float(req.strip('<'))
                    if number >= max_num:
                        check = False
                        if message == '':
                            message = f'{name} must be smaller than {max_num}'
                        else:
                            message += f', {name} must be smaller than {max_num}'
        except ValueError:
            check = False
            message = (f'{name} must be a number.')
            
        except NameError:
            check = False
            message = (f'{name} must be a number.')
            
        except TypeError:
            check = False
            message = (f'{name} must be a number.')
            
        if check == True:
            return(str_number, var_counter)
        else:
            print(message)
            return(None,var_counter)
    else:
        return(user_input,var_counter)
         
def get_Norm_Potential(Theta,mu,z,alpha,upsilon):
    modellist = mdl.modelpicker('Plasma_code/Models/',Theta,mu,z,alpha,upsilon)
    priority = 0
    for model in modellist:
        __import__(model.get_name())
        if model.priority() > priority:
            
            priority = model.priority()
            modelindex = modellist.index(model)

    print(modellist[modelindex])
    return modellist[modelindex].potential_finder()
#Define the periodic table
elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
#==============================MAIN CODE================================#

#Ask the user to input the species of the ions in the plasma or to override and imput a custom value of mu.
variable_counter = 0
choice = input('Do you want to use dimensionless variables (y/n); ').lower()
while choice != 'y' and choice != 'n':
     choice = input('Do you want to use dimensionless variables (y/n)').lower()
else: 
    if choice == 'y':
        mu = None
        while mu == None:
            mu,variable_counter = is_valid('mu',['is_num','>0'],variable_counter)

        z = None
        while z == None:
            z,variable_counter = is_valid('relative ion charge',['is_num','>0','is_int'],variable_counter)
        
        Theta = None
        while Theta == None:
            Theta, variable_counter = is_valid('Theta',['is_num','>=0'],variable_counter)
        
        alpha = None
        while alpha == None:
            alpha, variable_counter = is_valid('alpha',['is_num','>0'],variable_counter)

        upsilon = None
        while upsilon == None:
            upsilon, variable_counter = is_valid('upsilon',['is_num'],variable_counter)
        
        
    else:
        species = speciesinput()
        while (species in elementList) == False:
            print('This species does not exist')
            species = speciesinput()
    
        else:
            #Calculate the mu value.
            index = np.floor(elementList.index(species)/2)
            if index == 119: #Deuterium
                index = 1
                m_a = pt.deuterium.mass
            elif index == 120: #Tritium
                index = 1
                m_a = pt.tritium.mass
            else:   
                m_a = (elements[index].mass) #[kg]
            m_i = (m_a)*u #[kg]
            mu = np.sqrt(m_i/m_e) #Mu value
            proton_number = elements[index].number
            z_max = proton_number
        
        
        
        z = None
        while z == None:
            z,usless = is_valid('relative ion charge',['is_num','>0','is_int',f'<={z_max}'],-1,'1.60*10^-19 C')
        
        T_i = None
        while T_i == None:
            T_i,variable_counter = is_valid('ion temperature',['is_num','>0'],variable_counter,'Kelvin')
        
        T_e = None
        while T_e == None:
            T_e,variable_counter = is_valid('electron temperature',['is_num','>0'],variable_counter,'Kelvin')   
        
        n_e = None
        while n_e == None:
            n_e,variable_counter = is_valid('electron density',['is_num','>0'],variable_counter,'electrons per meter cubed')

        a = None
        while a == None:
            a,variable_counter = is_valid('dust radius',['is_num','>0'],variable_counter,'meters')
        
        v = None
        while v == None:
            v,variable_counter = is_valid('plasma flow speed',['is_num','>=0'],variable_counter,'meters per second')

if variable_counter == 0:
    if choice == 'n':
        Theta = T_i/T_e
        lambda_D = sp.sqrt(epsilon_0*k_B*T_e/(n_e*e**2))
        alpha = a/lambda_D
        if T_i != 0:
            upsilon = v / sp.sqrt(2*k_B*T_i/m_i)
        else:
            upsilon = 0


    Phi = get_Norm_Potential(Theta,mu,z,alpha,upsilon)
    #Return the normalised potential
    print('The normalised dust grain surface potential is:',Phi)

    #Return the potentail and the charge if available.
    if choice == 'n':
        phi = (Phi * k_B * T_e)/(e)
        Q = DH_potential_to_charge(a,phi,lambda_D)
        print('The dust grain surface potential is ' +str(phi)+ 'V')
        print('The charge on the dust grain is ' +str(Q)+ 'C')
else:
    print('I need a graph pls')
    