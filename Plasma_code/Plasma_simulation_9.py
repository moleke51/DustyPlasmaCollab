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
    elif word.lower() == override:
        word = override
    elif word.lower() == masteroverride:
        word = masteroverride
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


def is_valid(user_input,name,requirements):
    check = True
    message = ''
    if 'is_num' in requirements:
        try:
            number = float(user_input)
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
        return user_input
    else:
        print(message)

         
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
elementList.append(override)
elementList.append(masteroverride)

#==============================MAIN CODE================================#

#Ask the user to input the species of the ions in the plasma or to override and imput a custom value of mu.
variable_counter = 0
choice = input('Do you want to use dimensionless variables (y/n)').lower()
while choice != 'y' and choice != 'n':
     choice = input('Do you want to use dimensionless variables (y/n)').lower()
else: 
    if choice == 'y':
        counter = 0
        while counter == 0:
            is_valid(mu)
            '''
            try:
                mu = input('Enter the mu value; ')
                if mu.lower() == 'variable':
                    variable_counter = 1
                    counter = 1
                else:
                    mu = eval_input(mu) #Mu value
                    if mu > 0:
                        counter = 1
                    else:
                        print('Mu must be greater than zero.')
            except ValueError:
                print('Invalid mu value.')
            except NameError:
                print('Invalid mu value.')
            except TypeError:
                print('Invalid mu value.')
            except SyntaxError:
                print('Invalid mu value.')
            '''
        while counter == 1:
            try:
                z = input('Enter the relative charge on the plasma ions, [1.60*10^-19 C] ; ') # z is the ion relative charge
                if '.' in z:
                    print('The relative charge must be an integer.')
                else:
                    z = float(z)
                if z>0:
                    counter=2
                else:
                    print('The charge must be greater than zero.')
            except ValueError:
                print('Invalid charge value.')
            except NameError:
                print('Invalid charge value.')
            except TypeError:
                print('Invalid charge value.')
            except SyntaxError:
                print('Invalid charge value.')

        while counter == 2:
            try:
                #THETA VALUE
                Theta = eval_input(input('Set the Theta value; '))
                if Theta >= 0:
                    counter = 3
                else:
                    print('Theta must be greater than or equal to zero.')
                    
            except ValueError:
                print('Invalid Theta value.')
            except NameError:
                print('Invalid Theta value.')
            except TypeError:
                print('Invalid Theta value.')
            except SyntaxError:
                print('Invalid Theta value.')
        while counter == 3:
                try:
                    alpha = eval_input(input('Set the alpha value; '))
                    if alpha > 0:
                        counter = 4
                    else:
                        print('Alpha must be greater than zero.')
                except ValueError:
                    print('Invalid alpha value.')
                except NameError:
                    print('Invalid alpha value.')
                except TypeError:
                    print('Invalid alpha value.')
                except SyntaxError:
                    print('Invalid alpha value.')
        while counter == 4:
            try:
                upsilon = sp.absolute(eval_input(input('Set the upsilon value; ')))
                counter = 5
                print(str(upsilon).replace('e','*10^'))
            except ValueError:
                print('Invalid flow speed.')
            except NameError:
                print('Invalid flow speed.')
            except TypeError:
                print('Invalid flow speed.')
            except SyntaxError:
                print('Invalid flow speed.')  

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
        counter = 0
        while counter==0:
            try:
                #Ask the user for the relative charge on the ions in the plasma which should be an integer.
                z = input('Enter the relative charge on the plasma ions, [1.60*10^-19 C] ; ') # z is the ion relative charge
                if '.' in z:
                    print('The relative charge must be an integer.')

                else:
                    z = float(z)
                if (z > z_max) == True:
                    print('The maximum relative charge for', species, 'is +' +str(z_max))
                elif (z < 0) == True:
                    print('The relative charge value must be greater than 0.')
                    
                elif (z == 0) == True:
                    print('Species must be an ion.')
                    
                else:
                    counter = 1
            except ValueError:
                print('Invalid charge value.')
            except NameError:
                print('Invalid charge value.')
            except TypeError:
                print('Invalid charge value.')
            except SyntaxError:
                print('Invalid charge value.')
        while counter == 1:
            try:
                #Ask the user for the plasmas ion temperature in Kelvin, for temperatures in electron volts input eV after the value.
                T_i = input('Set ion temperature [K]; ').replace('eV','ev').replace('ev','*11594')
                T_i = eval_input(T_i)
                if (T_i < 0) == True:
                    print('The temperature must be greater or equal to zero')
                else:
                    counter = 2
            except ValueError:
                print('Invalid ion temperature.')
            except NameError:
                print('Invalid ion temperature.')
            except TypeError:
                print('Invalid ion temperature.')
            except SyntaxError:
                print('Invalid ion temperature.')
        while counter == 2:
            try:
                #Ask the user for the electron temperature in Kelvin, for temperatures in electron volts input eV after the value.
                T_e = input('Set electron temperature [K]; ').replace('eV','ev').replace('ev','*11594')
                T_e = eval_input(T_e)
                if (T_e <= 0) == True:
                    print('The temperature must be greater than zero')
                else:
                    counter = 3
                    Theta = T_i/T_e      
            except ValueError:
                print('Invalid electron temperature.') 
            except NameError:
                print('Invalid electron temperature.')
            except TypeError:
                print('Invalid electron temperature.')
            except SyntaxError:
                print('Invalid electron temperature.')
        while counter == 3:
            try:
                #The user will then be asked for an electron density in electrons per metre cubed, and a dust grain radius in metres.
                n_e = eval_input(input('Set the electron density [electrons][m^-3] ').replace('k','*10**3').replace('M','*10**6').replace('G','*10**9').replace('T','*10**12').replace('P','*10**15'))
                if (n_e <= 0) == True:
                    print('The electron density must be greater than 0')   
                else:
                    counter = 4
                    lambda_D = sp.sqrt(epsilon_0*k_B*T_e/(n_e*e**2))
            except ValueError:
                print('Invalid electron density')
            except NameError:
                print('Invalid electron density')
            except TypeError:
                print('Invalid electron density')
            except SyntaxError:
                print('Invalid electron density')
        while counter == 4:
            try:
                a = eval_input((input('Set the dust radius [m] ').replace('m','*10**-3').replace('u','*10**-6').replace('n','*10**-9')))
                if (a <= 0) == True:
                    print('The dust radius must be greater than 0')
                else:
                    counter = 5
                    alpha = a/lambda_D
            except ValueError:
                print('Invalid dust radius density')
            except NameError:
                print('Invalid dust radius density')
            except TypeError:
                print('Invalid dust radius density')
            except SyntaxError:
                print('Invalid dust radius density')
#Ask the user for the flow speed of the plasma in metres per second, for stationary plasmas input 0.
        while counter == 5:
            try:
                v = input('Set plasma flow speed (m/s); ')
                v = sp.absolute(eval_input(v))
                counter = 1
                if T_i != 0:
                    upsilon = v / sp.sqrt(2*k_B*T_i/m_i)
                else:
                    upsilon = 0
            except ValueError:
                print('Invalid flow speed.')
            except NameError:
                print('Invalid flow speed.')
            except TypeError:
                print('Invalid flow speed.')
            except SyntaxError:
                print('Invalid flow speed.')     
    
                





#The user may directly input a theta value by typing override. 





        


Phi = get_Norm_Potential(Theta,mu,z,alpha,upsilon)
#Return the normalised potential
print('The normalised dust grain surface potential is:',Phi)

#Return the potentail and the charge if available.
if choice == 'n':
    phi = (Phi * k_B * T_e)/(e)
    Q = DH_potential_to_charge(a,phi,lambda_D)
    print('The dust grain surface potential is ' +str(phi)+ 'V')
    print('The charge on the dust grain is ' +str(Q)+ 'C')
