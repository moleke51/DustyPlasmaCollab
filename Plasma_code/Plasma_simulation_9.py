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
        x = x.split('**',)
        a = x[0].split('*')
        b = x[-1].split('*')
        A = 1
        for i in range(len(a)-1):
            A *= float(a[i])
        for i in range(1,len(b)):
            A *= float(b[i])
        if len(x) == 3:
            c,d = x[1].split('*')
            B = (float(a[-1])**float(c))*(float(d)**float(b[0])) 
        else:
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
            return(np.linspace(minimum,maximum,100000).tolist(),var_counter+1)
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
         
def get_Norm_Potential(Theta,mu,z,alpha,upsilon,variable_counter,previous_model = None,previous_phi = None):
    modellist = mdl.modelpicker('Plasma_code/Models/',Theta,mu,z,alpha,upsilon)
    priority = 0
    for model in modellist:
        __import__(model.get_name())
        if model.priority() > priority:
            
            priority = model.priority()
            modelindex = modellist.index(model)
    
    if variable_counter == 0:
        print(modellist[modelindex])
        return modellist[modelindex].potential_finder()
    else:
        if modellist[modelindex].get_name() == 'ABR' and previous_model == 'ABR':
            return(previous_phi, previous_model)
        else:
            return(modellist[modelindex].potential_finder(), modellist[modelindex].get_name(),modellist[modelindex].get_colour())

def grapher(input_list,variable_index,Dimensionless,variable_counter=1):
    Phi_list = []
    mod_name = None
    mod_list =[]
    Phi = None
    mod_range = []
    colour_list = []
    if Dimensionless == False:
        #Physical
        title_list = ["ion temperature","electron temperature","electron density","relative ion charge","ion mass","dust radius","plasma flow speed"]
        label_list = [r"$T_{i} \ [$K$]$",r"$T_{e} \ [$K$]$",r"$n_{0} \ [$m^{-3}$]$",r"$z$",r"$m_{i} \ [$u$]$",r"$a \ [$m$]$",r"$v \ [$$ms^{-1}$$]$"]
        Theta = np.array(input_list[0]/input_list[1])
        mu = np.array(np.sqrt(input_list[4]/m_e))
        z = input_list[3]
        alpha = np.array(input_list[5]/np.sqrt((epsilon_0*k_B*input_list[1])/(input_list[2]*e**2)))
        upsilon = np.array(input_list[6]/np.sqrt(2*k_B*input_list[0]/input_list[4]))
        for i in range(len(input_list[variable_index])):
            Phi, mod_name, colour = get_Norm_Potential(Theta[i],mu[i],z[i],alpha[i],upsilon[i],variable_counter,mod_name,Phi)
            Phi_list.append(Phi)
            if mod_name not in mod_list:
                mod_list.append(mod_name)
                mod_range.append(i)
                colour_list.append(colour)
                print(mod_name)
        mod_range.append(len(input_list[variable_index])+1)     

        Phi = np.array(Phi_list)
        #un-normalise phi
        Phi = (Phi*k_B*input_list[1])/e
        for i in range(len(mod_list)):
            plt.plot(input_list[variable_index][mod_range[i]:mod_range[i+1]-1], Phi[mod_range[i]:mod_range[i+1]-1],color = colour_list[i],label = mod_list[i])
        plt.title(f"Variation of dust surface potential with {title_list[variable_index]}")
        plt.ylabel("Surface potential, " + r"$\phi$ [V]")
        plt.xlabel(f"{label_list[variable_index]}")
        plt.legend()
        plt.grid()
        plt.show()

    else:
        #Dimensionless
        title_list = ["theta","mu","relative ion charge","normalised dust radius","normalised plasma flow speed"]
        label_list = [r"$\Theta$",r"$\mu$",r"$\z$",r"$\alpha$",r"$\upsilon$"]
        for i in range(len(input_list[variable_index])):
            Phi, mod_name, colour = get_Norm_Potential(input_list[0][i],input_list[1][i],input_list[2][i],input_list[3][i],input_list[4][i],variable_counter,mod_name,Phi)
            Phi_list.append(Phi)
            if mod_name not in mod_list:
                mod_list.append(mod_name)
                mod_range.append(i)
                print(mod_name)
                colour_list.append(colour)
        mod_range.append(len(input_list[variable_index])+1)     
        Phi = np.array(Phi_list)
        for i in range(len(mod_list)):
            plt.plot(input_list[variable_index][mod_range[i]:mod_range[i+1]-1], Phi[mod_range[i]:mod_range[i+1]-1],color = colour_list[i],label = mod_list[i])
        plt.title(f"Variation of normalised dust surface potential with {title_list[variable_index]}")
        plt.ylabel("Normalised surface potential, " + r"$\Phi$")
        plt.xlabel(f"{label_list[variable_index]}")
        plt.legend()
        plt.xscale('log')
        plt.grid()
        plt.show()

#Define the periodic table
elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
#==============================MAIN CODE================================#

#Ask the user to input the species of the ions in the plasma or to override and imput a custom value of mu.
variable_counter = 0
choice = None 
while choice != 'y' and choice != 'n':
     choice = input('Do you want to use dimensionless variables (y/n); ').lower()
else: 
    if choice == 'y':
        choice = True
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
            alpha, variable_counter = is_valid('alpha',['is_num','>=0'],variable_counter)

        upsilon = None
        while upsilon == None:
            upsilon, variable_counter = is_valid('upsilon',['is_num'],variable_counter)
        
        
    else:
        choice = False
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
            a,variable_counter = is_valid('dust radius',['is_num','>=0'],variable_counter,'meters')
        
        v = None
        while v == None:
            v,variable_counter = is_valid('plasma flow speed',['is_num','>=0'],variable_counter,'meters per second')

if variable_counter == 0:
    if choice == False:
        Theta = T_i/T_e
        lambda_D = sp.sqrt(epsilon_0*k_B*T_e/(n_e*e**2))
        alpha = a/lambda_D
        if T_i != 0:
            upsilon = v / sp.sqrt(2*k_B*T_i/m_i)
        else:
            upsilon = 0


    Phi = get_Norm_Potential(Theta,mu,z,alpha,upsilon,variable_counter)
    #Return the normalised potential
    print('The normalised dust grain surface potential is:',Phi)

    #Return the potentail and the charge if available.
    if choice == False:
        phi = (Phi * k_B * T_e)/(e)
        Q = DH_potential_to_charge(a,phi,lambda_D)
        Z = Q/e
        print('The dust grain surface potential is ' +str(phi)+ ' V')
        print('The charge on the dust grain is ' +str(Q)+ ' C')
        print('The relative charge of the dust grain is ' + str(Z))
else:
    print("Producing your graph, please wait.")
    if choice == False:
        input_list = [T_i,T_e,n_e,z,m_i,a,v] #len = 6
    elif choice == True:
        input_list = [Theta,mu,z,alpha,upsilon] #len = 4
    for i in range(len(input_list)):
        if type(input_list[i]) == list:
            variable_index = i
    for i in range(len(input_list)):
        if i != variable_index:
            input_list[i] = input_list[i]*np.ones(len(input_list[variable_index]))
        else:
            input_list[i] = np.array(input_list[i])
    
    grapher(input_list,variable_index,choice)
