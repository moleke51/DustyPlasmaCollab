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
# upsilon = v * [(k_B * T_e)/m_i]^1/2: normalised flow speed
# Z: ion charge
# z = Z/e: normalised ion charge
#============================IMPORT PACKAGES============================#
import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.special import erf
#==============================CONSTANTS================================#
e = 1.60e-19 #[C]
epsilon_0 = 8.85e-12 #[F][m^-1]
k_B = 1.38e-23 #[m^2][kg][s^-2][K^-1]
m_e = 9.11e-31 #[kg]
u = 1.66e-27 #[kg]
Theta_critical = 0.05 #This is an arbitrary value
upsilon_critical = 0 #Until proved otherwise
#==============================FUNCTIONS===============================#

#This function gives the real value of the Lambert W0 function provided x > -1/e
def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64 :
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
def Limit_potential_finder(nu,mu,C,Theta,m_a): #nu and mu here are arbitrary constants

    x = nu*sp.log(m_a) + mu*sp.log(Theta) + C

    return x


#Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_OML,phi_TS,alpha_OML,alpha_TS,alpha):

    x = ((phi_TS - phi_OML)/(sp.log(alpha_TS) - sp.log(alpha_OML)))*sp.log((alpha)/(alpha_TS)) + phi_TS

    return x 


#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis (absolute values)
def OML_surface_potential_finder_abs(Theta,mu,z): #we are ingnoring absorption radii outside the dust grain

    x = sp.absolute((Theta/z) - realLambertW((mu*sp.sqrt(Theta)/z)*sp.exp(Theta/z)))

    return x

##SOML (Modified OML) model for normalised dust surface potential - eqn 2.138 in Thomas' thesis (absolute values)
def SOML_surface_potential_finder_abs(Theta,mu,z,upsilon):#u is the normalised plasma flow speed, normalised by sp.sqrt(2KT_i/m_i)

    s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
    s_2 = (sp.sqrt(sp.pi)*erf(upsilon))/(2*upsilon)

    x = sp.absolute((Theta*s_1)/(s_2) - realLambertW(((mu*z*sp.sqrt(Theta))/(s_2))*sp.exp((Theta*s_1)/(s_2))))
    
    return x

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis (absolute values)
def MOML_surface_potential_finder_abs(Theta,mu,z,gamma): #we are ingnoring absorption radii outside the dust grain

    x = sp.absolute(Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2)))

    return x

#SMOML (Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis (absolute values)
def SMOML_surface_potential_finder_abs(Theta,mu,z,gamma,upsilon):#u is the normalised plasma flow speed, normalised by sp.sqrt(2KT_i/m_i)

    s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
    s_2 = (sp.sqrt(sp.pi)*erf(upsilon))/(2*upsilon)

    x = sp.absolute((Theta*s_1)/(s_2) - realLambertW(((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))/(s_2))*sp.exp((Theta*s_1)/(s_2))) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2)))
    
    return x

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
def OML_surface_potential_finder(Theta,mu,z): #we are ingnoring absorption radii outside the dust grain

    x = (Theta/z) - realLambertW((mu*sp.sqrt(Theta)/z)*sp.exp(Theta/z))

    return x

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
def MOML_surface_potential_finder(Theta,mu,z,gamma): #we are ingnoring absorption radii outside the dust grain

    x = Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2))

    return x

##SOML (Modified OML) model for normalised dust surface potential - eqn 2.138 in Thomas' thesis
def SOML_surface_potential_finder(Theta,mu,z,upsilon):#u is the normalised plasma flow speed, normalised by sp.sqrt(2KT_i/m_i)

    s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
    s_2 = (sp.sqrt(sp.pi)*erf(upsilon))/(2*upsilon)

    x = (Theta*s_1)/(s_2) - realLambertW(((mu*z*sp.sqrt(Theta))/(s_2))*sp.exp((Theta*s_1)/(s_2)))
    
    return x

#SMOML (Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis
def SMOML_surface_potential_finder(Theta,mu,z,gamma,upsilon):#u is the normalised plasma flow speed, normalised by sp.sqrt(2KT_i/m_i)

    s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
    s_2 = (sp.sqrt(sp.pi)*erf(upsilon))/(2*upsilon)

    x = (Theta*s_1)/(s_2) - realLambertW(((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))/(s_2))*sp.exp((Theta*s_1)/(s_2))) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2))
    
    return x

#ABR model for normalised dust surface potential
def open_poly(chemical_symbol):
    fit = sp.loadtxt('Database/Polyfits/'+chemical_symbol+'polyfit.txt',dtype=float,unpack = True)
    fit = fit[::-1]
    x = sp.logspace(-4,4,10000)
    print(fit)
    poly = np.poly1d(fit)
    print(poly)
    y = poly(sp.log10(x))
    if len(species) <= 3:
        i = elementList.index(species)
        name = elementList[i+1].capitalize()
    plt.plot(x,y,color = 'purple',label = str(name))

def ABR_potential_finder(chemical_symbol,alpha):

    fit = sp.loadtxt('Database/Polyfits/'+chemical_symbol+'polyfit.txt',dtype=float,unpack = True)
    fit = fit[::-1]
    poly = np.poly1d(fit)
    y = poly(sp.log10(alpha))

    return y
    

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

#Define the periodic table
elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
elementList.append('override')
override = 'override'
#==============================MAIN CODE================================#

#Ask the user to input the species of the ions in the plasma or to override and input a custom value of mu.
species = speciesinput()
while (species in elementList) == False:
    print('This species does not exist')
    species = speciesinput()

else:
    if species == 'override':
        counter = 0
        while counter == 0:
            try:
                mu = float(eval(input('Enter the mu value; ').replace('^','**'))) #Mu value
                if mu > 0:
                    print(str(mu).replace('e','*10^'))
                    counter=1
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
    else:
        proton_number = 'pt.' + species + '.number'
        mass = 'pt.' + species + '.mass'
        m_a = eval(mass)
        m_i = (m_a)*u #[kg]
        print(species, str(m_i).replace('e','*10^'), ' kg')
        mu = sp.sqrt(m_i/m_e) #Mu value



Phi = ABR_potential_finder(species,1)
print(Phi)

plt.figure(6, figsize = [8,5])
open_poly(species)
plt.title('Dust surface potential variation with normalised dust grain radius - ABR')
plt.ylabel('Normalised surface potenital')
plt.xlabel('Normalised dust grain radius')
plt.grid()
plt.xscale('log')
plt.legend()

plt.show()
