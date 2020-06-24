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
# upsilon = v / [2*(k_B * T_e)/m_i]^1/2: normalised flow speed
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
def open_poly(chemical_symbol,z):
    fit = sp.loadtxt('Database/Polyfits/'+chemical_symbol+str(int(z))+'polyfit.txt',dtype=float,unpack = True)
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


    

#Ask the user for the charge on the ions in the plasma in Coulombs, charge can be added in multiples of e by inputting e after the value.
counter=0
if species.lower()=='override':
    while counter == 0:
        try:
            Z = float(eval((input('Enter the charge on the plasma ions, [C] ; ').replace('^','**')).replace('e',' * e'))) # Z is the ion charge
            if Z>0:
                print(str(Z).replace('e','*10^'),'C')
                counter=1
                z = Z/e
                print(z)
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
else:
    Z_max = eval(proton_number) * e
    while counter==0:
        try:
            Z = float(eval((input('Enter the charge on the plasma ions, [C] ; ').replace('^','**')).replace('e',' * e'))) # Z is the ion charge
            
            if (Z > Z_max) == True:
                print('The maximum charge for', species, 'is +' +str(Z_max).replace('e','*10^'),'C')
                
            elif (Z < 0) == True:
                print('The charge value must be greater than 0')
                
            elif (Z == 0) == True:
                print('Species must be an ion')
                
            else:
                print(str(Z).replace('e','*10^'),'C')
                z = Z/e
                print(z)
                counter = 1
        except ValueError:
            print('Invalid charge value.')
        except NameError:
            print('Invalid charge value.')
        except TypeError:
            print('Invalid charge value.')
        except SyntaxError:
            print('Invalid charge value.')


Theta = input('Set the Theta value (Ion temperature/Electron temperature)')
Theta = float(eval(Theta))



#BOUNDARY CONDITIONS

alpha_OML = 1.25*(Theta)**0.4
alpha_TS = 50
#OML used in this range
alpha_1 = sp.logspace(-3,sp.log10(alpha_OML),10)
Theta_OML = sp.array([Theta] * len(alpha_1))
#Linear fit model used in this range
alpha_2 = sp.logspace(sp.log10(alpha_OML), sp.log10(50),10)
#MOML used in this range
alpha_3 = sp.logspace(sp.log10(50),2,10)
Theta_MOML = sp.array([Theta] * len(alpha_3))

#Potential variation with distance from centre of dust grain and potential variation with distance from centre of dust grain for a flowing plasma

upsilon = sp.arange(1,6,1)
colour = ['red','green','darkorange','cyan','purple','pink']

#Non-flowing plasma
Phi_MOML = MOML_surface_potential_finder_abs(Theta_MOML,mu,z,5/3)
phi_MOML = Phi_MOML[1]
Phi_OML = OML_surface_potential_finder_abs(Theta_OML,mu,z)
phi_OML = Phi_OML[1]
Phi_func_of_alph = Linear_function(phi_OML,phi_MOML,alpha_OML,alpha_TS,alpha_2)


plt.figure(1, figsize= [12,5])
for i in range(0,len(upsilon)):
    Phi_SMOML_u = SMOML_surface_potential_finder_abs(Theta_MOML,mu,z,5/3,upsilon[i])
    phi_SMOML_u = Phi_SMOML_u[1]
    Phi_SOML_u = SOML_surface_potential_finder_abs(Theta_OML,mu,z,upsilon[i])
    phi_SOML_u = Phi_SOML_u[1]
    Phi_func_of_alph_flow = Linear_function(phi_SOML_u,phi_SMOML_u,alpha_OML,alpha_TS,alpha_2)
    plt.plot(alpha_1,Phi_SOML_u, color = colour[i])
    plt.plot(alpha_2, Phi_func_of_alph_flow , label='Flowing plasma, u = ' + str(upsilon[i]) , color = colour[i] )
    plt.plot(alpha_3, Phi_SMOML_u, color = colour[i])
    
 
plt.plot(alpha_1,Phi_OML,color = 'Blue')
plt.plot(alpha_2, Phi_func_of_alph , label='Static plasma' , color = 'Blue' )
plt.plot(alpha_3, Phi_MOML,color = 'Blue')
plt.title('LINEAR FIT model of Dust surface potential variation with normalised radius for flowing and static plasma')
plt.ylabel('Normalised surface potential, $\Phi_a$')
plt.xlabel('Normalsied radius, $\\alpha $ (Log scale)')
plt.xscale('log')
plt.grid()
plt.legend()

#OML FUNCTION - VARIATION WITH THETA

Theta__OML = sp.logspace(-2,1,100)
Phi__OML1 = OML_surface_potential_finder(Theta__OML,mu,z)
Phi__OML2 = OML_surface_potential_finder(Theta__OML,112,z)
Phi__OML3 = OML_surface_potential_finder(Theta__OML,204,z)

if len(species) <= 3:
    i = elementList.index(species)
    name = elementList[i+1].capitalize()

plt.figure(2, figsize= [8,5])
plt.plot(Theta__OML, Phi__OML1, label = str(name))
plt.plot(Theta__OML, Phi__OML2, label = 'Lithium')
plt.plot(Theta__OML, Phi__OML3, label = 'Sodium')
plt.title('Dust surface potential variation with Theta - OML')
plt.ylabel('Normalised surface potential, $\Phi_a$')
plt.xlabel('$\Theta$ (Log scale)')
plt.legend()
plt.xscale('log')
plt.grid()


#MOML FUNCTION - VARIATION WITH THETA

Theta__MOML = sp.logspace(-2,2,100)
Phi__MOML1 = MOML_surface_potential_finder(Theta__MOML,mu,z,5/3) #GAMMA = 1
Phi__MOML2 = MOML_surface_potential_finder(Theta__MOML,112,z,5/3) #GAMMA = 5/3
Phi__MOML3 = MOML_surface_potential_finder(Theta__MOML,204,z,5/3) #GAMMA = 3

plt.figure(3, figsize= [8,5])
plt.plot(Theta__MOML[1:], Phi__MOML1[1:], label= str(name))
plt.plot(Theta__MOML[1:], Phi__MOML2[1:], label='Lithium')
plt.plot(Theta__MOML[1:], Phi__MOML3[1:], label='Sodium')
plt.title('Dust surface potential variation with Theta - MOML')
plt.ylabel('Normalised surface potential, $\Phi_a$')
plt.xlabel('$\Theta$ (Log scale)')
plt.xscale('log')
plt.grid()
plt.legend()

Phi__MOML1 = MOML_surface_potential_finder(Theta__MOML,mu,z,1) #GAMMA = 1
Phi__MOML2 = MOML_surface_potential_finder(Theta__MOML,mu,z,5/3) #GAMMA = 5/3
Phi__MOML3 = MOML_surface_potential_finder(Theta__MOML,mu,z,3) #GAMMA = 3

plt.figure(4, figsize= [8,5])
plt.plot(Theta__MOML[1:], Phi__MOML1[1:], label= 'gamma = 1')
plt.plot(Theta__MOML[1:], Phi__MOML2[1:], label= 'gamma = 5/3')
plt.plot(Theta__MOML[1:], Phi__MOML3[1:], label= 'gamma = 3')
plt.title('Dust surface potential variation with Theta - MOML')
plt.ylabel('Normalised surface potential, $\Phi_a$')
plt.xlabel('$\Theta$ (Log scale)')
plt.xscale('log')
plt.grid()
plt.legend()

#SOML FUNCTION - VARIATION WITH THETA

Theta__SOML= sp.logspace(-2,1,100)
Phi__SOML1 = SOML_surface_potential_finder(Theta__SOML,43,z,1) #NORMALISED PLASMA SPEED = 1
Phi__SOML2 = SOML_surface_potential_finder(Theta__SOML,43,z,2) #NORMALISED PLASMA SPEED = 2
Phi__SOML4 = SOML_surface_potential_finder(Theta__SOML,43,z,4) #NORMALISED PLASMA SPEED = 4

plt.figure(5, figsize= [8,5])
plt.plot(Theta__SOML, Phi__SOML1,label = 'u = 1') 
plt.plot(Theta__SOML, Phi__SOML2,label = 'u = 2')
plt.plot(Theta__SOML, Phi__SOML4,label = 'u = 4')
plt.title('Dust surface potential variation with Theta - SOML')
plt.ylabel('Normalised surface potential, $\Phi_a$')
plt.xlabel('$\Theta$ (Log scale)')
plt.xscale('log')
plt.grid()
plt.legend()

#SMOML FUNCTION - VARIATION WITH THETA


Phi__SMOML1 = SMOML_surface_potential_finder(Theta__SOML,43,z,5/3,1) #NORMALISED PLASMA SPEED = 1
Phi__SMOML2 = SMOML_surface_potential_finder(Theta__SOML,43,z,5/3,2) #NORMALISED PLASMA SPEED = 2
Phi__SMOML4 = SMOML_surface_potential_finder(Theta__SOML,43,z,5/3,4) #NORMALISED PLASMA SPEED = 4

plt.figure(6, figsize= [8,5])
plt.plot(Theta__SOML, Phi__SMOML1, label = r'$\upsilon$ = '+'1')
plt.plot(Theta__SOML, Phi__SMOML2, label = r'$\upsilon$ = '+'2')
plt.plot(Theta__SOML, Phi__SMOML4, label = r'$\upsilon$ = '+'4')
plt.title('Dust surface potential variation with Theta - SMOML')
plt.ylabel('Normalised surface potential, $\Phi_a$')
plt.xlabel('$\Theta$ (Log scale)')
plt.xscale('log')
plt.grid()
plt.legend()

#ABR FUNCTION - VARIATION WITH ALPHA

plt.figure(7, figsize = [8,5])
open_poly(species,z)
plt.title('Dust surface potential variation with normalised dust grain radius - ABR')
plt.ylabel('Normalised surface potenital, $\Phi_a$')
plt.xlabel('Normalised dust grain radius, $\\alpha$')
plt.grid()
plt.xscale('log')
plt.legend()

plt.show()
