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
#def Limit_potential_finder(nu,mu,C,Theta,m_a): #nu and mu here are arbitrary constants

    #x = nu*sp.log(m_a) + mu*sp.log(Theta) + C

    #return x


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

#Ask the user to input the species of the ions in the plasma or to override and imput a custom value of mu.
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
                m_i = mu*m_e - z*m_e
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
                m_i -= z*m_e
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

#Ask the user for the plasmas ion temperature in Kelvin, for temperatures in electron volts input eV after the value.
#The user may directly input a theta value by typing override. 
counter = 0
thetaOverride = 0
while counter == 0:
    try:
        T_i = input('Set ion temperature (Kelvin); ').replace('^','**').replace('eV',' * e/k_B')
        T_i = float(eval(T_i))
        while (T_i <= 0) == True:
            print('This temperature violates 3rd law of Thermodynamics')
            T_i = float(eval(input('Set ion temperature (Kelvin); ').replace('^','**')))
        else:
            counter = 1
            print(str(T_i).replace('e','*10^'),'K')
    except ValueError:
        if T_i.lower() == 'override':
            counter = 1
            while counter == 1:
                try:
                    #THETA VALUE
                    Theta = float(eval(input('Set the Theta value; ').replace('^','**')))
                    if Theta > 0:
                        counter = 2
                        thetaOverride = 1
                    else:
                        print('Theta must be greater than zero.')
                        
                except ValueError:
                    print('Invalid Theta value.')
                except NameError:
                    print('Invalid Theta value.')
                except TypeError:
                    print('Invalid Theta value.')
                except SyntaxError:
                    print('Invalid Theta value.')
        else:
            print('Invalid ion temperature.')
    except NameError:
        print('Invalid ion temperature.')
    except TypeError:
        print('Invalid ion temperature.')
    except SyntaxError:
        print('Invalid ion temperature.')
#Ask the user for the electron temperature in Kelvin, for temperatures in electron volts input eV after the value.
#The user will then be asked for an electron density in electrons per metre cubed, and a dust grain radius in metres.
#The user may skip this step by inputting override, this will allow an alpha value to be directly inputted.

counter = 0
while counter == 0:
    try:
        T_e = input('Set electron temperature (Kelvin); ').replace('^','**').replace('eV',' * e/k_B')
        T_e = float(eval(T_e))
        if (T_e <= 0) == True:
            print('This temperature violates 3rd law of Thermodynamics')
            
        else:
            print(str(T_e).replace('e','*10^'),'K')
            counter = 1
            if thetaOverride == 0:   
                Theta = T_i/T_e #THETA VALUE
                print('Theta is '+str(Theta).replace('e','*10^'))
            
                #ASK ELECTRON DENSITY
            while counter == 1:
                try:
                    n_e = float(eval(input('Set the electron density [electrons][m^-3] ').replace('^','**')))
                    if (n_e <= 0) == True:
                        print('The electron density must be greater than 0')   
                    else:
                        print(str(n_e).replace('e','*10^'),'electrons m^-3')
                        counter = 2
                        lambda_D = sp.sqrt(epsilon_0*k_B*T_e/(n_e*e**2))
                        print(str(lambda_D).replace('e','*10^'),'m')
                        while counter == 2:
                            try:
                                a = float(eval(input('Set the dust radius [m] ').replace('^','**')))
                                if (a <= 0) == True:
                                    print('The dust radius must be greater than 0')
                                else:
                                    print(str(a).replace('e','*10^'),'m')
                                    counter = 3
                                    alpha = a/lambda_D
                                    print(str(alpha).replace('e','*10^'))
                            except ValueError:
                                print('Invalid dust radius density')
                            except NameError:
                                print('Invalid dust radius density')
                            except TypeError:
                                print('Invalid dust radius density')
                            except SyntaxError:
                                print('Invalid dust radius density')
                except ValueError:
                    print('Invalid electron density')
                except NameError:
                    print('Invalid electron density')
                except TypeError:
                    print('Invalid electron density')
                except SyntaxError:
                    print('Invalid electron density')
    except ValueError:
        if T_e.lower() == 'override':
            counter = 1
            thetaOverride = 1
            while counter == 1:
                try:
                    alpha = float(eval(input('Set the alpha value; ').replace('^','**')))
                    if alpha > 0:
                        counter = 2
                        print(str(alpha).replace('e','*10^'))
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

        else:
            print('Invalid electron temperature.') 
    except NameError:
        print('Invalid electron temperature.')
    except TypeError:
        print('Invalid electron temperature.')
    except SyntaxError:
        print('Invalid electron temperature.')
        
#Ask the user for the flow speed of the plasma in metres per second, for stationary plasmas input 0.
#The user may directly input an upsilon value by inputting override.
counter = 0
if thetaOverride == 0:
    while counter == 0:
        try:
            v = input('Set plasma flow speed (m/s); ').replace('^','**')
            v = sp.absolute(float(eval(v)))
            counter = 1
            upsilon = v / sp.sqrt(k_B*T_e/m_i)
            print(str(upsilon).replace('e','*10^'))
            if v == 0:
                print('This is a static plasma')
            else:
                print('We will use a flowing plasma model, with speed', v ,'m/s')
            
        except ValueError:
            print('Invalid flow speed.')
        except NameError:
            print('Invalid flow speed.')
        except TypeError:
            print('Invalid flow speed.')
        except SyntaxError:
            print('Invalid flow speed.')     
else:
        while counter == 0:
            try:
                upsilon = sp.absolute(float(eval(input('Set the upsilon value; ').replace('^','**'))))
                counter = 1
                print(str(upsilon).replace('e','*10^'))
            except ValueError:
                print('Invalid flow speed.')
            except NameError:
                print('Invalid flow speed.')
            except TypeError:
                print('Invalid flow speed.')
            except SyntaxError:
                print('Invalid flow speed.')  



print('mu =',str(mu).replace('e','*10^'),'z =',str(z).replace('e','*10^'),'Theta =',str(Theta).replace('e','*10^'),'alpha =',str(alpha).replace('e','*10^'),'upsilon =',str(upsilon).replace('e','*10^'))


#Warm of cold test
if Theta >= Theta_critical:
    #Small, medium or large test
    if alpha > 1.25*(Theta**0.4): #Medium or big
        if alpha < 50: #Medium
            #Flowing test
            if upsilon > upsilon_critical:
                print('We need a model for medium flowing plasmas')
                #We are assuming that SOML and SMOML break down at the same points as OML and MOML

                alpha_OML = 1.25*(Theta)**0.4
                alpha_TS = 50
                Phi_SMOML = SMOML_surface_potential_finder_abs(Theta,mu,z,5/3)
                Phi_SOML = SOML_surface_potential_finder_abs(Theta,mu,z)
                Phi_func_of_alph_flow = Linear_function(Phi_SOML,Phi_SMOML,alpha_OML,alpha_TS,alpha)
                #Convert this to a negative value
                Phi = -1 * Phi_func_of_alph_flow
                print(Phi)
                #Un-normalise the potential
                phi = (Phi * k_B * T_e)/(e)
                print(phi)
            else:
                print('Use the numerical fit for static plasmas')

                alpha_OML = 1.25*(Theta)**0.4
                alpha_TS = 50
                Phi_MOML = MOML_surface_potential_finder_abs(Theta,mu,z,5/3)
                Phi_OML = OML_surface_potential_finder_abs(Theta,mu,z)
                Phi_func_of_alph = Linear_function(Phi_OML,Phi_MOML,alpha_OML,alpha_TS,alpha)
                #Convert this to a negative value
                Phi = -1 * Phi_func_of_alph
                print(Phi)
                #Un-normalise the potential
                phi = (Phi * k_B * T_e)/(e)
                print(phi)

        else: #Large
            if upsilon > upsilon_critical:
                print('Use SMOML')     
                Phi= SMOML_surface_potential_finder(Theta,mu,z,5/3,upsilon)
                print(Phi)
               
            else:
                print('Use MOML')
                Phi = MOML_surface_potential_finder(Theta,mu,z,5/3)
                print(Phi)

    else: #Small
        if upsilon > upsilon_critical:
            print('Use SOML')     
            Phi = SOML_surface_potential_finder(Theta,mu,z,upsilon)
            print(Phi)
        else:
            print('Use OML')
            Phi = OML_surface_potential_finder(Theta,mu,z)
            print(Phi)
else:
    #Flowing test
    if upsilon > upsilon_critical:
        print('We need a flowing ABR')     
    else:
        print('Use ABR')

#Un-normalise the potential

if thetaOverride == 0:
    phi = (Phi * k_B * T_e)/(e)
    print(phi)

