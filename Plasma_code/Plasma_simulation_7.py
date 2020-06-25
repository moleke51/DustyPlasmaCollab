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
# Z: ion charge
# z = Z/e: normalised ion charge
#============================IMPORT PACKAGES============================#
import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.optimize import fsolve
#=============================SET THE PATH==============================#
#Enter the path to the program
path = '/Users/georgedoran/Google Drive/Dusty_Bois/'
#==============================CONSTANTS================================#
#We define some of the universal constants
e = 1.60e-19 #[C] The charge on an electron
epsilon_0 = 8.85e-12 #[F][m^-1] The permitivity of free space
k_B = 1.38e-23 #[m^2][kg][s^-2][K^-1] The Boltzmann constant
m_e = 9.11e-31 #[kg] The mass of an electron
u = 1.66e-27 #[kg] The mass of a nucleon 
Theta_critical = 10**(-5) #This is an arbitrary value
#==============================MESSAGES================================#
#Define the messages that will be returned for each of the models used 
fmessage = 'We will use the {} model'
message_OML = fmessage.format('OML')
message_MOML = fmessage.format('MOML')
message_SOML = fmessage.format('SOML')
message_SMOML = fmessage.format('SMOML')
message_ABR = fmessage.format('ABR')
message_ABR_flow = 'ABR can not account for flow, therefore we will us the flowing OML solutions. '
message_med_static = fmessage.format('numerical fit of OML')
message_med_flow = fmessage.format('numerical fit of flowing OML')
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

#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
def OML_surface_potential_finder(Theta,mu,z): 
    x = (Theta/z) - realLambertW((mu*sp.sqrt(Theta)/z)*sp.exp(Theta/z))
    return x
#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis (absolute values)
def OML_surface_potential_finder_abs(Theta,mu,z): 
    x = sp.absolute(OML_surface_potential_finder(Theta,mu,z))
    return x

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
def MOML_surface_potential_finder(Theta,mu,z,gamma): 
    x = Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2))
    return x
#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis (absolute values)
def MOML_surface_potential_finder_abs(Theta,mu,z,gamma): 
    x = sp.absolute(MOML_surface_potential_finder(Theta,mu,z,gamma))
    return x

#SOML (Shifted OML) model for normalised dust surface potential - eqn 2.138 in Thomas' thesis
def SOML_surface_potential_finder(Theta,mu,z,upsilon): 
    s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
    s_2 = (sp.sqrt(sp.pi)*sps.erf(upsilon))/(2*upsilon)
    x = (Theta*s_1)/(s_2) - realLambertW(((mu*z*sp.sqrt(Theta))/(s_2))*sp.exp((Theta*s_1)/(s_2)))
    return x
#SOML (Shifted OML) model for normalised dust surface potential - eqn 2.138 in Thomas' thesis (absolute values)
def SOML_surface_potential_finder_abs(Theta,mu,z,upsilon): 
    x = sp.absolute(SOML_surface_potential_finder(Theta,mu,z,upsilon))
    return x

#SMOML (Shifted Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis
def SMOML_surface_potential_finder(Theta,mu,z,gamma,upsilon):
    s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
    s_2 = (sp.sqrt(sp.pi)*sps.erf(upsilon))/(2*upsilon)
    x = (Theta*s_1)/(s_2) - realLambertW(((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))/(s_2))*sp.exp((Theta*s_1)/(s_2))) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2))
    return x

#SMOML (Shifted Modified OML) model for normalised dust surface potential - eqn 2.140 in Thomas' thesis (absolute values)
def SMOML_surface_potential_finder_abs(Theta,mu,z,gamma,upsilon):
    x = sp.absolute(SMOML_surface_potential_finder(Theta,mu,z,gamma,upsilon))
    return x

#ABR numerical model for normalised dust grain surface 
def ABR_potential_finder(chemical_symbol,alpha,z):
    fit = sp.loadtxt(path + 'Database/Polyfits/'+chemical_symbol+str(int(z))+'polyfit.txt',dtype=float,unpack = True)
    fit = fit[::-1]
    poly = np.poly1d(fit)
    y = poly(sp.log10(alpha))
    return y
#ABR model for normalised dust surface potential
def open_poly(chemical_symbol,z):
    fit = sp.loadtxt(path + 'Database/Polyfits/'+chemical_symbol+str(int(z))+'polyfit.txt',dtype=float,unpack = True)
    fit = fit[::-1]
    x = sp.logspace(-4,4,10000)
    poly = np.poly1d(fit)
    y = poly(sp.log10(x))
    if len(species) <= 3:
        i = elementList.index(species)
        name = elementList[i+1].capitalize()
    plt.plot(x,y,color = 'purple',label = str(name))

#Define the differential equation ABR 9. 
def JDE(H,rho, J, z):
    return( H[1], (-2/rho * H[1]) + (J/(sp.sqrt(H[0]*z) * rho**2)) -sp.exp(-H[0]))

#Define equation ABR 13 for the boundary potential Phi_b.   
def boundary_func(Phi_b, J, z,gamma = 10000):
    sol = sp.sqrt(z)*((4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3)) - J/gamma 
    return(sol)

#Define the normalised current J from equation ABR 12
def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*sp.pi)**0.5))*sp.exp(-Phi)
    return(j)

#Define a function to solve equation ABR (9).
# Guess Phi_a 
def ABR_potential_solver(mu,alpha_lim,z,J,iter,gamma = 10000, resolution = 0.0001):
    m_e = 9.11e-31 #kg
    Nu = mu/sp.sqrt((4*sp.pi))
    #Make an initial guess for Phi at the boundary of quasi_neutrality validity.
    Phi_b_initial_guess = 0.25
    #Solve numerically for Phi_b
    Phi_b_solution = fsolve(boundary_func, Phi_b_initial_guess, args = (J,z ,gamma))
    Phi_b = Phi_b_solution[0]
    #Calculate rho at the boundary. 
    rho_b = (sp.sqrt(J) * sp.exp(Phi_b/2))/((Phi_b*z)**(1/4))
    #Calculate the first derivative of Phi with respect to rho evaluated at the boundary.
    dPhi_drho_b = ((2*rho_b/J) * (Phi_b**(3/2))/(Phi_b - 1/2) * sp.exp(-Phi_b))*z
    H0 = [Phi_b,dPhi_drho_b]
    rhos = sp.linspace(rho_b,alpha_lim,10000000)
    #Integrate to find Phi from rho_b to alpha_lim.
    Phis = integrate.odeint(JDE,H0,rhos,args=(J,z))
    Phi_sol = Phis[:,0]
    #Caluate Phi_a for each value of rho and the given value of J.
    Phi_a = sp.log(((rhos**2)*Nu)/J)
    for i in range (0,len(Phi_a)):
        #Test if Phi_a is Phi_sol, if it is keep it.
        if sp.absolute((Phi_a[i] - Phi_sol[i])/Phi_a[i]) < resolution[iter]:
            if i != 0:
                turning_point = sp.absolute(((norm_J_current(rhos[i],Phi_sol[i],mu) - J)//1) - ((norm_J_current(rhos[i-1],Phi_sol[i-1],mu) - J)//1))
                if turning_point:
                    alpha = rhos[i]
                    Phi = Phi_sol[i]
                    print(alpha,Phi)
                    return(alpha,Phi)
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
    elif word.lower() == 'override':
        word = 'override'
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
#Define the function to choose and implement the models based on the inputs.
def Potential_Finder(mu,z,alpha,Theta,upsilon,species):
    #Warm of cold test
    if Theta >= Theta_critical:
        #Small, medium or large test
        if alpha > 1.25*(Theta**0.4): #Medium or big
            if alpha < 50: #Medium
                #Flowing test
                if upsilon > 0:
                    #Use the FNF model
                    alpha_OML = 1.25*(Theta)**0.4
                    alpha_TS = 50
                    Phi_SMOML = SMOML_surface_potential_finder_abs(Theta,mu,z,5/3,upsilon)
                    Phi_SOML = SOML_surface_potential_finder_abs(Theta,mu,z,upsilon)
                    Phi_med_flow = Linear_function(Phi_SOML,Phi_SMOML,alpha_OML,alpha_TS,alpha)
                    #Convert this to a negative value
                    Phi = -1 * Phi_med_flow
                    message = message_med_flow
                    return(Phi,message)

                else:
                    #Use the SNF model
                    alpha_OML = 1.25*(Theta)**0.4
                    alpha_TS = 50
                    Phi_MOML = MOML_surface_potential_finder_abs(Theta,mu,z,5/3)
                    Phi_OML = OML_surface_potential_finder_abs(Theta,mu,z)
                    Phi_func_of_alph = Linear_function(Phi_OML,Phi_MOML,alpha_OML,alpha_TS,alpha)
                    #Convert this to a negative value
                    Phi = -1 * Phi_func_of_alph
                    message = message_med_static
                    return(Phi,message)
            else: #Large
                if upsilon > 0:
                    #Use SMOML.
                    Phi= SMOML_surface_potential_finder(Theta,mu,z,5/3,upsilon)
                    message = message_SMOML
                    return(Phi,message)
                
                else:
                    #Use MOML
                    Phi = MOML_surface_potential_finder(Theta,mu,z,5/3)
                    message = message_MOML
                    return(Phi,message)

        else: #Small
            if upsilon > 0:
                #Use SOML
                Phi = SOML_surface_potential_finder(Theta,mu,z,upsilon)
                message = message_SOML
                return(Phi,message)

            else:
                #Use OML
                Phi = OML_surface_potential_finder(Theta,mu,z)
                message = message_OML
                return(Phi,message)
    else:
        #Flowing test
        if upsilon > 0:
            message = message_ABR_flow  

            if alpha > 1.25*(Theta**0.4): #Medium or big
                if alpha < 50: #Medium
                    #Use FNF
                    alpha_OML = 1.25*(Theta)**0.4
                    alpha_TS = 50
                    Phi_SMOML = SMOML_surface_potential_finder_abs(Theta,mu,z,5/3,upsilon)
                    Phi_SOML = SOML_surface_potential_finder_abs(Theta,mu,z,upsilon)
                    Phi_med_flow = Linear_function(Phi_SOML,Phi_SMOML,alpha_OML,alpha_TS,alpha)
                    #Convert this to a negative value
                    Phi = -1 * Phi_med_flow
                    message += message_med_flow
                    return(Phi,message)    
                else: #Large
                    #Use SMOML
                    Phi= SMOML_surface_potential_finder(Theta,mu,z,5/3,upsilon)
                    message += message_SMOML
                    return(Phi,message)
            else: #Small
                #Use SOML
                Phi = SOML_surface_potential_finder(Theta,mu,z,upsilon)
                message += message_SOML
                return(Phi,message)
        
        else:
            #Use ABR
            #Convert the species into its chemical symbol and check for override. 
            if len(species) >=3 and species != 'override':
                index = elementList.index(species) - 1
                species = elementList[index]
            if species == 'override':
                species = 'mu'+str(int(mu))+'z'
            #If the data exists in the database, extract the data and use it to find the value.
            try:    
                Phi = ABR_potential_finder(species,alpha,z)
                Phi = - Phi
                message = message_ABR
                return(Phi,message)
            except OSError:
                #If the data doesn't exist in the database, then produce the data and use it to find the value. 
                print('These plasma conditions are not availible in the database, please wait while we generate the solution, this can take a while.')
                ABR_gamma = 10000
                alpha_lim = 0.00000000001
                J = sp.logspace(-7,7,100)
                ABR_alpha = sp.zeros(len(J))
                ABR_Phi_a = sp.zeros(len(J))
                res = sp.zeros(len(J))

                #Define the resolution values for different J.
                for i in range(0, len(J)):
                    if J[i] < 1e-6 :
                        res[i] = 1
                    elif J[i] < 1e-5:
                        res[i] = 0.1
                    elif J[i] < 1e-4:
                        res[i] = 0.01
                    elif J[i] < 1e-2:
                        res[i] = 0.001
                    else:
                        res[i] = 0.0001
                for i in range(0,len(J)):
                    ABR_alpha[i],ABR_Phi_a[i] = ABR_potential_solver(mu,alpha_lim,z,J[i], i ,ABR_gamma,res)
                #Save the data 
                with open(path + 'Database/ABRoutputs/' + species + str(int(z)) +"data.txt","w") as f:
                    for (distance,potential) in zip(ABR_alpha,ABR_Phi_a):
                        f.write(f"{distance},{potential}\n")
                #Produce the polyfit and save the coefficients.
                Phi_fit = np.poly1d(np.polyfit(sp.log10(ABR_alpha), ABR_Phi_a, 10))
                with open(path + 'Database/Polyfits/'+species + str(int(z))+ "polyfit.txt","w") as f:
                    for i in range(len(Phi_fit)+1):
                        f.write(str(Phi_fit[i])+'\n')
                #Return the potential
                Phi = ABR_potential_finder(species,alpha,z)
                Phi = - Phi
                message = message_ABR
                return(Phi,message)


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
                mu = eval_input(input('Enter the mu value; ')) #Mu value
                if mu > 0:
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
        #Calculate the mu value.
        index = np.floor(elementList.index(species)/2)
        index2 = 0
        if index == 119:
            index = 1
            index2 = 1
        elif index == 120:
            index = 1
            index2 = 2
        
            
        m_a = (elements[index].mass) + index2*(elements[0].mass)#[kg]
        m_i = (m_a)*u #[kg]
        mu = np.sqrt(m_i/m_e) #Mu value
        proton_number = elements[index].number
        #print(index,m_a,proton_number)

    

#Ask the user for the relative charge on the ions in the plasma in integer multiples of the electron charge.
counter=0
if species.lower()=='override':
    while counter == 0:
        try:
            z = input('Enter the relative charge on the plasma ions, [1.60*10^-19 C] ; ') # z is the ion relative charge
            if '.' in z:
                print('The relative charge must be an integer.')
            else:
                z = float(z)
            if z>0:
                counter=1
                m_i = mu*m_e - z*m_e
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
    z_max = proton_number
    while counter==0:
        try:
            z = input('Enter the relative charge on the plasma ions, [1.60*10^-19 C] ; ') # z is the ion relative charge
            if '.' in z:
                print('The relative charge must be an integer.')
            else:
                z = float(z)
            #z = eval_input(input('Enter the relative charge on the plasma ions, [1.60*10^-19 C] ; ')) 
            
            if (z > z_max) == True:
                print('The maximum charge for', species, 'is +' +str(z_max))
            elif (z < 0) == True:
                print('The charge value must be greater than 0.')
                
            elif (z == 0) == True:
                print('Species must be an ion.')
                
            else:
                m_i -= z*m_e
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
        T_i = input('Set ion temperature [K]; ').replace('eV','ev').replace('ev','*11594')
        T_i = eval_input(T_i)
        if (T_i < 0) == True:
            print('The temperature must be greater or equal to zero')
        else:
            counter = 1
    except ValueError:
        if T_i.lower() == 'override':
            counter = 1
            while counter == 1:
                try:
                    #THETA VALUE
                    Theta = eval_input(input('Set the Theta value; '))
                    if Theta >= 0:
                        counter = 2
                        thetaOverride = 1
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
        T_e = input('Set electron temperature [K]; ').replace('eV','ev').replace('ev','*11594')
        T_e = eval_input(T_e)
        if (T_e <= 0) == True:
            print('The temperature must be greater than zero')
        else:
            
            counter = 1
            if thetaOverride == 0:   
                Theta = T_i/T_e #THETA VALUE
            #ASK ELECTRON DENSITY
            while counter == 1:
                try:
                    n_e = eval_input(input('Set the electron density [electrons][m^-3] ').replace('k','*10**3').replace('M','*10**6').replace('G','*10**9').replace('T','*10**12').replace('P','*10**15'))
                    if (n_e <= 0) == True:
                        print('The electron density must be greater than 0')   
                    else:
                        counter = 2
                        lambda_D = sp.sqrt(epsilon_0*k_B*T_e/(n_e*e**2))
                        while counter == 2:
                            try:
                                a = eval_input((input('Set the dust radius [m] ').replace('m','*10**-3').replace('u','*10**-6').replace('n','*10**-9')))
                                if (a <= 0) == True:
                                    print('The dust radius must be greater than 0')
                                else:
                                    counter = 3
                                    alpha = a/lambda_D
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
                    alpha = eval_input(input('Set the alpha value; '))
                    if alpha > 0:
                        counter = 2
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
else:
        while counter == 0:
            try:
                upsilon = sp.absolute(eval_input(input('Set the upsilon value; ')))
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

#Return the normalised potential
Phi, model = Potential_Finder(mu,z,alpha,Theta,upsilon,species)
print(model)
print('The normalised dust grain surface potential is:',Phi)

#Return the potentail and the charge if available.
if thetaOverride == 0:
    phi = (Phi * k_B * T_e)/(e)
    Q = DH_potential_to_charge(a,phi,lambda_D)
    print('The dust grain surface potential is ' +str(phi)+ 'V')
    print('The charge on the dust grain is ' +str(Q)+ 'C')
