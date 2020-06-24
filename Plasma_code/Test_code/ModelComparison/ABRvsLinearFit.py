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

Path = '/Users/doganakpinar/Google Drive /Dusty_Bois/'
#Path = '/Users/georgedoran/Google Drive/Dusty_Bois/ModelComparison/'
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


#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis (absolute values)
def MOML_surface_potential_finder_abs(Theta,mu,z,gamma): #we are ingnoring absorption radii outside the dust grain

    x = sp.absolute(Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2)))

    return x


#OML model for normalised dust surface potential - eqn 2.107 in Thomas' thesis
def OML_surface_potential_finder(Theta,mu,z): #we are ingnoring absorption radii outside the dust grain

    x = (Theta/z) - realLambertW((mu*sp.sqrt(Theta)/z)*sp.exp(Theta/z))

    return x

#MOML (Modified OML) model for normalised dust surface potential - eqn 2.130 in Thomas' thesis
def MOML_surface_potential_finder(Theta,mu,z,gamma): #we are ingnoring absorption radii outside the dust grain

    x = Theta - realLambertW((sp.sqrt(2*sp.pi*Theta*(1+gamma*Theta)))*sp.exp(Theta)) + 0.5*sp.log((2*sp.pi*(1+gamma*Theta))/((z**2)*(mu)**2))

    return x


#ABR model for normalised dust surface potential
def open_poly(chemical_symbol,z,alpha):
    fit = sp.loadtxt(Path + 'DustInPlasmaTestFunctions/Database/Polyfits/'+chemical_symbol+str(int(z))+'polyfit.txt',dtype=float,unpack = True)
    fit = fit[::-1]
    poly = np.poly1d(fit)
    y = poly(sp.log10(alpha))
    
    return(y)
    

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


def static_OML_Finder(mu,z,alpha,Theta):
    #Warm of cold test
    
        #Small, medium or large test
    if alpha > 1.25*(Theta**0.4): #Medium or big
        if alpha < 50: #Medium   
            alpha_OML = 1.25*(Theta)**0.4
            alpha_TS = 50
            Phi_MOML = MOML_surface_potential_finder_abs(Theta,mu,z,5/3) +0.5
            Phi_OML = OML_surface_potential_finder_abs(Theta,mu,z)
            Phi_func_of_alpha = Linear_function(Phi_OML,Phi_MOML,alpha_OML,alpha_TS,alpha)
            #Convert this to a negative value
            Phi = -1 * Phi_func_of_alpha
            return(Phi)
        else: #Large  
            Phi = MOML_surface_potential_finder(Theta,mu,z,5/3) -0.5
            return(Phi)

    else: #Small
        Phi = OML_surface_potential_finder(Theta,mu,z)
        return(Phi)

            
        
        

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


#BOUNDARY CONDITIONS
Theta = sp.logspace(-8,-4,5)
#Theta = sp.logspace(-10,-4,7)
print(Theta)
colour = ['purple','blue','green','darkorange','red','pink']


alpha = sp.logspace(-4,4,10000)
plt.figure(1, figsize = [8,5])
ABR = open_poly(species,z,alpha)
#Linear fit function for a static plasma
for i in range(len(Theta)):
    Phi = sp.zeros(len(alpha))
    y = sp.linspace(0,0.7,len(alpha))
    ALPHA = 1.25 *(Theta[i]**0.4)
    x = sp.ones(len(alpha))*ALPHA
    plt.plot(x,y,'--',color=colour[i])
    for j in range(len(alpha)):
        Phi[j] = -static_OML_Finder(mu,z,alpha[j],Theta[i])
    #plt.plot(alpha,Phi,label = '$\Theta$ = ' + str(Theta[i]),color = colour[i])
    res = sp.zeros(len(alpha))
    for k in range(len(alpha)):
        res[k] = sp.sqrt((Phi[k]-ABR[k])**2)
    plt.plot(alpha,res,label = '$\Theta$ = ' + str(Theta[i]).replace('e','x10^'),color = colour[i])
x = 50*sp.ones(len(alpha))
plt.plot(x,y,'--',color = 'black')  

plt.title('The residuals between ABR and the OML family against normalised dust radius')
plt.ylabel('The difference in the Normalised surface potenital predictions, $\Delta\Phi_a$ ')
plt.xlabel('Normalised dust grain radius, $\\alpha$')
plt.grid()
plt.xscale('log')
plt.legend(loc = 'upper right')


plt.figure(2, figsize = [8,5])
for i in range(len(Theta)):
    Phi = sp.zeros(len(alpha))
    
    for j in range(len(alpha)):
        Phi[j] = -static_OML_Finder(mu,z,alpha[j],Theta[i])
    plt.plot(alpha,Phi,label = 'OML family $\Theta$ = ' + str(Theta[i]).replace('e','x10^'),color = colour[i])

plt.plot(alpha,ABR,label = 'ABR',color = 'black')



plt.title('ABR and OML family linear fit normalised potential against normalised dust radius')
plt.ylabel('Normalised surface potenital, $\Phi_a$ ')
plt.xlabel('Normalised dust grain radius, $\\alpha$')
plt.grid()
plt.xscale('log')
plt.legend(loc = 'lower right')


plt.show()