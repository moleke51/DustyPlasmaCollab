import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.special import erf
from scipy.optimize import fsolve


gamma = 10000
F = 1
alpha_lim = 0.00000000001
z = 1
m_e = 9.11e-31 #kg
J = sp.logspace(-7,7,100)
alpha = sp.zeros(len(J))
Phi_a = sp.zeros(len(J))
res = sp.zeros(len(J))

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



def plasma_criterion(phi_b):

    x = (4*(phi_b)**(3/2)*((2*phi_b)-3)*((2*phi_b)+1))/((2*phi_b)-1)**3

    return x

##Define the differential equation ABR 9. 
def JDE(H,rho, J, F):
    return( F*H[1], (-2/rho * H[1]) + (J/(sp.sqrt(H[0]) * rho**2)) -sp.exp(-H[0]))

#Define equation ABR 13 for the boundary potential Phi_b.   
def boundary_func(Phi_b, J, gamma):
    
    sol = (4 * (Phi_b**(3/2)) * (2*Phi_b -3) * (2*Phi_b +1))/((2*Phi_b -1)**3) - J/gamma 
    
    return(sol)

def norm_J_current(alpha,Phi,mu):
    j = (alpha**2)*(mu/((4*sp.pi)**0.5))*sp.exp(-Phi)
    return(j)

#Define a function to solve equation ABR (9).
# Guess Phi_a 
def ABR_potential_finder(mu,alpha_lim,z,J,iter,gamma = 10000, resolution = 0.0001):

    m_e = 9.11e-31 #kg
    Nu = mu/sp.sqrt((4*sp.pi))
    F=1
    Phi_b_initial_guess = 0.25
    Phi_b_solution = fsolve(boundary_func, Phi_b_initial_guess, args = (J, gamma))
    Phi_b = Phi_b_solution[0]
    rho_b = (sp.sqrt(J) * sp.exp(Phi_b/2))/(Phi_b**(1/4))
    #print(rho_b)
    dPhi_drho_b = (2*rho_b/J) * (Phi_b**(3/2))/(Phi_b - 1/2) * sp.exp(-Phi_b)
    H0 = [Phi_b,dPhi_drho_b]
    rhos = sp.linspace(rho_b,alpha_lim,10000000)
    Phis = integrate.odeint(JDE,H0,rhos,args=(J,F))
    Phi_sol = Phis[:,0]
    Phi_a = sp.log(((rhos**2)*Nu)/J)
    for i in range (0,len(Phi_a)):
        if sp.absolute((Phi_a[i] - Phi_sol[i])/Phi_a[i]) < resolution[iter]:
            if i != 0:
                turning_point = sp.absolute(((norm_J_current(rhos[i],Phi_sol[i],mu) - J)//1) - ((norm_J_current(rhos[i-1],Phi_sol[i-1],mu) - J)//1))
                if turning_point:
                    alpha = rhos[i]
                    Phi = Phi_sol[i]
                    print(str.format('{0:.4f}', alpha), str.format('{0:.4f}', Phi) )
                    return(alpha,Phi)


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

print(mu,z)

for i in range(0, len(J)):
    if J[i] < 1e-6 :
        res[i] = 1
    elif J[i] < 1e-4:
        res[i] = 0.01
    elif J[i] < 1e-3:
        res[i] = 0.001
    else:
        res[i] = 0.0001
for i in range(0,len(J)):
    alpha[i],Phi_a[i] = ABR_potential_finder(mu,alpha_lim,z,J[i], i ,gamma,res)

with open(species + "data.txt","w") as f:
    for (distance,potential) in zip(alpha,Phi_a):
        f.write(f"{distance},{potential}\n")

plt.plot(alpha, Phi_a, label = 'Hydrogen plasma')
plt.title('ABR')
plt.ylabel('Normalised surface potential')
plt.xlabel('Normalised dust grain radius')
plt.xscale('log')
plt.grid()
plt.legend()
plt.show()


