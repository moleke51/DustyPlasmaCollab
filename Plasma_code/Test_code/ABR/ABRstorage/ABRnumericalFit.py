import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from scipy.special import erf
from scipy.optimize import fsolve

u = 1.66e-27
m_e = 9.11e-31


def element_graph(chemical_symbol):
    alpha, Phi_a = sp.loadtxt(chemical_symbol+'data.txt',dtype=float,delimiter=',',unpack=True)
    plt.plot(alpha, Phi_a, label = chemical_symbol)
    
    return(alpha, Phi_a)

def ABR_limit(mu):
    return(-(0.4189 - sp.log(mu)))


def ABR_poly(element,z,order):
    alpha, Phi_a = sp.loadtxt('Database/ABRoutputs/'+element+str(int(z))+'data.txt',dtype=float,delimiter=',',unpack=True)
    Phi_fit = np.poly1d(np.polyfit(sp.log10(alpha), Phi_a, order))
    plt.plot(alpha, Phi_fit(sp.log10(alpha)))
    
    with open('Database/Polyfits/'+element + str(int(z))+ "polyfit.txt","w") as f:
        for i in range(len(Phi_fit)+1):
            f.write(str(Phi_fit[i])+'\n')



elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
elementList.append('override')
override = 'override'

chem_symbols = []
for elem in elementList:
    if elem == 'n':
        continue
    if len(elem) <=2:
        chem_symbols.append(elem)
        if elem == 'Xe':
            break

mu = sp.zeros(len(chem_symbols))
lims = sp.zeros(len(chem_symbols))


for i in range(len(chem_symbols)):
    
    ABR_poly(chem_symbols[i],1,10)
    proton_number = 'pt.' + chem_symbols[i] + '.number'
    mass = 'pt.' + chem_symbols[i] + '.mass'
    m_a = eval(mass)
    m_i = (m_a)*u #[kg]
    mu[i] = sp.sqrt(m_i/m_e) #Mu value    
    lims[i] = ABR_limit(mu[i])


#a,p = element_graph('H')
#loga = sp.log10(a)
#plt.plot(a,p,'x',label = 'Data')
#logfit = np.poly1d(np.polyfit(loga, p, 10))
#plt.plot(a,logfit(loga))
plt.title('ABR')
plt.ylabel('Normalised surface potential')
plt.xlabel('Normalised dust grain radius')
plt.xscale('log')
plt.grid()
plt.legend()
plt.show()