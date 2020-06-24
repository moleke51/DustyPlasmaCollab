import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import periodictable as pt
u = 1.66e-27
m_e = 9.11e-31
plt.figure(1)
def element_graph(chemical_symbol):
    alpha, Phi_a = sp.loadtxt(chemical_symbol+'1data.txt',dtype=float,delimiter=',',unpack=True)
    plt.plot(alpha, Phi_a, label = chemical_symbol)

def ABR_limit(mu):
    return(-(0.4189 - sp.log(mu)))

def open_poly(chemical_symbol,colour,labels):
    fit = sp.loadtxt('/Users/georgedoran/Google Drive/Dusty_Bois/Database/Polyfits/'+chemical_symbol+'1polyfit.txt',dtype=float,unpack = True)
    fit = fit[::-1]
    x = sp.logspace(-4,3.5,10000)
    poly = np.poly1d(fit)
    y = poly(sp.log10(x))
    plt.plot(x,y,color = colour,label = labels)


plt.figure(num=1,figsize=[8,5])

elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())
elementList.append('override')
override = 'override'
chem_symbols = []
labels = []
for elem in elementList:
    if elem == 'H' or elem == 'Ar' or elem == 'Xe':
        chem_symbols.append(elem)
        q = elementList.index(elem)
        labels.append(elementList[q+1].capitalize())

    else:
        continue

    if elem == 'n':
        continue
    if len(elem) <=2:
        #chem_symbols.append(elem)
        if elem == 'Xe':
            break

mu = sp.zeros(len(chem_symbols))
lims = sp.zeros(len(chem_symbols))
colours = ['blue','orange','green']
for i in range(len(chem_symbols)):
    
    open_poly(chem_symbols[i],colours[i],labels[i])
    proton_number = 'pt.' + chem_symbols[i] + '.number'
    mass = 'pt.' + chem_symbols[i] + '.mass'
    m_a = eval(mass)
    m_i = (m_a)*u #[kg]
    
    mu[i] = sp.sqrt(m_i/m_e) #Mu value    
    lims[i] = ABR_limit(mu[i])
x = sp.logspace(-4,3.5,10000)
y = sp.ones(len(x))

for i in range(0,len(mu)):
    plt.plot(x,lims[i]*y,'--',color = colours[i])

plt.title('Dust surface potential variation with normalised dust grain radius - ABR')
plt.ylabel('Normalised surface potential, $\\Phi_a$ ')
plt.xlabel('Normalised dust grain radius, $\\alpha$')
plt.xscale('log')
plt.grid()
plt.legend(loc='lower right',ncol =1)
plt.show()

