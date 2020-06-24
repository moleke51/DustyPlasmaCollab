import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
from scipy import integrate
#erf
upsilon = sp.linspace(0,10, 10000)
s_1 = ((sp.sqrt(sp.pi))*(1+2*(upsilon**2))*sps.erf(upsilon))/(4*upsilon) + 0.5*sp.exp(-(upsilon**2))
s_2 = (sp.sqrt(sp.pi)*sps.erf(upsilon))/(2*upsilon)

#y = sps.erf(x)
#yprime = sp.exp(-x**2)
plt.plot(upsilon,s_1,color='red',label = 'S1')
plt.plot(upsilon,s_2, color='darkorange',label='S2')
plt.plot(upsilon, s_1/s_2, color = 'green', label = 'S1/S2')
plt.legend()
plt.grid()
plt.show()

