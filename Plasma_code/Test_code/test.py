import scipy as sp
import matplotlib.pyplot as plt
import periodictable as pt

x = sp.logspace(-20,-5,1000)
y = x**0.4

plt.plot(x,y)
plt.grid()
plt.xscale('log')
plt.show()