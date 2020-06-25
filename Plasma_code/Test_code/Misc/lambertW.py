import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt 


def realLambertW(x):
    if type(x) is float or type(x) is int or type(x) is np.float64:
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




x = sp.arange(-0.367879,10,0.01)
y = realLambertW(x)
plt.plot(x,y)
plt.grid()
plt.show()