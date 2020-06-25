import scipy as sp
import matplotlib.pyplot as plt

from datetime import datetime as dt
start = dt.now()


def f(x,y,z):
    return z

def g(x,y,z):
    return -2*z -y

def RK(x0,y0,z0,X,N):
    h = (X-x0)/N
    x = X - N*h
    y = y0
    z = z0
    for n in range(N):
        k1 = f(x,y,z)
        l1 = g(x,y,z)

        k2 = f(x+h/2,y+k1*h/2,z+l1*h/2)
        l2 = g(x+h/2,y+k1*h/2,z+l1*h/2)

        k3 = f(x+h/2,y+k2*h/2,z+l2*h/2)
        l3 = g(x+h/2,y+k2*h/2,z+l2*h/2)

        k4 = f(x+h,y+k3*h,z+l3*h)
        l4 = g(x+h,y+k3*h,z+l3*h)

        x += h
        y += 1/6 * (k1 + 2*k2 + 2*k3 + k4)*h
        z += 1/6 * (l1 + 2*l2 + 2*l3 + l4)*h
                
    return(x,y)

print(RK(0,1,0,10,100))


end = dt.now()
print(f'Code completed in {(end-start).seconds}s')