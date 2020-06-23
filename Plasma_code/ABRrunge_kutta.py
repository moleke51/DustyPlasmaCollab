import scipy as sp
import matplotlib.pyplot as plt

from datetime import datetime as dt
start = dt.now()


def ABR_f(x,y,z):
    return z

def ABR_g(x,y,z,A):
    return -(2/x)*z + (A/(x**2))*(y**(-0.5)) - sp.exp(-y)

def ABR_RK(x0,y0,z0,X,N,A):
    h = (X-x0)/N
    x = X - N*h
    y = y0
    z = z0
    for n in range(N):
        k1 = ABR_f(x,y,z)
        l1 = ABR_g(x,y,z,A)

        k2 = ABR_f(x+h/2,y+k1*h/2,z+l1*h/2)
        l2 = ABR_g(x+h/2,y+k1*h/2,z+l1*h/2,A)

        k3 = ABR_f(x+h/2,y+k2*h/2,z+l2*h/2)
        l3 = ABR_g(x+h/2,y+k2*h/2,z+l2*h/2,A)

        k4 = ABR_f(x+h,y+k3*h,z+l3*h)
        l4 = ABR_g(x+h,y+k3*h,z+l3*h,A)

        x += h
        y += 1/6 * (k1 + 2*k2 + 2*k3 + k4)*h
        z += 1/6 * (l1 + 2*l2 + 2*l3 + l4)*h
                
    return(x,y)

alpha, Phi_a = ABR_RK(7.03,4.10e-4,-2.34e-4,10**(-6),100000,1)
J_ = alpha**2 * 43/(sp.sqrt(4*sp.pi)) * sp.exp(-Phi_a)
print(alpha, Phi_a)
print(J_)
end = dt.now()
print(f'Code completed in {(end-start).seconds}s')