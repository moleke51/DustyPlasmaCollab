import matplotlib.pyplot as plt
import scipy as sp
alpha = [0.0001,0.0004,0.0017,0.0382,0.2235,1.5005,9.6315,52.0330,246.3514,1099.6963,4807.4252]
Phi_a = [0.0015,0.0063,0.0261,0.3614,0.9640,1.8420,2.6299,3.0730,3.2519,3.3135,3.3335]
def func(x,a,b,c,d):
    tanh = a*sp.tanh((10**((x-b))))+c
    return(tanh)
x = sp.logspace(-4,4,1000)
y = func(x,3.5,1,-0.4,0.5)
#plt.plot(x,y)

with open("Hdata.txt","w") as f:
    for (distance,potential) in zip(alpha,Phi_a):
        f.write(f"{distance},{potential}\n")



J_group = sp.logspace(-7,7,12)
print(J_group)
plt.plot(alpha, Phi_a,label = 'Hydrogen plasma')
plt.title('ABR')
plt.ylabel('Normalised surface potential')
plt.xlabel('Normalised dust grain radius')
plt.xscale('log')
plt.grid()
plt.legend()
plt.show()

