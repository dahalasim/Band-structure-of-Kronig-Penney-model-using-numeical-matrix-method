#Code to generate KP potential
import numpy as np
import matplotlib.pyplot as plt

nb=5
v0=100
l=nb
b=1/6
x_r=[-1/2+r for r in range(1,nb+1)]

def potential(x,x_r,b,v0):
    potential=0
    for i in range(nb):
        if x_r[i]-b/2<=x<=x_r[i]+b/2:
            potential=v0
    return potential

x=np.linspace(0,l,1000)
y=[potential(x,x_r,b,v0) for x in x]

plt.figure(figsize=(8,6))
plt.plot(x,y,color='r')
plt.xlabel('x',fontsize=18)
plt.ylabel('Potential',fontsize=15)
plt.title('Kronig-Penney Model', fontsize=15)
plt.grid()
plt.show()





 