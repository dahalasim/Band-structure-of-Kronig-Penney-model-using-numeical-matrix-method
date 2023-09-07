#Code to generate band structure
import numpy as np
import matplotlib.pyplot as plt

nb=10
l=nb
v0=100
b=1/6
x_r=[-1/2+r for r in range(1,nb+1)]
N=100

def f(k,x,l):
    return np.sin(k*np.pi*x/l)/(np.pi*k)

def Fnn(n,x,l):
    return x/l-f(2*n,x,l)

def Fmn(m,n,x,l):
    return f(m-n,x,l)-f(m+n,x,l)

def hnn(n,s,b,l):
    return Fnn(n,s+b/2,l)-Fnn(n,s-b/2,l)

def hmn(m,n,s,b,l):
    return Fmn(m,n,s+b/2,l)-Fmn(m,n,s-b/2,l)

def Hnn(n,l,v0,x_r,b):
    result=(n*np.pi/l)**2
    for i in range(len(x_r)):
        result+=v0*hnn(n,x_r[i],b,l)
    return result

def Hmn(m,n,v0,x_r,b):
    result=0
    for i in range(len(x_r)):
        result+=v0*hmn(m,n,x_r[i],b,l)
    return result

def Hamiltonian(N,l,v0,x_r,b):
    matrix=[]
    for m in range(1,N+1):
        row=[]
        for n in range(1,N+1):
            if m==n:
                row.append(Hnn(n,l,v0,x_r,b))
            else:
                row.append(Hmn(m,n,v0,x_r,b))
        matrix.append(row)
    return matrix

Hamiltonian_matrix=Hamiltonian(N,l,v0,x_r,b)
u,v=np.linalg.eig(Hamiltonian_matrix)
sorted_indices=np.argsort(u)
energies=u[sorted_indices]

fig,ax=plt.subplots(figsize=(8,6))

for n in range(nb*4+1):
    ax.scatter((n+1)*np.pi/l,energies[n],color='r')

ax.set_xlabel('Wave vector',fontsize=15)
ax.set_ylabel('Energy',fontsize=15)
ax.set_xlim([0,4*np.pi+0.5])
ax.set_ylim([0,170])
plt.show()























