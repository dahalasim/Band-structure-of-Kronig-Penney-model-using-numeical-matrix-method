#Code to compare numerical vs analytical results
import numpy as np
from cmath import sqrt,sin,cos,cosh,sinh
import matplotlib.pyplot as plt

def g(a,b,v0,E):
    result = None    
    
    w = sqrt(E)
    z = sqrt(v0-E)
    
    zw_term = (z*z - w*w) /(2*z*w)
    trace = zw_term*sin(w*a)*sinh(z*b) + cos(w*a)*cosh(z*b)
    
    if abs(trace) < 1:
        result = np.real_if_close(np.arccos(trace))
    
    return(result)
   
plt.figure(figsize=(8,6))


def band_struc(v0 ,a ,b ):
    energies = np.linspace(0.1,170,100000)
    
    k = np.array([g(a,b,v0,E) for E in energies])
    neg_k = np.array([-k if k is not None else None for k in k])

    plt.plot(k,energies,'r');
    plt.plot(neg_k,energies,'r');
    


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

def bandstruc_second(v0,a,b,energies):
    energy_array = energies
    k = np.array([g(a, b,v0, E) for E in energy_array])
    neg_k = np.array([-k if k is not None else None for k in k])

    plt.scatter(k, energy_array, facecolors='none', edgecolors='black', s=50) 
    plt.scatter(neg_k, energy_array, facecolors='none', edgecolors='black', s=50)
    


band_struc(100,5/6,1/6)
bandstruc_second(100,5/6,1/6,energies)

plt.axis([-np.pi,np.pi, 0,170])
plt.xlabel('Wavevector(k)')
plt.ylabel('Energy')
plt.title('Numerical vs analytical solution(in reduced zone scheme)')
plt.show()
