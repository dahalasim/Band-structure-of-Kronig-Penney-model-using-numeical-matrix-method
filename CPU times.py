#Code to calculate CPU time
import numpy as np
import matplotlib.pyplot as plt
import time

nb=10 
l=nb
v0=100
N=200
x_r=[-1/2+r for r in range(1,nb+1)]
b=1/6

def f(k,x,l):
    return np.sin(k*np.pi*x/l)/(k*np.pi)

def Fnn(x,l,n):
    return x/l-f(2*n,x,l) 

def Fmn(x,l,m,n):
    return f(m-n,x,l)-f(m+n,x,l)

def hnn(s,b,l,n):
    return Fnn(s+b/2,l,n)-Fnn(s-b/2,l,n)

def hmn(s,b,l,m,n):
    return Fmn(s+b/2,l,m,n)-Fmn(s-b/2,l,m,n)

def Hnn(n,l,x_r,v0,b):
    result=(n*np.pi/l)**2
    for i in range(len(x_r)):
        result+=v0*hnn(x_r[i],b,l,n)
    return result

def Hmn(m,n,l,x_r,v0,b):
    result=0
    for i in range(len(x_r)):
        result+=v0*hmn(x_r[i],b,l,m,n)
    return result

def Hamiltonian(N,l,x_r,v0,b):
    matrix=[]
    for m in range(1,N+1):
        row=[]
        for n in range(1,N+1):
            if m==n:
                row.append(Hnn(n,l,x_r,v0,b))
            else:
                row.append(Hmn(m,n,l,x_r,v0,b))
        matrix.append(row)
    return matrix

time_matrix=time.time()
Hamiltonian_matrix=Hamiltonian(N,l,x_r,v0,b)
time_matrix=time.time()-time_matrix

time_eigensystem=time.time()
u,v=np.linalg.eig(Hamiltonian_matrix)
time_eigensystem=time.time()-time_eigensystem

total_time=time_matrix+time_eigensystem

print(f'Time to build Hamiltonian matrix = {round(time_matrix,3)} seconds')
print(f'Time to solve the eigensytem = {round(time_eigensystem,3)} seconds')
print(f'Total time = {round(total_time,3)} seconds\n\n')
